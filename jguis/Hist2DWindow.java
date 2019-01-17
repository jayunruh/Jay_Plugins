/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import jalgs.algutils;
import jalgs.jstatistics;
import jalgs.jfit.linleastsquares;

import java.awt.Button;
import java.awt.Checkbox;
import java.awt.CheckboxGroup;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Hist2DWindow extends Panel implements ActionListener,AdjustmentListener,MouseListener,MouseMotionListener,ItemListener{
	public final static int H=612;
	public final static int WR=800;
	public Dimension totalSize=new Dimension();
	private TextField scalehistval,dispminval,dispmaxval,threshval;
	private TextField xminval,xmaxval,yminval,ymaxval;
	private TextField roiwidthval,roiheightval,roixval,roiyval;
	private TextField sval,s0val,offval,backval,sliceval;
	private Label dispminlabel,dispmaxlabel,scalehistlabel,scorelabel,xlabel,ylabel,slicelabel,xavglabel,yavglabel,threshlabel;
	private Label xstdevlabel,ystdevlabel,zstdevlabel,zavglabel;
	private Label roiwidthlabel,roiheightlabel,roixlabel,roiylabel,slabel,s0label,offlabel,backlabel;
	private Button smooth_button,revert_button,savehist_button,saveimg_button,saveyimg_button,unmix_button,unmix_button2;
	CheckboxGroup roishapegroup;
	Checkbox roisquare,roioval;
	CheckboxGroup plottypegroup;
	Checkbox plotnorm,plotdiff;
	CheckboxGroup smoothtypegroup;
	Checkbox medsmooth,gassmooth,binsmooth;
	Checkbox ascalecheck;
	Checkbox threshfirstcheck;
	Calibration cal;
	Scrollbar slice_slider;
	Image dispimg,histimg,tempimg;
	Graphics tempgraphics;
	float xmin,xmax,ymin,ymax,dispmin,dispmax,multiplier,s,off,s0,score,back;
	float xavg,yavg,thresh,xstdev,ystdev,zavg,zstdev;
	int width,height,roishape,roix,roiy,roiwidth,roiheight,slices,currslice;
	int oldxpos,oldypos,histmax,plottype,pct,smoothtype,ascale;
	float[] xpix,ypix,zpix,histxpix,histypix,histzpix,disppix;
	float[] smoothxpix,smoothypix,smoothzpix;
	int[] mask;
	boolean inroi,threshfirst;
	int imageheight=256;
	public ImageStack datastack;

	public static void launch_frame(String title,Hist2DWindow panel){
		final Frame f=new Frame(title);
		f.setLocation(10,10);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				Component[] comps=f.getComponents();
				for(int i=0;i<comps.length;i++){
					comps[i].setVisible(false);
				}
				f.dispose();
			}
		});

		f.setLayout(null);
		Insets ins=f.getInsets();
		panel.totalSize.height=Hist2DWindow.H+ins.bottom+ins.top+100;
		panel.totalSize.width=Hist2DWindow.WR+ins.left+ins.right;
		panel.setBounds(ins.left+5,ins.top+25,panel.totalSize.width,panel.totalSize.height-50);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(panel.totalSize);
		f.setVisible(true);
		panel.requestFocus();
		// WindowManager.addWindow(f);
	}

	public void init(ImagePlus ximp,ImagePlus yimp,ImagePlus dispimp,int calltype){
		init(ximp,yimp,null,dispimp,calltype);
	}
	
	public void init(ImagePlus ximp,ImagePlus yimp,ImagePlus zimp,ImagePlus dispimp,int calltype){
		init(ximp,yimp,zimp,null,dispimp,calltype);
	}

	public void init(ImagePlus ximp,ImagePlus yimp,ImagePlus zimp,ImagePlus dataimp,ImagePlus dispimp,int calltype){
		//calltype is 1: N&B, 2: acceptor photobleaching fret, and 3: ratiometric fret
		//want to add a spectral phasor option (4)
		//need a placeholder for profile data
		init_options();
		if(calltype==4) {
			ascale=0;
			xmin=-1.0f;
			xmax=1.0f;
			ymin=-1.0f;
			ymax=1.0f;
		}
		setLayout(null);
		// initialize all of the variables
		width=ximp.getWidth();
		int newwidth=width;
		if(newwidth>256){
			newwidth=256;
		}
		height=ximp.getHeight();
		ImageStack xstack=ximp.getStack();
		ImageStack ystack=yimp.getStack();
		ImageStack dispstack=dispimp.getStack();
		ImageStack zstack=null;
		if(zimp!=null)
			zstack=zimp.getStack();
		if(dataimp!=null) datastack=dataimp.getStack();
		cal=dispimp.getCalibration();
		slices=xstack.getSize();
		xpix=new float[width*height*slices];
		ypix=new float[width*height*slices];
		disppix=new float[width*height*slices];
		zpix=null;
		if(zstack!=null)
			zpix=new float[width*height*zstack.getSize()];
		for(int i=0;i<slices;i++){
			System.arraycopy((float[])xstack.getPixels(i+1),0,xpix,i*width*height,width*height);
			System.arraycopy((float[])ystack.getPixels(i+1),0,ypix,i*width*height,width*height);
			System.arraycopy((float[])dispstack.getPixels(i+1),0,disppix,i*width*height,width*height);
		}
		if(zstack!=null){
			for(int i=0;i<zstack.getSize();i++){
				System.arraycopy((float[])zstack.getPixels(i+1),0,zpix,i*width*height,width*height);
			}
		}
		revert();
		dispmin=getminlesszero(disppix);
		dispmax=getmaxgreatzero(disppix);
		currslice=0;
		addMouseMotionListener(this);
		addMouseListener(this);

		xminval=new TextField(""+xmin,10);
		xminval.setBounds(90,10+imageheight+50+256+20,80,20);
		xminval.addActionListener(this);
		xmaxval=new TextField(""+xmax,10);
		xmaxval.setBounds(90+245,10+imageheight+50+256+20,80,20);
		xmaxval.addActionListener(this);
		yminval=new TextField(""+ymin,10);
		yminval.addActionListener(this);
		yminval.setBounds(10,10+imageheight+40+256,80,20);
		ymaxval=new TextField(""+ymax,10);
		ymaxval.setBounds(10,10+imageheight+40,80,20);
		ymaxval.addActionListener(this);
		smooth_button=new Button("Smooth");
		smooth_button.setBounds(10,10+imageheight+70,60,30);
		smooth_button.addActionListener(this);
		revert_button=new Button("Revert");
		revert_button.setBounds(10,10+imageheight+110,60,30);
		revert_button.addActionListener(this);
		savehist_button=new Button("Save");
		savehist_button.setBounds(10,10+imageheight+150,60,30);
		savehist_button.addActionListener(this);
		saveimg_button=new Button("Save Image");
		saveimg_button.setBounds(10,50,80,30);
		saveimg_button.addActionListener(this);
		String saveyimglabel="Save Y Image";
		if(calltype==4) saveyimglabel="Save Profile";
		saveyimg_button=new Button(saveyimglabel);
		saveyimg_button.setBounds(100+256+10+110+10+100,10+imageheight+50+260,80,20);
		saveyimg_button.addActionListener(this);
		unmix_button=new Button("Unmix Histogram");
		unmix_button.setBounds(100+256+10+110,10+imageheight+50+260,100,20);
		unmix_button.addActionListener(this);
		unmix_button2=new Button("Unmix Histogram2");
		unmix_button2.setBounds(100+256+10+110,10+imageheight+50+290,100,20);
		unmix_button2.addActionListener(this);

		roiwidthlabel=new Label("Roi Width");
		roiwidthlabel.setBounds(100+256+10,10+imageheight+50,100,20);
		roiwidthval=new TextField(""+roiwidth,15);
		roiwidthval.setBounds(100+256+10,10+imageheight+50+30,80,20);
		roiwidthval.addActionListener(this);
		roiheightlabel=new Label("Roi Height");
		roiheightlabel.setBounds(100+256+10,10+imageheight+50+60,100,20);
		roiheightval=new TextField(""+roiheight,15);
		roiheightval.setBounds(100+256+10,10+imageheight+50+90,80,20);
		roiheightval.addActionListener(this);
		roixlabel=new Label("Roi X");
		roixlabel.setBounds(100+256+10+110,10+imageheight+50,100,20);
		roixval=new TextField(""+(roix+roiwidth/2),15);
		roixval.setBounds(100+256+10+110,10+imageheight+50+30,80,20);
		roixval.addActionListener(this);
		roiylabel=new Label("Roi Y");
		roiylabel.setBounds(100+256+10+110,10+imageheight+50+60,100,20);
		roiyval=new TextField(""+(roiy+roiheight/2),15);
		roiyval.setBounds(100+256+10+110,10+imageheight+50+90,80,20);
		roiyval.addActionListener(this);

		if(calltype>1){
			s=1.0f;
			s0=0.0f;
			off=0.0f;
			back=0.0f;
		}
		slabel=new Label("Slope");
		slabel.setBounds(100+256+10+110,10+imageheight+50+120,100,20);
		sval=new TextField(""+s,15);
		sval.setBounds(100+256+10+110,10+imageheight+50+150,80,20);
		sval.addActionListener(this);
		s0label=new Label("Zero Var");
		s0label.setBounds(100+256+10+110,10+imageheight+50+180,100,20);
		s0val=new TextField(""+s0,15);
		s0val.setBounds(100+256+10+110,10+imageheight+50+210,80,20);
		s0val.addActionListener(this);
		offlabel=new Label("Offset");
		offlabel.setBounds(100+256+10+110,10+imageheight+50+240,100,20);
		offval=new TextField(""+off,15);
		offval.setBounds(100+256+10+110,10+imageheight+50+270,80,20);
		offval.addActionListener(this);
		backlabel=new Label("Background");
		backlabel.setBounds(100+256+10+110,10+imageheight+50+300,100,20);
		backval=new TextField(""+back,15);
		backval.setBounds(100+256+10+110,10+imageheight+50+330,80,20);
		backval.addActionListener(this);

		score=0.0f;
		scorelabel=new Label("Score = "+score);
		scorelabel.setBounds(100+256+10+110+10+100,10+imageheight+50,150,20);
		xlabel=new Label("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
		xlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+30,100,20);
		ylabel=new Label("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		ylabel.setBounds(100+256+10+110+10+100,10+imageheight+50+60,100,20);
		xavg=0.0f;
		xavglabel=new Label("x avg = "+xavg);
		xavglabel.setBounds(100+256+10+110+10+100,10+imageheight+50+90,100,20);
		yavg=0.0f;
		yavglabel=new Label("y avg = "+yavg);
		yavglabel.setBounds(100+256+10+110+10+100,10+imageheight+50+120,100,20);
		zavg=0.0f;
		zavglabel=new Label("z avg = "+zavg);
		zavglabel.setBounds(100+256+10+110+10+100,10+imageheight+50+150,100,20);
		xstdev=0.0f;
		xstdevlabel=new Label("x stdev = "+xstdev);
		xstdevlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+180,120,20);
		ystdev=0.0f;
		ystdevlabel=new Label("y stdev = "+ystdev);
		ystdevlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+210,120,20);
		zstdev=0.0f;
		zstdevlabel=new Label("z stdev = "+zstdev);
		zstdevlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+240,120,20);

		roishapegroup=new CheckboxGroup();
		roisquare=new Checkbox("Square roi",roishapegroup,(roishape==0)?true:false);
		roisquare.setBounds(100+256+10,10+imageheight+50+120,100,20);
		roisquare.addItemListener(this);
		roioval=new Checkbox("Oval roi",roishapegroup,(roishape==1)?true:false);
		roioval.setBounds(100+256+10,10+imageheight+50+150,100,20);
		roioval.addItemListener(this);
		plottypegroup=new CheckboxGroup();
		String templabel="var vs. I";
		if(calltype==2){
			templabel="Af vs. Bf";
		}
		if(calltype==3){
			templabel="Acc vs. Don";
		}
		if(calltype==4){
			templabel="S vs. G";
		}
		plotnorm=new Checkbox(templabel,plottypegroup,(plottype==0)?true:false);
		plotnorm.setBounds(10,10+imageheight+180,80,20);
		plotnorm.addItemListener(this);
		templabel="B/S vs. I/S";
		if(calltype==2){
			templabel="Ratio vs. Bf";
		}
		if(calltype==3){
			templabel="Ratio vs. Don";
		}
		if(calltype==4){
			templabel="Ratio vs. G";
		}
		plotdiff=new Checkbox(templabel,plottypegroup,(plottype==1)?true:false);
		plotdiff.setBounds(10,10+imageheight+200,80,20);
		plotdiff.addItemListener(this);
		smoothtypegroup=new CheckboxGroup();
		medsmooth=new Checkbox("Median Sm",smoothtypegroup,(smoothtype==0)?true:false);
		medsmooth.setBounds(10,10+imageheight+240,80,20);
		medsmooth.addItemListener(this);
		gassmooth=new Checkbox("Mean Sm",smoothtypegroup,(smoothtype==1)?true:false);
		gassmooth.setBounds(10,10+imageheight+260,80,20);
		gassmooth.addItemListener(this);
		binsmooth=new Checkbox("Bin Sm",smoothtypegroup,(smoothtype==2)?true:false);
		binsmooth.setBounds(10,10+imageheight+220,80,20);
		binsmooth.addItemListener(this);
		ascalecheck=new Checkbox("Autoscale",(ascale==1)?true:false);
		ascalecheck.setBounds(200,10+imageheight+50+256+30,80,20);
		ascalecheck.addItemListener(this);

		// xmaxval.setBounds(90+245,10+imageheight+50+256+20,80,20);

		multiplier=1.0f;
		scalehistlabel=new Label("Scale Histogram");
		scalehistlabel.setBounds(100+256+10,10+imageheight+50+180,100,20);
		scalehistval=new TextField(""+multiplier);
		scalehistval.setBounds(100+256+10,10+imageheight+50+205,80,20);
		scalehistval.addActionListener(this);
		threshlabel=new Label("Threshhold");
		threshlabel.setBounds(100+newwidth+10,10,150,20);
		threshval=new TextField(""+thresh,10);
		threshval.setBounds(100+newwidth+10,30,80,20);
		threshval.addActionListener(this);
		threshfirstcheck=new Checkbox("Thresh by 1st?",true);
		threshfirstcheck.setBounds(100+newwidth+10+150,60,100,20);
		threshfirstcheck.addItemListener(this);
		dispminlabel=new Label("Minimum Intensity");
		dispminlabel.setBounds(100+newwidth+10,60,150,20);
		dispminval=new TextField(""+dispmin,10);
		dispminval.setBounds(100+newwidth+10,90,80,20);
		dispminval.addActionListener(this);
		dispmaxlabel=new Label("Maximum Intensity");
		dispmaxlabel.setBounds(100+newwidth+10,120,150,20);
		dispmaxval=new TextField(""+dispmax,10);
		dispmaxval.setBounds(100+newwidth+10,150,80,20);
		dispmaxval.addActionListener(this);

		slicelabel=new Label("Slice");
		slicelabel.setBounds(100+newwidth+10+150,90,80,20);
		sliceval=new TextField(""+currslice,10);
		sliceval.setBounds(100+newwidth+10+150,120,80,20);
		sliceval.addActionListener(this);
		slice_slider=new Scrollbar(Scrollbar.HORIZONTAL,currslice,5,0,slices-1+5);
		slice_slider.setBounds(100+newwidth+10+150,150,150,20);
		slice_slider.addAdjustmentListener(this);

		add(xminval);
		add(xmaxval);
		add(yminval);
		add(ymaxval);
		add(scalehistval);
		add(dispminval);
		add(dispminlabel);
		add(threshlabel);
		add(threshval);
		add(threshfirstcheck);
		add(dispmaxval);
		add(dispmaxlabel);
		add(scalehistlabel);
		add(roisquare);
		add(roioval);
		add(roiwidthval);
		add(roiheightval);
		add(roixval);
		add(roiyval);
		add(roiwidthlabel);
		add(roiheightlabel);
		add(smooth_button);
		add(revert_button);

		if(calltype==1){
			add(sval);
			add(s0val);
			add(offval);
			add(backlabel);
			add(backval);
			add(slabel);
			add(s0label);
			add(offlabel);
		}

		add(plotnorm);
		add(plotdiff);
		add(medsmooth);
		add(gassmooth);
		add(binsmooth);
		add(roixlabel);
		add(roiylabel);
		add(saveimg_button);
		add(savehist_button);
		add(scorelabel);
		add(xlabel);
		add(ylabel);
		add(saveyimg_button);
		add(unmix_button);
		add(unmix_button2);
		add(ascalecheck);
		add(slicelabel);
		add(sliceval);
		add(slice_slider);
		add(xavglabel);
		add(yavglabel);
		add(xstdevlabel);
		add(ystdevlabel);
		add(zavglabel);
		add(zstdevlabel);
		histimg=create_2D_histogram();
		update_mask();
		dispimg=create_masked_image();
		repaint();
		set_options();
	}

	public void paint(Graphics g){
		g.drawImage(dispimg,100,10,this);
		g.drawImage(histimg,100,50+imageheight+10,this);
		g.setClip(100,imageheight+50+10,256,256);
		g.setColor(Color.red);
		if(roishape==0){
			g.drawRect(100+roix,(50+imageheight+10+256)-(roiy+roiheight),roiwidth,roiheight);
		}else{
			g.drawOval(100+roix,(50+imageheight+10+256)-(roiy+roiheight),roiwidth,roiheight);
		}
		g.setColor(Color.green);
		float intercept=s0-off*s;
		float y1,y2;
		if(plottype==0){
			y1=s*(xmin+off)+intercept-s0;
			y2=s*(xmax+off)+intercept-s0;
		}else{
			y1=1.0f;
			y2=1.0f;
		}
		g.setColor(Color.green);
		g.drawLine(100,imageheight+50+10+256-(int)(((y1-ymin)/(ymax-ymin))*256.0),100+256,imageheight+50+10+256-(int)(((y2-ymin)/(ymax-ymin))*256.0));
		g.setClip(null);
	}

	public void update(Graphics g){
		paint(g);
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==roisquare||e.getSource()==roioval){
			if(roisquare.getState()==true){
				roishape=0;
			}else{
				roishape=1;
			}
		}
		if(e.getSource()==plotnorm||e.getSource()==plotdiff){
			if(plotnorm.getState()==true){
				plottype=0;
			}else{
				plottype=1;
			}
			histimg=create_2D_histogram();
			update_mask();
			dispimg=create_masked_image();
			xlabel.setText("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		}
		if(e.getSource()==medsmooth||e.getSource()==gassmooth||e.getSource()==binsmooth){
			if(medsmooth.getState()){
				smoothtype=0;
			}else{
				if(gassmooth.getState()){
					smoothtype=1;
				}else{
					smoothtype=2;
				}
			}
		}
		if(e.getSource()==ascalecheck){
			if(ascalecheck.getState()==true){
				ascale=1;
				histimg=create_2D_histogram();
				update_mask();
				dispimg=create_masked_image();
				xlabel.setText("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
				ylabel.setText("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
			}else{
				ascale=0;
			}
		}
		if(e.getSource()==threshfirstcheck){
			if(threshfirstcheck.getState()==true){
				threshfirst=true;
			}else{
				threshfirst=false;
			}
			histimg=create_2D_histogram();
			update_mask();
			dispimg=create_masked_image();
			xlabel.setText("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		}
		repaint();
		set_options();
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==smooth_button){
			smooth();
		}
		if(e.getSource()==revert_button){
			revert();
		}
		if(e.getSource()==savehist_button){
			save_histogram();
		}
		if(e.getSource()==saveimg_button){
			save_masked_image();
		}
		if(e.getSource()==saveyimg_button){
			if(datastack!=null) save_profile();
			else save_histypix();
		}
		if(e.getSource()==unmix_button){
			if(datastack!=null) save_unmixed();
		}
		if(e.getSource()==unmix_button2){
			if(datastack!=null) unmix_geom();
		}
		xmin=Float.parseFloat(xminval.getText());
		xmax=Float.parseFloat(xmaxval.getText());
		ymin=Float.parseFloat(yminval.getText());
		ymax=Float.parseFloat(ymaxval.getText());
		dispmin=Float.parseFloat(dispminval.getText());
		dispmax=Float.parseFloat(dispmaxval.getText());
		multiplier=Float.parseFloat(scalehistval.getText());
		roiwidth=(int)Float.parseFloat(roiwidthval.getText());
		roiheight=(int)Float.parseFloat(roiheightval.getText());
		roix=(int)(Float.parseFloat(roixval.getText())-0.5f*roiwidth);
		roiy=(int)(Float.parseFloat(roiyval.getText())-0.5f*roiheight);
		xlabel.setText("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
		ylabel.setText("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		s=Float.parseFloat(sval.getText());
		s0=Float.parseFloat(s0val.getText());
		off=Float.parseFloat(offval.getText());
		back=Float.parseFloat(backval.getText());
		thresh=Float.parseFloat(threshval.getText());
		currslice=(int)Float.parseFloat(sliceval.getText());
		if(currslice>=slices){
			currslice=slices-1;
		}
		sliceval.setText(""+currslice);
		slice_slider.setValue(currslice);
		histimg=create_2D_histogram();
		update_mask();
		dispimg=create_masked_image();
		repaint();
		set_options();
	}

	public void adjustmentValueChanged(AdjustmentEvent e){
		currslice=slice_slider.getValue();
		sliceval.setText(""+currslice);
		histimg=create_2D_histogram();
		update_mask();
		dispimg=create_masked_image();
		repaint();
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
		int xpos=e.getX();
		int xvalue=xpos-100;
		int ypos=e.getY();
		int yvalue=(256+imageheight+10+50)-ypos;
		if(roishape==0){
			// roi is rectangular
			if(xvalue>roix&&xvalue<(roix+roiwidth)&&yvalue>roiy&&yvalue<(roiy+roiheight)){
				inroi=true;
			}
		}
		if(roishape==1){
			// roi is oval
			// calculate the center of the ellipse and the positions of the foci
			float centerx=roix+0.5f*roiwidth;
			float centery=roiy+0.5f*roiheight;
			float f1x,f1y,f2x,f2y,major;
			if(roiwidth>=roiheight){
				f1y=centery;
				f2y=centery;
				f1x=centerx+(float)Math.sqrt((roiwidth*roiwidth-roiheight*roiheight)/4.0f);
				f2x=centerx-(float)Math.sqrt((roiwidth*roiwidth-roiheight*roiheight)/4.0f);
				major=roiwidth;
			}else{
				f1x=centerx;
				f2x=centerx;
				f1y=centery+(float)Math.sqrt((roiheight*roiheight-roiwidth*roiwidth)/4.0f);
				f2y=centery-(float)Math.sqrt((roiheight*roiheight-roiwidth*roiwidth)/4.0f);
				major=roiheight;
			}
			float distance1,distance2;
			distance1=(float)Math.sqrt((xvalue-f1x)*(xvalue-f1x)+(yvalue-f1y)*(yvalue-f1y));
			distance2=(float)Math.sqrt((xvalue-f2x)*(xvalue-f2x)+(yvalue-f2y)*(yvalue-f2y));
			// IJ.log("f1x "+f1x+"f2x "+f2x);
			if((distance1+distance2)<=major){
				inroi=true;
			}
		}
		oldxpos=xpos;
		oldypos=ypos;
	}

	public void mouseReleased(MouseEvent e){
		inroi=false;
		set_options();
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		if(inroi){
			int xpos=e.getX();
			int ypos=e.getY();
			roix+=xpos-oldxpos;
			roiy-=(ypos-oldypos);
			roixval.setText(""+(roix+roiwidth/2));
			roiyval.setText(""+(roiy+roiheight/2));
			xlabel.setText("x: "+(((roix+roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+(((roiy+roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
			oldxpos=xpos;
			oldypos=ypos;
			update_mask();
			dispimg=create_masked_image();
			repaint();
		}
	}

	float getminlesszero(float[] farr){
		float min;
		min=0.0f;
		for(int i=0;i<farr.length;i++){
			if(farr[i]<min&&!Float.isInfinite(farr[i])){
				min=farr[i];
			}
		}
		return min;
	}

	float getmaxgreatzero(float[] farr){
		float max;
		max=0.0f;
		for(int i=0;i<farr.length;i++){
			if(farr[i]>max&&!Float.isInfinite(farr[i])){
				max=farr[i];
			}
		}
		return max;
	}

	Image create_masked_image(){
		int[] color_pixels=new int[width*height];
		int dumint;
		for(int i=0;i<width*height;i++){
			dumint=(int)(((disppix[i+currslice*width*height]-dispmin)/(dispmax-dispmin))*256.0);
			if(threshfirst){
				if(disppix[i]<thresh){
					dumint=0;
				}
			}else{
				if(disppix[i+currslice*width*height]<thresh){
					dumint=0;
				}
			}
			if(dumint<0){
				dumint=0;
			}
			if(dumint>255){
				dumint=255;
			}
			if(disppix[i+currslice*width*height]>=dispmin){
				if(mask[i]==0){
					color_pixels[i]=0xff000000+(dumint<<16)+(dumint<<8)+dumint;
				}else{
					color_pixels[i]=0xffff0000;
				}
			}else{
				color_pixels[i]=0xff000000;
			}
		}
		// if the image is larger than 256, scale it down
		int scalefactor=1+(int)((width-1)/256.0f);
		dumint=1+(int)((height-1)/256.0f);
		if(dumint>scalefactor){
			scalefactor=dumint;
		}
		int newwidth=width/scalefactor;
		int newheight=height/scalefactor;
		if(scalefactor>1){
			int[] newcolorpixels=new int[newwidth*newheight];
			for(int i=0;i<newheight;i++){
				for(int j=0;j<newwidth;j++){
					boolean masked=false;
					float avg=0.0f;
					local_search: for(int k=0;k<scalefactor;k++){
						for(int l=0;l<scalefactor;l++){
							dumint=j*scalefactor+l+(i*scalefactor+k)*width;
							if(mask[dumint]==1){
								masked=true;
								break local_search;
							}
							avg+=(float)(color_pixels[dumint]&0x000000ff)/(float)(scalefactor*scalefactor);
						}
					}
					dumint=(int)avg;
					if(masked){
						newcolorpixels[i*newwidth+j]=0xffff0000;
					}else{
						newcolorpixels[i*newwidth+j]=0xff000000+(dumint<<16)+(dumint<<8)+dumint;
					}
				}
			}
			color_pixels=newcolorpixels;
		}
		ColorProcessor cp=new ColorProcessor(newwidth,newheight,color_pixels);
		return cp.createImage();
	}

	void update_mask(){
		int dumint1,dumint2;
		// clear the mask image
		mask=new int[width*height];
		for(int i=0;i<width*height;i++){
			mask[i]=0;
		}
		// calculate the mask based on histogram hit testing
		score=0.0f;
		xavg=0.0f;
		yavg=0.0f;
		zavg=0.0f;
		xstdev=0.0f;
		ystdev=0.0f;
		zstdev=0.0f;
		float totpixels=0.0f;
		for(int i=0;i<width*height;i++){
			dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
			dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
			boolean abovethresh=true;
			if(threshfirst){
				if(disppix[i]<thresh){
					abovethresh=false;
				}
			}else{
				if(disppix[i+currslice*width*height]<thresh){
					abovethresh=false;
				}
			}
			if(((dumint1<256&&dumint2<256)&&(dumint1>0&&dumint2>0))&&abovethresh){
				if(roishape==0){
					// roi is rectangular
					if((dumint1>roix&&dumint1<(roix+roiwidth))&&(dumint2>roiy&&dumint2<(roiy+roiheight))){
						mask[i]=1;
						score+=1.0;
						xavg+=histxpix[i];
						yavg+=histypix[i];
						xstdev+=histxpix[i]*histxpix[i];
						ystdev+=histypix[i]*histypix[i];
						if(histzpix!=null){
							zavg+=histzpix[i];
							zstdev+=histzpix[i]*histzpix[i];
						}
					}
				}
				if(roishape==1){
					// roi is oval
					// calculate the center of the ellipse and the positions of
					// the foci
					float centerx=roix+0.5f*roiwidth;
					float centery=roiy+0.5f*roiheight;
					float f1x,f1y,f2x,f2y,major;
					if(roiwidth>=roiheight){
						f1y=centery;
						f2y=centery;
						f1x=centerx+(float)Math.sqrt((roiwidth*roiwidth-roiheight*roiheight)/4.0f);
						f2x=centerx-(float)Math.sqrt((roiwidth*roiwidth-roiheight*roiheight)/4.0f);
						major=roiwidth;
					}else{
						f1x=centerx;
						f2x=centerx;
						f1y=centery+(float)Math.sqrt((roiheight*roiheight-roiwidth*roiwidth)/4.0f);
						f2y=centery-(float)Math.sqrt((roiheight*roiheight-roiwidth*roiwidth)/4.0f);
						major=roiwidth;
					}
					float distance1,distance2;
					distance1=(float)Math.sqrt((dumint1-f1x)*(dumint1-f1x)+(dumint2-f1y)*(dumint2-f1y));
					distance2=(float)Math.sqrt((dumint1-f2x)*(dumint1-f2x)+(dumint2-f2y)*(dumint2-f2y));
					if((distance1+distance2)<=major){
						mask[i]=1;
						score+=1.0;
						xavg+=histxpix[i];
						yavg+=histypix[i];
						xstdev+=histxpix[i]*histxpix[i];
						ystdev+=histypix[i]*histypix[i];
						if(histzpix!=null){
							zavg+=histzpix[i];
							zstdev+=histzpix[i]*histzpix[i];
						}
					}
				}
			}
			if(abovethresh){
				totpixels+=1.0f;
			}
		}
		xavg/=score;
		yavg/=score;
		zavg/=score;
		xstdev/=score;
		ystdev/=score;
		zstdev/=score;
		xstdev-=xavg*xavg;
		ystdev-=yavg*yavg;
		zstdev-=zavg*zavg;
		xstdev*=(1.0f-1.0f/score);
		ystdev*=(1.0f-1.0f/score);
		zstdev*=(1.0f-1.0f/score);
		xstdev=(float)Math.sqrt(xstdev);
		ystdev=(float)Math.sqrt(ystdev);
		zstdev=(float)Math.sqrt(zstdev);
		score/=totpixels;
		scorelabel.setText("Score = "+score);
		xavglabel.setText("x avg = "+xavg);
		yavglabel.setText("y avg = "+yavg);
		xstdevlabel.setText("x stdev = "+xstdev);
		ystdevlabel.setText("y stdev = "+ystdev);
		zavglabel.setText("z avg = "+zavg);
		zstdevlabel.setText("z stdev = "+zstdev);
	}

	Image create_2D_histogram(){
		int dumint1,dumint2;
		update_histpix();
		int[][] histogram=new int[256][256];
		// calculate the histogram
		for(int i=0;i<width*height;i++){
			dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
			dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
			boolean abovethresh=true;
			if(threshfirst){
				if(disppix[i]<thresh){
					abovethresh=false;
				}
			}else{
				if(disppix[i+currslice*width*height]<thresh){
					abovethresh=false;
				}
			}
			if(abovethresh){
				if((dumint1<256&&dumint2<256)&&(dumint1>0&&dumint2>0)){
					histogram[dumint2][dumint1]++;
				}
			}
		}
		// transfer it to a single array for output
		float[] tempfloat=new float[256*256];
		histmax=0;
		for(int i=0;i<256;i++){
			for(int j=0;j<256;j++){
				tempfloat[(255-i)*256+j]=histogram[i][j];
				if(histogram[i][j]>histmax){
					histmax=histogram[i][j];
				}
			}
		}
		FloatProcessor fp=new FloatProcessor(256,256,tempfloat,null);
		fp.setMinAndMax(0.0,(histmax)/multiplier);
		update_mask();
		return fp.createImage();
	}

	void smooth(){
		if(smoothtype!=2){
			// here we use a median or gaussian smoothing method
			// the weighting is 2 for the pixel and 1 for its neighbors
			// the pixels around the edge are not smoothed
			int i,j,k,l;
			float dumflt;
			float[] pixelvals=new float[10];
			// first smooth the x image
			for(int m=0;m<slices;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothxpix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothxpix[i*width+j+m*width*height];
						}
						if(smoothtype==0){
							smoothxpix[i*width+j+m*width*height]=median(pixelvals);
						}else{
							smoothxpix[i*width+j+m*width*height]=avg(pixelvals);
						}
					}
				}
			}
			// now smooth the y image
			for(int m=0;m<slices;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothypix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothypix[i*width+j+m*width*height];
						}
						if(smoothtype==0){
							smoothypix[i*width+j+m*width*height]=median(pixelvals);
						}else{
							smoothypix[i*width+j+m*width*height]=avg(pixelvals);
						}
					}
				}
			}
			if(smoothzpix!=null){
				for(int m=0;m<slices;m++){
					for(i=1;i<(height-1);i++){
						for(j=1;j<(width-1);j++){
							for(k=0;k<3;k++){
								for(l=0;l<3;l++){
									pixelvals[3*k+l]=smoothzpix[(i-(k-1))*width+j-(l-1)+m*width*height];
								}
								pixelvals[9]=smoothzpix[i*width+j+m*width*height];
							}
							if(smoothtype==0){
								smoothzpix[i*width+j+m*width*height]=median(pixelvals);
							}else{
								smoothzpix[i*width+j+m*width*height]=avg(pixelvals);
							}
						}
					}
				}
			}
		}else{
			GenericDialog gd=new GenericDialog("Bin Size");
			int binsize=2;
			gd.addNumericField("Bin Size?",binsize,0);
			gd.showDialog();
			if(!gd.wasCanceled()){
				binsize=(int)gd.getNextNumber();
				for(int m=0;m<slices;m++){
					for(int i=0;i<height/binsize;i++){
						for(int j=0;j<width/binsize;j++){
							float sumx=0.0f;
							float sumy=0.0f;
							float sumz=0.0f;
							int npix=0;
							for(int k=0;k<binsize;k++){
								for(int l=0;l<binsize;l++){
									boolean abovethresh=true;
									if(threshfirst){
										if(disppix[j*binsize+l+(i*binsize+k)*width]<thresh){
											abovethresh=false;
										}
									}else{
										if(disppix[j*binsize+l+(i*binsize+k)*width+m*width*height]<thresh){
											abovethresh=false;
										}
									}
									if(abovethresh){
										npix++;
										sumx+=smoothxpix[j*binsize+l+(i*binsize+k)*width+m*width*height];
										sumy+=smoothypix[j*binsize+l+(i*binsize+k)*width+m*width*height];
										if(smoothzpix!=null)
											sumz+=smoothzpix[j*binsize+l+(i*binsize+k)*width+m*width*height];
									}
								}
							}
							if(npix>0){
								sumx/=npix;
								sumy/=npix;
								sumz/=npix;
							}
							for(int k=0;k<binsize;k++){
								for(int l=0;l<binsize;l++){
									smoothxpix[j*binsize+l+(i*binsize+k)*width+m*width*height]=sumx;
									smoothypix[j*binsize+l+(i*binsize+k)*width+m*width*height]=sumy;
									if(smoothzpix!=null)
										smoothzpix[j*binsize+l+(i*binsize+k)*width+m*width*height]=sumz;
								}
							}
						}
					}
				}
			}
		}
		histimg=create_2D_histogram();
		update_mask();
		dispimg=create_masked_image();
		repaint();
	}

	void revert(){
		int i;
		smoothxpix=xpix.clone();
		smoothypix=ypix.clone();
		if(zpix!=null)
			smoothzpix=zpix.clone();
	}

	float median(float[] values){
		int i,length,dumint,imin,j;
		float min,dumfloat;
		length=values.length;
		dumint=1+(int)(length/2.0f);
		// sort the vector one at a time
		for(j=0;j<dumint;j++){
			imin=j;
			min=values[imin];
			for(i=imin+1;i<length;i++){
				if(values[i]<min){
					imin=i;
					min=values[imin];
				}
			}
			dumfloat=values[j];
			values[j]=values[imin];
			values[imin]=dumfloat;
		}
		if((dumint*2)>length){
			return values[dumint-1];
		}else{
			min=values[dumint];
			for(i=dumint+1;i<length;i++){
				if(values[i]<min){
					min=values[i];
				}
			}
			return 0.5f*values[dumint-1]*min;
		}
	}

	float avg(float[] values){
		int length;
		float dumfloat=0.0f;
		length=values.length;
		for(int i=0;i<length;i++){
			dumfloat+=values[i]/length;
		}
		return dumfloat;
	}

	void save_histogram(){
		GenericDialog gdsave=new GenericDialog("Options");
		boolean saveashist=true;
		gdsave.addCheckbox("Save as histogram (will sum slices)?",saveashist);
		gdsave.showDialog();
		if(gdsave.wasCanceled()){
			return;
		}
		saveashist=gdsave.getNextBoolean();
		if(!saveashist){
			int ascaled,oldslice;
			ImageStack out_stack=new ImageStack(256,256);
			if(ascale==1){
				ascaled=1;
			}else{
				ascaled=0;
			}
			ascale=0;
			oldslice=currslice;
			for(int k=0;k<slices;k++){
				currslice=k;
				update_histpix();
				int[][] histogram=new int[256][256];
				// calculate the histogram
				for(int i=0;i<width*height;i++){
					int dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
					int dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
					boolean abovethresh=true;
					if(threshfirst){
						if(disppix[i]<thresh){
							abovethresh=false;
						}
					}else{
						if(disppix[i+currslice*width*height]<thresh){
							abovethresh=false;
						}
					}
					if(abovethresh){
						if((dumint1<256&&dumint2<256)&&(dumint1>0&&dumint2>0)){
							histogram[dumint2][dumint1]++;
						}
					}
				}
				// transfer it to a single array for output
				float[] tempfloat=new float[256*256];
				for(int i=0;i<256;i++){
					for(int j=0;j<256;j++){
						tempfloat[(255-i)*256+j]=histogram[i][j];
					}
				}
				out_stack.addSlice("",tempfloat);
			}
			if(ascaled==1){
				ascale=1;
			}
			currslice=oldslice;
			ImagePlus imp=new ImagePlus("2D histogram",out_stack);
			Calibration cal=imp.getCalibration();
			cal.pixelWidth=(double)(xmax-xmin)/256.0f;
			cal.pixelHeight=(double)(ymax-ymin)/256.0f;
			cal.xOrigin=-(double)xmin/cal.pixelWidth;
			cal.yOrigin=ymin/cal.pixelHeight+255.0;
			cal.setInvertY(true);
			// IJ.log("ymax = "+ymax);
			cal.setUnit("nb");
			imp.setRoi(new Roi(roix,256-roiy-roiheight,roiwidth,roiheight));
			imp.show();
		}else{
			int ascaled,oldslice;
			if(ascale==1){
				ascaled=1;
			}else{
				ascaled=0;
			}
			ascale=0;
			oldslice=currslice;
			float[] outhistxpix=new float[slices*width*height];
			float[] outhistypix=new float[slices*width*height];
			int counter=0;
			for(int k=0;k<slices;k++){
				currslice=k;
				update_histpix();
				// calculate the histogram
				for(int i=0;i<width*height;i++){
					boolean abovethresh=true;
					if(threshfirst){
						if(disppix[i]<thresh){
							abovethresh=false;
						}
					}else{
						if(disppix[i+currslice*width*height]<thresh){
							abovethresh=false;
						}
					}
					if(abovethresh){
						outhistxpix[counter]=histxpix[i];
						outhistypix[counter]=histypix[i];
						counter++;
					}
				}
			}
			if(ascaled==1){
				ascale=1;
			}
			currslice=oldslice;
			float[] outhistxpix2=new float[counter];
			float[] outhistypix2=new float[counter];
			System.arraycopy(outhistxpix,0,outhistxpix2,0,counter);
			System.arraycopy(outhistypix,0,outhistypix2,0,counter);
			PlotWindow2DHist pw=new PlotWindow2DHist("2D Histogram","x","y",outhistxpix2,outhistypix2,null);
			pw.draw();
			float[] limits=pw.getLimits();
			limits[0]=xmin;
			limits[1]=xmax;
			limits[2]=ymin;
			limits[3]=ymax;
			pw.setLimits(limits);
			pw.intautoscale();
			//pw.draw();
		}
	}

	void update_histpix(){
		histxpix=new float[width*height];
		histypix=new float[width*height];
		histzpix=null;
		if(zpix!=null)
			histzpix=new float[width*height];
		if(plottype==0){
			for(int i=0;i<width*height;i++){
				histxpix[i]=smoothxpix[i+currslice*width*height]-off;
				histypix[i]=smoothypix[i+currslice*width*height]-s0;
				if(zpix!=null)
					histzpix[i]=smoothzpix[i+currslice*width*height];
			}
		}else{
			for(int i=0;i<width*height;i++){
				histxpix[i]=(smoothxpix[i+currslice*width*height]-off-back)/s;
				histypix[i]=(smoothypix[i+currslice*width*height]-s0-back*s)/(histxpix[i]*s*s);
				if(zpix!=null)
					histzpix[i]=smoothzpix[i+currslice*width*height];
			}
		}
		if(ascale==1){
			xmax=getmaxgreatzero(histxpix)+0.1f;
			xmaxval.setText(""+xmax);
			ymax=getmaxgreatzero(histypix)+0.1f;
			ymaxval.setText(""+ymax);
			xmin=getminlesszero(histxpix);
			xminval.setText(""+xmin);
			ymin=getminlesszero(histypix);
			yminval.setText(""+ymin);
		}
	}

	void save_masked_image(){
		int[] color_pixels=new int[width*height];
		float[] float_pix=new float[width*height];
		System.arraycopy(disppix,currslice*width*height,float_pix,0,width*height);
		int dumint;
		for(int i=0;i<width*height;i++){
			dumint=(int)(((disppix[i+currslice*width*height]-dispmin)/(dispmax-dispmin))*256.0);
			if(threshfirst){
				if(disppix[i]<thresh){
					dumint=0;
				}
			}else{
				if(disppix[i+currslice*width*height]<thresh){
					dumint=0;
				}
			}
			if(dumint<0){
				dumint=0;
			}
			if(dumint>255){
				dumint=255;
			}
			if(disppix[i+currslice*width*height]>=dispmin){
				if(mask[i]==0){
					//color_pixels[i]=0xff000000+(dumint<<16)+(dumint<<8)+dumint; //grayscale image
					color_pixels[i]=0xff000000; //black
				}else{
					color_pixels[i]=0xffff0000; //red
				}
			}else{
				color_pixels[i]=0xff000000; //black
			}
		}
		ImagePlus imp2=new ImagePlus("Masked Image",new FloatProcessor(width,height,float_pix,null));
		ImageRoi maskroi=new ImageRoi(0,0,new ColorProcessor(width,height,color_pixels));
		maskroi.setZeroTransparent(true);
		imp2.setCalibration(cal);
		imp2.show();
		imp2.setOverlay(new Overlay(maskroi));
	}

	void save_histypix(){
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Create Heate Map?",true);
		gd.addNumericField("Minimum Value",0.9,5,15,null);
		gd.addNumericField("Maximum Value",1.3,5,15,null);
		gd.addCheckbox("Output Calibration Bar",true);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		boolean heatmap=gd.getNextBoolean();
		float heatmin=(float)gd.getNextNumber();
		float heatmax=(float)gd.getNextNumber();
		boolean outcal=gd.getNextBoolean();
		int tempcurrslice=currslice;
		ImageStack histystack=new ImageStack(width,height);
		for(int j=0;j<slices;j++){
			float[] filteredypix=new float[width*height];
			currslice=j;
			update_histpix();
			for(int i=0;i<width*height;i++){
				boolean abovethresh=true;
				if(threshfirst){
					if(disppix[i]<thresh){
						abovethresh=false;
					}
				}else{
					if(disppix[i+currslice*width*height]<thresh){
						abovethresh=false;
					}
				}
				if(abovethresh){
					if((histxpix[i]<xmax&&histxpix[i]>xmin)&&(histypix[i]<ymax&&histypix[i]>ymin)){
						filteredypix[i]=histypix[i];
					}else{
						filteredypix[i]=0.0f;
					}
				}
			}
			histystack.addSlice("",filteredypix);
		}
		if(!heatmap){
			ImagePlus imphisty=new ImagePlus("Y Image",histystack);
			imphisty.setCalibration(cal);
			imphisty.show();
		}else{
			int[] heatimage=new int[width*height];
			float[] temp=(float[])histystack.getPixels(currslice+1);
			for(int i=0;i<width*height;i++){
				float temp2=(temp[i]-heatmin)/(heatmax-heatmin);
				if(temp[i]<=0){
					heatimage[i]=0xff000000;
				}else{
					heatimage[i]=getgreenyellowred(temp2);
				}
			}
			ImagePlus imphisty=new ImagePlus("Heat Map",new ColorProcessor(width,height,heatimage));
			imphisty.setCalibration(cal);
			imphisty.show();
			if(outcal)
				show_cal_bar(heatmin,heatmax);
		}
		currslice=tempcurrslice;
	}

	int getgreenyellowred(float val){
		// here value is between zero and 1
		if(val<=0f){
			return 0xff00ff00;
		}else{
			if(val>1f){
				return 0xffff0000;
			}else{
				if(val<0.5f){
					int red=(int)(255.0f*val/0.5f);
					return 0xff000000+(red<<16)+0xff00;
				}else{
					int green=(int)(255.0f*(0.5f-val)/0.5f);
					return 0xffff0000+(green<<8);
				}
			}
		}
	}

	void show_cal_bar(float heatmin,float heatmax){
		ColorProcessor cp=new ColorProcessor(60,160);
		cp.setColor(Color.WHITE);
		cp.fill();
		cp.setColor(Color.BLACK);
		cp.drawRect(4,10,23,142);
		for(int i=0;i<140;i++){
			float val=i/139f;
			Color temp=new Color(getgreenyellowred(val));
			cp.setColor(temp);
			cp.drawLine(5,150-i,25,150-i);
		}
		float heatmiddle=0.5f*(heatmax+heatmin);
		cp.setFont(new Font("SansSerif",Font.BOLD,12));
		cp.setColor(Color.BLACK);
		cp.drawString(IJ.d2s(heatmin,1),30,157);
		cp.drawString(IJ.d2s(heatmiddle,1),30,87);
		cp.drawString(IJ.d2s(heatmax,1),30,17);
		new ImagePlus("Heat Map Calibration",cp).show();
	}
	
	void save_profile(){
		if(datastack!=null){
			int dtype=algutils.get_array_type(datastack.getPixels(1));
			int stacklength=datastack.getSize()/slices;
			float[] decay=new float[stacklength];
			if(dtype==2){
				for(int i=0;i<width*height;i++){
					if(mask[i]==1){
						for(int j=0;j<stacklength;j++){
							decay[j]+=((float[])datastack.getPixels(currslice*stacklength+j+1))[i];
						}
					}
				}
			} else if(dtype==1) {
				for(int i=0;i<width*height;i++){
					if(mask[i]==1){
						for(int j=0;j<stacklength;j++){
							float temp=(float)(((short[])datastack.getPixels(currslice*stacklength+j+1))[i]&0xffff);
							decay[j]+=temp;
						}
					}
				}
			} else {
				for(int i=0;i<width*height;i++){
					if(mask[i]==1){
						for(int j=0;j<stacklength;j++){
							float temp=(float)(((byte[])datastack.getPixels(currslice*stacklength+j+1))[i]&0xff);
							decay[j]+=temp;
						}
					}
				}
			}
			(new PlotWindow4("Masked Profile","channel","Intensity",decay)).draw();
		}
	}
	
	void save_unmixed(){
		//start by getting the spectra
		ImageWindow[] iw=jutils.selectPlots(false,1,new String[]{"Reference Spectra"});
		if(iw==null) return;
		float[][] spectra=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
		float[][] unmixed=PU.unmix_phasor(spectra,this);
		ImagePlus unimp=new ImagePlus("Unmixed",jutils.array2stack(unmixed,width,height));
		unimp.setOpenAsHyperStack(true);
		unimp.setDimensions(unmixed.length/slices,slices,1);
		new CompositeImage(unimp,CompositeImage.COLOR).show();
	}
	
	void unmix_geom(){
		GenericDialog gd2=new GenericDialog("Number_of_components");
		gd2.addNumericField("How_Many_Components?",3,0);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		int ncomp=(int)gd2.getNextNumber();
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<ncomp;i++){
			gd.addNumericField("G"+(i+1),0.0,5,15,null);
			gd.addNumericField("S"+(i+1),0.0,5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()) return;
		float[][] positions=new float[ncomp][2];
		for(int i=0;i<ncomp;i++){
			positions[i][0]=(float)gd.getNextNumber();
			positions[i][1]=(float)gd.getNextNumber();
		}
		Object[] temp=PU.unmix_phasor_geom3(positions,this);
		float[][] unmixed=(float[][])temp[0];
		float[][] fractions=(float[][])temp[1];
		ImagePlus unimp=new ImagePlus("Unmixed",jutils.array2stack(unmixed,width,height));
		unimp.setOpenAsHyperStack(true);
		unimp.setDimensions(unmixed.length/slices,slices,1);
		new CompositeImage(unimp,CompositeImage.COLOR).show();
		new ImagePlus("Fractions",jutils.array2stack(fractions,256,256)).show();
	}
	
	void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"histogram_2D_N_B_jru_v1.jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			roix=Integer.parseInt(d.readLine());
			roiy=Integer.parseInt(d.readLine());
			roiwidth=Integer.parseInt(d.readLine());
			roiheight=Integer.parseInt(d.readLine());
			plottype=Integer.parseInt(d.readLine());
			smoothtype=Integer.parseInt(d.readLine());
			roishape=Integer.parseInt(d.readLine());
			ascale=Integer.parseInt(d.readLine());
			threshfirst=(Integer.parseInt(d.readLine())==0)?false:true;
			xmin=Float.parseFloat(d.readLine());
			xmax=Float.parseFloat(d.readLine());
			ymin=Float.parseFloat(d.readLine());
			ymax=Float.parseFloat(d.readLine());
			s=Float.parseFloat(d.readLine());
			s0=Float.parseFloat(d.readLine());
			off=Float.parseFloat(d.readLine());
			back=Float.parseFloat(d.readLine());
			d.close();
		}catch(IOException e){
			roix=0;
			roiy=0;
			roiwidth=30;
			roiheight=30;
			plottype=0;
			smoothtype=1;
			roishape=0;
			ascale=1;
			threshfirst=true;
			xmin=0.0f;
			xmax=100.0f;
			ymin=0.0f;
			ymax=1000000.0f;
			s=200.0f;
			s0=1000.0f;
			off=2500.0f;
			back=0.0f;
		}
		return;
	}

	void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){
				a.mkdir();
			}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"histogram_2D_N_B_jru_v1.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+roix+"\n");
			d.write(""+roiy+"\n");
			d.write(""+roiwidth+"\n");
			d.write(""+roiheight+"\n");
			d.write(""+plottype+"\n");
			d.write(""+smoothtype+"\n");
			d.write(""+roishape+"\n");
			d.write(""+ascale+"\n");
			d.write(""+(threshfirst?1:0)+"\n");
			d.write(""+xmin+"\n");
			d.write(""+xmax+"\n");
			d.write(""+ymin+"\n");
			d.write(""+ymax+"\n");
			d.write(""+s+"\n");
			d.write(""+s0+"\n");
			d.write(""+off+"\n");
			d.write(""+back+"\n");
			d.close();
		}catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

}
