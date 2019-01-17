/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import ij.measure.*;
import java.io.*;
import jguis.*;

public class histogram_2D_N_B_jru_v1 implements PlugIn {

	public void run(String arg) {

		int i;
		final int numimages=3;
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}

		GenericDialog gd = new GenericDialog("Choose Images");
		boolean calcavgvar=true;
		gd.addCheckbox("Calculate Avg and Var?",calcavgvar);
		boolean avgall=true;
		gd.addCheckbox("Avg all frames?",avgall);
		int nframes=100;
		gd.addNumericField("# Frames (for Avg Var Calc)?",nframes,0);
		gd.addChoice("Image Stack (for Avg Var Calc)",titles,titles[0]);
		gd.addChoice("Avg image",titles,titles[0]);
		gd.addChoice("Var image",titles,titles[0]);
		gd.addChoice("display image",titles,titles[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int[] index=new int[4];
		calcavgvar=gd.getNextBoolean();
		avgall=gd.getNextBoolean();
		nframes=(int)gd.getNextNumber();
		for(i=0;i<4;i++){
			index[i]=gd.getNextChoiceIndex();
		}
		ImagePlus impstack=WindowManager.getImage(wList[index[0]]);
		ImagePlus imp1 = WindowManager.getImage(wList[index[1]]);
		ImagePlus imp2 = WindowManager.getImage(wList[index[2]]);
		ImagePlus imp3 = WindowManager.getImage(wList[index[3]]);
		if(calcavgvar){
			if(avgall){nframes=impstack.getStack().getSize();}
			ImagePlus[] imps=stack2avgvar(impstack.getStack(),nframes);
			imp1=imps[0];
			imp2=imps[1];
			imp3=imps[2];
		}
		
		/*final CustomWindow cw = new CustomWindow();
		cw.init(imp1,imp2,imp3);

		final  Frame f = new Frame("Interactive 2D histogram");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				f.dispose();
			}
		});

		f.add(cw);
		f.pack(); 
		f.setResizable(false);
		Insets ins = f.getInsets();
		cw.totalSize.height = CustomWindow.H + ins.bottom + ins.top + 65;
		cw.totalSize.width  = CustomWindow.WR + ins.left + ins.right;
		f.setSize(cw.totalSize);
		f.setVisible(true);
		cw.requestFocus();*/
		final Hist2DWindow cw=new Hist2DWindow();
		cw.init(imp1,imp2,imp3,1);
		Hist2DWindow.launch_frame("Interactive 2D Histogram",cw);
	}

	ImagePlus[] stack2avgvar(ImageStack stack,int nframes){
		int width=stack.getWidth();
		int height=stack.getHeight();
		int slices=stack.getSize();
		int result_slices=1;
		ImageProcessor ip=stack.getProcessor(1);
		if(slices!=nframes){
			result_slices=(int)((float)slices/(float)nframes);
		}
		ImageStack avg_stack=new ImageStack(width,height); ImageStack var_stack=new ImageStack(width,height);
		for(int i=0;i<result_slices;i++){
			float[] avg=new float[width*height];
			float[] var=new float[width*height];
			for(int j=0;j<nframes;j++){
				if(ip instanceof FloatProcessor){
					float[] temp=(float[])stack.getPixels(j+i*nframes+1);
					for(int k=0;k<(width*height);k++){
						avg[k]+=temp[k]/(float)nframes;
						var[k]+=(temp[k]*temp[k])/(float)nframes;
					}
				}
				if(ip instanceof ShortProcessor){
					short[] temp=(short[])stack.getPixels(j+i*nframes+1);
					for(int k=0;k<(width*height);k++){
						float temp2=temp[k]&0xffff;
						avg[k]+=temp2/(float)nframes;
						var[k]+=(temp2*temp2)/(float)nframes;
					}
				}
			}
			for(int j=0;j<(width*height);j++){
				var[j]-=avg[j]*avg[j];
			}
			avg_stack.addSlice("",(Object)avg);
			var_stack.addSlice("",(Object)var);
		}
		ImagePlus[] imps=new ImagePlus[3];
		imps[0]=new ImagePlus("avg",avg_stack);
		imps[2]=imps[0];
		imps[1]=new ImagePlus("var",var_stack);
		return imps;
	}

}

class CustomWindow extends Panel implements ActionListener, AdjustmentListener, MouseListener, MouseMotionListener, ItemListener {

	public final static int H = 612;
	public final static int WR = 800;
	public Dimension totalSize = new Dimension();
	private TextField scalehistval,dispminval,dispmaxval,threshval;
	private TextField xminval,xmaxval,yminval,ymaxval;
	private TextField roiwidthval, roiheightval,roixval,roiyval;
	private TextField sval,s0val,offval,backval,sliceval;
	private Label dispminlabel,dispmaxlabel,scalehistlabel,scorelabel,xlabel,ylabel,slicelabel,xavglabel,yavglabel,threshlabel;
	private Label xstdevlabel,ystdevlabel;
	private Label roiwidthlabel,roiheightlabel,roixlabel,roiylabel,slabel,s0label,offlabel,backlabel;
	private Button smooth_button,revert_button,savehist_button,saveimg_button,saveyimg_button;
	CheckboxGroup roishapegroup;
	Checkbox roisquare, roioval;
	CheckboxGroup plottypegroup;
	Checkbox plotnorm,plotdiff;
	CheckboxGroup smoothtypegroup;
	Checkbox medsmooth,gassmooth,binsmooth;
	Checkbox ascalecheck;
	Checkbox threshfirstcheck;
	Calibration cal;
	Scrollbar slice_slider;
	Image dispimg, histimg,tempimg;
	Graphics tempgraphics;
	float xmin,xmax,ymin,ymax,dispmin,dispmax,multiplier,s,off,s0,score,back;
	float xavg,yavg,thresh,xstdev,ystdev;
	int width,height,roishape,roix,roiy,roiwidth,roiheight,slices,currslice;
	int oldxpos,oldypos,histmax,plottype,pct,smoothtype,ascale;
	float[] xpix,ypix,histxpix,histypix,disppix;
	float[] smoothxpix,smoothypix;
	int[] mask;
	boolean inroi,threshfirst;
	int imageheight=256;

	void init(ImagePlus ximp, ImagePlus yimp, ImagePlus dispimp){
		init_options();
		setLayout(null);
		//initialize all of the variables
		width=ximp.getWidth();
		int newwidth=width;
		if(newwidth>256){newwidth=256;}
		height=ximp.getHeight();
		ImageStack xstack=ximp.getStack();
		ImageStack ystack=yimp.getStack();
		ImageStack dispstack=dispimp.getStack();
		cal=dispimp.getCalibration();
		slices=xstack.getSize();
		xpix=new float[width*height*slices];
		ypix=new float[width*height*slices];
		disppix=new float[width*height*slices];
		for(int i=0;i<slices;i++){
			float[] temp=(float[])xstack.getPixels(i+1);
			for(int j=0;j<width*height;j++){xpix[j+i*width*height]=temp[j];}
			temp=(float[])ystack.getPixels(i+1);
			for(int j=0;j<width*height;j++){ypix[j+i*width*height]=temp[j];}
			temp=(float[])dispstack.getPixels(i+1);
			for(int j=0;j<width*height;j++){disppix[j+i*width*height]=temp[j];}
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
		saveyimg_button=new Button("Save Y Image");
		saveyimg_button.setBounds(100+256+10+110+10+100,10+imageheight+50+200,80,20);
		saveyimg_button.addActionListener(this);

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
		xlabel=new Label("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
		xlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+30,100,20);
		ylabel=new Label("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		ylabel.setBounds(100+256+10+110+10+100,10+imageheight+50+60,100,20);
		xavg=0.0f;
		xavglabel=new Label("x avg = "+xavg);
		xavglabel.setBounds(100+256+10+110+10+100,10+imageheight+50+90,100,20);
		yavg=0.0f;
		yavglabel=new Label("y avg = "+yavg);
		yavglabel.setBounds(100+256+10+110+10+100,10+imageheight+50+120,100,20);
		xstdev=0.0f;
		xstdevlabel=new Label("x stdev = "+xstdev);
		xstdevlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+150,120,20);
		ystdev=0.0f;
		ystdevlabel=new Label("y stdev = "+ystdev);
		ystdevlabel.setBounds(100+256+10+110+10+100,10+imageheight+50+180,120,20);

		roishapegroup=new CheckboxGroup();
		roisquare=new Checkbox("Square roi",roishapegroup,(roishape==0) ? true : false);
		roisquare.setBounds(100+256+10,10+imageheight+50+120,100,20);
		roisquare.addItemListener(this);
		roioval=new Checkbox("Oval roi",roishapegroup,(roishape==1) ? true : false);
		roioval.setBounds(100+256+10,10+imageheight+50+150,100,20);
		roioval.addItemListener(this);
		plottypegroup=new CheckboxGroup();
		plotnorm=new Checkbox("var vs. I",plottypegroup,(plottype==0) ? true : false);
		plotnorm.setBounds(10,10+imageheight+180,80,20);
		plotnorm.addItemListener(this);
		plotdiff=new Checkbox("B/S vs. I/S",plottypegroup,(plottype==1) ? true : false);
		plotdiff.setBounds(10,10+imageheight+200,80,20);
		plotdiff.addItemListener(this);
		smoothtypegroup=new CheckboxGroup();
		medsmooth=new Checkbox("Median Sm",smoothtypegroup,(smoothtype==0) ? true : false);
		medsmooth.setBounds(10,10+imageheight+240,80,20);
		medsmooth.addItemListener(this);
		gassmooth=new Checkbox("Mean Sm",smoothtypegroup,(smoothtype==1) ? true : false);
		gassmooth.setBounds(10,10+imageheight+260,80,20);
		gassmooth.addItemListener(this);
		binsmooth=new Checkbox("Bin Sm",smoothtypegroup,(smoothtype==2) ? true : false);
		binsmooth.setBounds(10,10+imageheight+220,80,20);
		binsmooth.addItemListener(this);
		ascalecheck=new Checkbox("Autoscale",(ascale==1) ? true : false);
		ascalecheck.setBounds(150,10+imageheight+50+256+50,80,20);
		ascalecheck.addItemListener(this);

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

		add(xminval); add(xmaxval); add(yminval); add(ymaxval);
		add(scalehistval); add(dispminval); add(dispminlabel); add(threshlabel); add(threshval); add(threshfirstcheck);
		add(dispmaxval); add(dispmaxlabel); add(scalehistlabel);
		add(roisquare); add(roioval); add(roiwidthval); add(roiheightval); add(roixval); add(roiyval);
		add(roiwidthlabel); add(roiheightlabel); add(smooth_button); add(revert_button);
		add(sval); add(s0val); add(offval); add(plotnorm); add(plotdiff); add(medsmooth); add(gassmooth);  add(binsmooth);
		add(roixlabel); add(roiylabel); add(slabel); add(s0label); add(offlabel);
		add(saveimg_button); add(savehist_button); add(scorelabel); add(xlabel); add(ylabel); add(saveyimg_button);
		add(backlabel); add(backval); add(ascalecheck); add(slicelabel); add(sliceval); add(slice_slider);
		add(xavglabel); add(yavglabel); add(xstdevlabel); add(ystdevlabel);
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
		if(roishape==0){g.drawRect(100+roix,(50+imageheight+10+256)-(roiy+roiheight),roiwidth,roiheight);}
		else{g.drawOval(100+roix,(50+imageheight+10+256)-(roiy+roiheight),roiwidth,roiheight);}
		g.setColor(Color.green);
		float intercept = s0-off*s;
		float y1,y2;
		if(plottype==0){
			y1=s*(xmin+off)+intercept-s0;
			y2=s*(xmax+off)+intercept-s0;
		} else {
			y1=1.0f;
			y2=1.0f;
		}
		g.setColor(Color.green);
		g.drawLine(100,imageheight+50+10+256-(int)(((y1-ymin)/(ymax-ymin))*256.0),100+256,imageheight+50+10+256-(int)(((y2-ymin)/(ymax-ymin))*256.0));
		g.setClip(null);
	}

	public void update(Graphics g)
	{
		paint(g);
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==roisquare || e.getSource()==roioval){
			if(roisquare.getState()==true){roishape=0;}
			else{roishape=1;}
		}
		if(e.getSource()==plotnorm || e.getSource()==plotdiff){
			if(plotnorm.getState()==true){plottype=0;}
			else{plottype=1;}
			histimg=create_2D_histogram(); 
			update_mask();
			dispimg=create_masked_image();
			xlabel.setText("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		}
		if(e.getSource()==medsmooth || e.getSource()==gassmooth || e.getSource()==binsmooth){
			if(medsmooth.getState()){smoothtype=0;}
			else{
				if(gassmooth.getState()){smoothtype=1;}
				else{smoothtype=2;}
			}
		}
		if(e.getSource()==ascalecheck){
			if(ascalecheck.getState()==true){
				ascale=1;
				histimg=create_2D_histogram();
				update_mask();
				dispimg=create_masked_image();
				xlabel.setText("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
				ylabel.setText("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
			}
			else{ascale=0;}
		}
		if(e.getSource()==threshfirstcheck){
			if(threshfirstcheck.getState()==true){
				threshfirst=true;
			}
			else{
				threshfirst=false;
			}
			histimg=create_2D_histogram();
			update_mask();
			dispimg=create_masked_image();
			xlabel.setText("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		}
		repaint();
		set_options();
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==smooth_button){smooth();}
		if(e.getSource()==revert_button){revert();}
		if(e.getSource()==savehist_button){save_histogram();}
		if(e.getSource()==saveimg_button){save_masked_image();}
		if(e.getSource()==saveyimg_button){save_histypix();}
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
		xlabel.setText("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
		ylabel.setText("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
		s=Float.parseFloat(sval.getText());
		s0=Float.parseFloat(s0val.getText());
		off=Float.parseFloat(offval.getText());
		back=Float.parseFloat(backval.getText());
		thresh=Float.parseFloat(threshval.getText());
		currslice=(int)Float.parseFloat(sliceval.getText());
		if(currslice>=slices){currslice=slices-1;}
		sliceval.setText(""+currslice);
		slice_slider.setValue(currslice);
		histimg=create_2D_histogram(); 
		update_mask();
		dispimg=create_masked_image();
		repaint();
		set_options();
	}

	public void adjustmentValueChanged(AdjustmentEvent e)
	{ 
		currslice=slice_slider.getValue();
		sliceval.setText(""+currslice);
		histimg=create_2D_histogram(); 
		update_mask();
		dispimg=create_masked_image();
		repaint();
	}

	public void mouseClicked(MouseEvent e){}

	public void mouseEntered(MouseEvent e){}

	public void mouseExited(MouseEvent e){}

	public void mousePressed(MouseEvent e){
		int xpos=e.getX();
		int xvalue=xpos-100;
		int ypos=e.getY();
		int yvalue=(256+imageheight+10+50)-ypos;
		if(roishape==0){
			//roi is rectangular
			if(xvalue>roix && xvalue<(roix+roiwidth) && yvalue>roiy && yvalue<(roiy+roiheight)){
				inroi=true;
		}}
		if(roishape==1){
			//roi is oval
			//calculate the center of the ellipse and the positions of the foci
			float centerx=(float)roix+0.5f*(float)roiwidth;
			float centery=(float)roiy+0.5f*(float)roiheight;
			float f1x,f1y,f2x,f2y,major;
			if(roiwidth>=roiheight){
				f1y=centery; f2y=centery;
				f1x=centerx+(float)Math.sqrt((float)(roiwidth*roiwidth-roiheight*roiheight)/4.0f);
				f2x=centerx-(float)Math.sqrt((float)(roiwidth*roiwidth-roiheight*roiheight)/4.0f);
				major=(float)roiwidth;
			}
			else{
				f1x=centerx; f2x=centerx;
				f1y=centery+(float)Math.sqrt((float)(roiheight*roiheight-roiwidth*roiwidth)/4.0f);
				f2y=centery-(float)Math.sqrt((float)(roiheight*roiheight-roiwidth*roiwidth)/4.0f);
				major=(float)roiheight;
			}
			float distance1,distance2;
			distance1=(float)Math.sqrt(((float)xvalue-f1x)*((float)xvalue-f1x)+((float)yvalue-f1y)*((float)yvalue-f1y));
			distance2=(float)Math.sqrt(((float)xvalue-f2x)*((float)xvalue-f2x)+((float)yvalue-f2y)*((float)yvalue-f2y));
			//IJ.log("f1x "+f1x+"f2x "+f2x);
			if((distance1+distance2)<=major){
				inroi=true;
		}}
		oldxpos=xpos;
		oldypos=ypos;
	}

	public void mouseReleased(MouseEvent e){
		inroi=false;
		set_options();
	}

	public void mouseMoved(MouseEvent e){}

	public void mouseDragged(MouseEvent e){
		if(inroi){
			int xpos=e.getX();
			int ypos=e.getY();
			roix+=xpos-oldxpos;
			roiy-=(ypos-oldypos);
			roixval.setText(""+(roix+roiwidth/2));
			roiyval.setText(""+(roiy+roiheight/2));
			xlabel.setText("x: "+((((float)roix+(float)roiwidth/2.0f)/256.0f)*(xmax-xmin)+xmin));
			ylabel.setText("y: "+((((float)roiy+(float)roiheight/2.0f)/256.0f)*(ymax-ymin)+ymin));
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
			if(farr[i]<min){min=farr[i];}
		}
		return min;
	}
	
	float getmaxgreatzero(float[] farr){
		float max;
		max=0.0f;
		for(int i=0;i<farr.length;i++){
			if(farr[i]>max){max=farr[i];}
		}
		return max;
	}
	
	Image create_masked_image(){
		int[] color_pixels = new int[width*height];
		int dumint;
		for(int i=0;i<width*height;i++){
			dumint=(int)(((disppix[i+currslice*width*height]-dispmin)/(dispmax-dispmin))*256.0);
			if(threshfirst){
				if(disppix[i]<thresh){dumint=0;}
			} else {
				if(disppix[i+currslice*width*height]<thresh){dumint=0;}
			}
			if(dumint<0){dumint=0;}
			if(dumint>255){dumint=255;}
			if(disppix[i+currslice*width*height]>=dispmin){
				if(mask[i]==0){color_pixels[i]=0xff000000 + (dumint<<16) + (dumint<<8) + dumint;}
				else{color_pixels[i]=0xffff0000;}
			}
			else{color_pixels[i]=0xff000000;}
		}
		//if the image is larger than 256, scale it down
		int scalefactor=1+(int)((float)(width-1)/256.0f);
		dumint=1+(int)((float)(height-1)/256.0f);
		if(dumint>scalefactor){scalefactor=dumint;}
		int newwidth=width/scalefactor;
		int newheight=height/scalefactor;
		if(scalefactor>1){
			int[] newcolorpixels=new int[newwidth*newheight];
			for(int i=0;i<newheight;i++){
				for(int j=0;j<newwidth;j++){
					boolean masked=false;
					float avg=0.0f;
					local_search:
					for(int k=0;k<scalefactor;k++){
						for(int l=0;l<scalefactor;l++){
							dumint=j*scalefactor+l+(i*scalefactor+k)*width;
							if(mask[dumint]==1){masked=true; break local_search;}
							avg+=(float)(color_pixels[dumint]&0x000000ff)/(float)(scalefactor*scalefactor);
						}
					}
					dumint=(int)avg;
					if(masked){newcolorpixels[i*newwidth+j]=0xffff0000;}
					else{newcolorpixels[i*newwidth+j]=0xff000000 + (dumint<<16) + (dumint<<8) + dumint;}
				}
			}
			color_pixels=newcolorpixels;
		}
		ColorProcessor cp = new ColorProcessor(newwidth,newheight,color_pixels);
		return cp.createImage();
	}

	void update_mask(){
		int dumint1,dumint2;
		//clear the mask image
		mask = new int[width*height];
		for(int i=0;i<width*height;i++){mask[i]=0;}
		//calculate the mask based on histogram hit testing
		score=0.0f;
		xavg=0.0f;
		yavg=0.0f;
		xstdev=0.0f;
		ystdev=0.0f;
		float totpixels=0.0f;
		for(int i=0;i<width*height;i++){
			dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
			dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
			boolean abovethresh=true;
			if(threshfirst){
				if(disppix[i]<thresh){abovethresh=false;}
			} else {
				if(disppix[i+currslice*width*height]<thresh){abovethresh=false;}
			}
			if(((dumint1<256 && dumint2<256) && (dumint1>0 && dumint2>0)) && abovethresh){
				if(roishape==0){
					//roi is rectangular
					if((dumint1>roix && dumint1<(roix+roiwidth)) && (dumint2>roiy && dumint2<(roiy+roiheight))){
						mask[i]=1;
						score+=1.0;
						xavg+=histxpix[i];
						yavg+=histypix[i];
						xstdev+=histxpix[i]*histxpix[i];
						ystdev+=histypix[i]*histypix[i];
				}}
				if(roishape==1){
					//roi is oval
					//calculate the center of the ellipse and the positions of the foci
					float centerx=(float)roix+0.5f*(float)roiwidth;
					float centery=(float)roiy+0.5f*(float)roiheight;
					float f1x,f1y,f2x,f2y,major;
					if(roiwidth>=roiheight){
						f1y=centery; f2y=centery;
						f1x=centerx+(float)Math.sqrt((float)(roiwidth*roiwidth-roiheight*roiheight)/4.0f);
						f2x=centerx-(float)Math.sqrt((float)(roiwidth*roiwidth-roiheight*roiheight)/4.0f);
						major=(float)roiwidth;
					}
					else{
						f1x=centerx; f2x=centerx;
						f1y=centery+(float)Math.sqrt((float)(roiheight*roiheight-roiwidth*roiwidth)/4.0f);
						f2y=centery-(float)Math.sqrt((float)(roiheight*roiheight-roiwidth*roiwidth)/4.0f);
						major=(float)roiwidth;
					}
					float distance1,distance2;
					distance1=(float)Math.sqrt(((float)dumint1-f1x)*((float)dumint1-f1x)+((float)dumint2-f1y)*((float)dumint2-f1y));
					distance2=(float)Math.sqrt(((float)dumint1-f2x)*((float)dumint1-f2x)+((float)dumint2-f2y)*((float)dumint2-f2y));
					if((distance1+distance2)<=major){
						mask[i]=1;
						score+=1.0;
						xavg+=histxpix[i];
						yavg+=histypix[i];
						xstdev+=histxpix[i]*histxpix[i];
						ystdev+=histypix[i]*histypix[i];
				}}
			}
			if(abovethresh){totpixels+=1.0f;}
		}
		xavg/=score;
		yavg/=score;
		xstdev/=score;
		ystdev/=score;
		xstdev-=xavg*xavg;
		ystdev-=yavg*yavg;
		xstdev*=(1.0f-1.0f/score);
		ystdev*=(1.0f-1.0f/score);
		xstdev=(float)Math.sqrt(xstdev);
		ystdev=(float)Math.sqrt(ystdev);
		score/=totpixels;
		scorelabel.setText("Score = "+score);
		xavglabel.setText("x avg = "+xavg);
		yavglabel.setText("y avg = "+yavg);
		xstdevlabel.setText("x stdev = "+xstdev);
		ystdevlabel.setText("y stdev = "+ystdev);
	}

	Image create_2D_histogram(){
		int dumint1,dumint2;
		update_histpix();
		int[][] histogram=new int[256][256];
		//calculate the histogram
		for(int i=0;i<width*height;i++){
			dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
			dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
			boolean abovethresh=true;
			if(threshfirst){
				if(disppix[i]<thresh){abovethresh=false;}
			} else{
				if(disppix[i+currslice*width*height]<thresh){abovethresh=false;}
			}
			if(abovethresh){
				if((dumint1<256 && dumint2<256) && (dumint1>0 && dumint2>0)){
					histogram[dumint2][dumint1]++;}
			}
		}
		//transfer it to a single array for output
		float[] tempfloat = new float[256*256];
		histmax=0;
		for(int i=0;i<256;i++){
			for(int j=0;j<256;j++){
				tempfloat[(255-i)*256+j]=(float)histogram[i][j];
				if(histogram[i][j]>histmax){histmax=histogram[i][j];}
			}
		}
		FloatProcessor fp = new FloatProcessor(256,256,tempfloat,null);
		fp.setMinAndMax(0.0,(double)(((float)histmax)/multiplier));
		update_mask();
		return fp.createImage();
	}

	void smooth(){
		if(smoothtype!=2){
			//here we use a median or gaussian smoothing method
			//the weighting is 2 for the pixel and 1 for its neighbors
			//the pixels around the edge are not smoothed
			int i,j,k,l;
			float dumflt;
			float[] pixelvals=new float[10];
			//first smooth the x image
			for(int m=0;m<slices;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothxpix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothxpix[i*width+j+m*width*height];
						}
						if(smoothtype==0){smoothxpix[i*width+j+m*width*height]=median(pixelvals);}
						else{smoothxpix[i*width+j+m*width*height]=avg(pixelvals);}
					}
				}
			}
			//now smooth the y image
			for(int m=0; m<slices;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothypix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothypix[i*width+j+m*width*height];
						}
						if(smoothtype==0){smoothypix[i*width+j+m*width*height]=median(pixelvals);}
						else{smoothypix[i*width+j+m*width*height]=avg(pixelvals);}
					}
				}
			}
		} else {
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
							int npix=0;
							for(int k=0;k<binsize;k++){
								for(int l=0;l<binsize;l++){
									boolean abovethresh=true;
									if(threshfirst){
										if(disppix[j*binsize+l+(i*binsize+k)*width]<thresh){abovethresh=false;}
									} else {
										if(disppix[j*binsize+l+(i*binsize+k)*width+m*width*height]<thresh){abovethresh=false;}
									}
									if(abovethresh){
										npix++;
										sumx+=smoothxpix[j*binsize+l+(i*binsize+k)*width+m*width*height];
										sumy+=smoothypix[j*binsize+l+(i*binsize+k)*width+m*width*height];
									}
								}
							}
							if(npix>0){
								sumx/=(float)npix;
								sumy/=(float)npix;
							}
							for(int k=0;k<binsize;k++){
								for(int l=0;l<binsize;l++){
									smoothxpix[j*binsize+l+(i*binsize+k)*width+m*width*height]=sumx;
									smoothypix[j*binsize+l+(i*binsize+k)*width+m*width*height]=sumy;
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
		smoothxpix=new float[width*height*slices];
		smoothypix=new float[width*height*slices];
		for(i=0;i<width*height*slices;i++){
			smoothxpix[i]=xpix[i];
			smoothypix[i]=ypix[i];
		}
	}

	float median(float[] values){
		int i,length,dumint,imin,j;
		float min,dumfloat;
		length=values.length;
		dumint=1+(int)((float)length/2.0f);
		//sort the vector one at a time
		for(j=0;j<dumint;j++){
			imin=j;
			min=values[imin];
			for(i=imin+1;i<length;i++){
				if(values[i]<min){imin=i; min=values[imin];}
			}
			dumfloat=values[j]; values[j]=values[imin]; values[imin]=dumfloat;
		}
		if((dumint*2)>length){return values[dumint-1];}
		else{
			min=values[dumint];
			for(i=dumint+1;i<length;i++){
				if(values[i]<min){min=values[i];}
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
		gdsave.showDialog(); if(gdsave.wasCanceled()){return;}
		saveashist=gdsave.getNextBoolean();
		if(!saveashist){
			int ascaled,oldslice;
			ImageStack out_stack=new ImageStack(256,256);
			if(ascale==1){ascaled=1;}
			else{ascaled=0;}
			ascale=0; oldslice=currslice;
			for(int k=0;k<slices;k++){
				currslice=k;
				update_histpix();
				int[][] histogram=new int[256][256];
				//calculate the histogram
				for(int i=0;i<width*height;i++){
					int dumint1=(int)(((histxpix[i]-xmin)/(xmax-xmin))*256.0);
					int dumint2=(int)(((histypix[i]-ymin)/(ymax-ymin))*256.0);
					boolean abovethresh=true;
					if(threshfirst){
						if(disppix[i]<thresh){abovethresh=false;}
					} else{
						if(disppix[i+currslice*width*height]<thresh){abovethresh=false;}
					}
					if(abovethresh){
						if((dumint1<256 && dumint2<256) && (dumint1>0 && dumint2>0)){
							histogram[dumint2][dumint1]++;}
					}
				}
				//transfer it to a single array for output
				float[] tempfloat = new float[256*256];
				for(int i=0;i<256;i++){
					for(int j=0;j<256;j++){
						tempfloat[(255-i)*256+j]=(float)histogram[i][j];
					}
				}
				out_stack.addSlice("",(Object)tempfloat);
			}
			if(ascaled==1){ascale=1;}
			currslice=oldslice;
			ImagePlus imp = new ImagePlus("2D histogram",out_stack);
			Calibration cal=imp.getCalibration();
			cal.pixelWidth=(double)(xmax-xmin)/256.0f;
			cal.pixelHeight=(double)(ymax-ymin)/256.0f;
			cal.xOrigin=-(double)xmin/cal.pixelWidth;
			cal.yOrigin=(double)ymin/cal.pixelHeight+255.0;
			cal.setInvertY(true);
			//IJ.log("ymax = "+ymax);
			cal.setUnit("nb");
			imp.setRoi(new Roi(roix,256-roiy-roiheight,roiwidth,roiheight));
			imp.show();
		} else {
			int ascaled,oldslice;
			if(ascale==1){ascaled=1;}
			else{ascaled=0;}
			ascale=0; oldslice=currslice;
			float[] outhistxpix=new float[slices*width*height];
			float[] outhistypix=new float[slices*width*height];
			int counter=0;
			for(int k=0;k<slices;k++){
				currslice=k;
				update_histpix();
				//calculate the histogram
				for(int i=0;i<width*height;i++){
					boolean abovethresh=true;
					if(threshfirst){
						if(disppix[i]<thresh){abovethresh=false;}
					} else{
						if(disppix[i+currslice*width*height]<thresh){abovethresh=false;}
					}
					if(abovethresh){
						outhistxpix[counter]=histxpix[i];
						outhistypix[counter]=histypix[i];
						counter++;
					}
				}
			}
			if(ascaled==1){ascale=1;}
			currslice=oldslice;
			float[] outhistxpix2=new float[counter];
			float[] outhistypix2=new float[counter];
			System.arraycopy(outhistxpix,0,outhistxpix2,0,counter);
			System.arraycopy(outhistypix,0,outhistypix2,0,counter);
			PlotWindow2DHist pw=new PlotWindow2DHist("2D N B Histogram","I/S","B/S",outhistxpix2,outhistypix2,null);
			float[] limits=pw.getLimits();
			limits[0]=xmin; limits[1]=xmax; limits[2]=ymin; limits[3]=ymax;
			pw.setLimits(limits);
			pw.intautoscale();
			pw.draw();
		}
	}

	void update_histpix(){
		histxpix=new float[width*height];
		histypix=new float[width*height];
		if(plottype==0){
			for(int i=0;i<width*height;i++){
				histxpix[i]=smoothxpix[i+currslice*width*height]-off;
				histypix[i]=smoothypix[i+currslice*width*height]-s0;
			}
		} else {
			for(int i=0;i<width*height;i++){
				histxpix[i]=(smoothxpix[i+currslice*width*height]-off-back)/s;
				histypix[i]=(smoothypix[i+currslice*width*height]-s0-back*s)/(histxpix[i]*s*s);
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
		int[] color_pixels = new int[width*height];
		int dumint;
		for(int i=0;i<width*height;i++){
			dumint=(int)(((disppix[i+currslice*width*height]-dispmin)/(dispmax-dispmin))*256.0);
			if(threshfirst){
				if(disppix[i]<thresh){dumint=0;}
			} else {
				if(disppix[i+currslice*width*height]<thresh){dumint=0;}
			}
			if(dumint<0){dumint=0;}
			if(dumint>255){dumint=255;}
			if(disppix[i+currslice*width*height]>=dispmin){
				if(mask[i]==0){color_pixels[i]=0xff000000 + (dumint<<16) + (dumint<<8) + dumint;}
				else{color_pixels[i]=0xffff0000;}
			}
			else{color_pixels[i]=0xff000000;}
		}
		ColorProcessor cp = new ColorProcessor(width,height,color_pixels);
		ImagePlus imp2 = new ImagePlus("Masked Image",cp);
		imp2.setCalibration(cal);
		imp2.show();
	}

	void save_histypix(){
		int tempcurrslice=currslice;
		ImageStack histystack=new ImageStack(width,height);
		for(int j=0;j<slices;j++){
			float[] filteredypix=new float[width*height];
			currslice=j;
			update_histpix();
			for(int i=0;i<width*height;i++){
				boolean abovethresh=true;
				if(threshfirst){
					if(disppix[i]<thresh){abovethresh=false;}
				} else{
					if(disppix[i+currslice*width*height]<thresh){abovethresh=false;}
				}
				if(abovethresh){
					if((histxpix[i]<xmax && histxpix[i]>xmin) && (histypix[i]<ymax && histypix[i]>ymin)){
						filteredypix[i]=histypix[i];}
					else{filteredypix[i]=0.0f;}
				}
			}
			histystack.addSlice("",(Object)filteredypix);
		}
		ImagePlus imphisty=new ImagePlus("Y Image",histystack);
		imphisty.setCalibration(cal);
		imphisty.show();
		currslice=tempcurrslice;
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
			threshfirst = (Integer.parseInt(d.readLine())==0) ? false : true;
			xmin=Float.parseFloat(d.readLine());
			xmax=Float.parseFloat(d.readLine());
			ymin=Float.parseFloat(d.readLine());
			ymax=Float.parseFloat(d.readLine());
			s=Float.parseFloat(d.readLine());
			s0=Float.parseFloat(d.readLine());
			off=Float.parseFloat(d.readLine());
			back=Float.parseFloat(d.readLine());
			d.close();
		}
		catch(IOException e){
			roix=0;
			roiy=0;
			roiwidth=30;
			roiheight=30;
			plottype=0;
			smoothtype=1;
			roishape=0;
			ascale=1;
			threshfirst = true;
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
			if(!a.exists()){a.mkdir();}
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
			d.write(""+(threshfirst ? 1:0)+"\n");
			d.write(""+xmin+"\n");
			d.write(""+xmax+"\n");
			d.write(""+ymin+"\n");
			d.write(""+ymax+"\n");
			d.write(""+s+"\n");
			d.write(""+s0+"\n");
			d.write(""+off+"\n");
			d.write(""+back+"\n");
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

}
