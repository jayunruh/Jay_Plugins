/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import jalgs.jstatistics;
import jalgs.jfit.PALMutils;

import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
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

public class PALMWindow extends Panel implements ActionListener,AdjustmentListener,MouseListener,MouseMotionListener,ItemListener{

	public final static int H=612;
	public final static int WR=800;
	public Dimension totalSize=new Dimension();
	public int mw=4;
	public float maxerr;
	private TextField xminval,xmaxval,yminval,ymaxval,c2minval,c2maxval,ampminval,ampmaxval;
	private TextField dispminval,dispmaxval,gainval,resmultval;
	private Label activitylabel,c2minlabel,c2maxlabel,ampminlabel,ampmaxlabel,dispminlabel;
	private Label gainlabel,dispmaxlabel,resmultlabel;
	private Button rescalezoombutton,rescaleintensitybutton,saveimagebutton;
	private Button nudgeleftbutton,nudgerightbutton,nudgedownbutton,nudgeupbutton;
	private Button scaleupbutton,scaledownbutton;
	Checkbox showthreshcheck,addnoisecheck,showorigcheck;
	Image dispimg;
	Graphics tempgraphics;
	float xmin,xmax,ymin,ymax,dispmin,dispmax,c2min,c2max,ampmin,ampmax,baseline,stdev;
	float scaling,gain,resmult;
	float[][] molecules;
	int width,height,nummolecules,startx,starty,roiwidth;
	int imageheight=512;
	float[] fdispimg,fthreshimg;
	boolean showthresh,inroi,addnoise,showorig;

	public static void launch_frame(PALMWindow panel){
		final Frame f=new Frame("Interactive PALM Plot");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				f.dispose();
			}
		});
		f.add(panel);
		f.pack();
		f.setResizable(false);
		Insets ins=f.getInsets();
		panel.totalSize.height=PALMWindow.H+ins.bottom+ins.top+65;
		panel.totalSize.width=PALMWindow.WR+ins.left+ins.right;
		// panel.setBounds(ins.top+5,ins.left+5,panel.totalSize.width,panel.totalSize.height);
		f.setSize(panel.totalSize);
		f.setVisible(true);
		panel.requestFocus();
	}

	public void init(float[] molecules1){
		setLayout(null);
		width=(int)molecules1[0];
		height=width;
		baseline=molecules1[1];
		stdev=molecules1[2];
		maxerr=2.0f*stdev;
		gain=molecules1[3];
		resmult=1.0f;
		int temp=molecules1.length/mw;
		nummolecules=temp-1;
		molecules=new float[mw][nummolecules];
		// the molecules array has x,y,amp,c2 and maybe other parameters
		for(int i=1;i<temp;i++){
			for(int j=0;j<mw;j++)
				molecules[j][i-1]=molecules1[i*mw+j];
		}

		// initialize all of the variables
		xmin=0.0f;
		xmax=width-1;
		ymin=0.0f;
		ymax=height-1;
		scaling=512.0f/width;

		c2min=jstatistics.getstatistic("Min",molecules[3],null);
		c2max=jstatistics.getstatistic("Max",molecules[3],null);
		ampmin=jstatistics.getstatistic("Min",molecules[2],null);
		ampmax=jstatistics.getstatistic("Max",molecules[2],null);
		ampmin=get_gaus_integral(ampmin,stdev)/gain;
		ampmax=get_gaus_integral(ampmax,stdev)/gain;

		showthresh=false;
		addnoise=true;
		inroi=false;

		addMouseMotionListener(this);
		addMouseListener(this);

		xminval=new TextField(""+xmin,10);
		xminval.setBounds(10,10+imageheight+10,80,20);
		xminval.addActionListener(this);
		xmaxval=new TextField(""+xmax,10);
		xmaxval.setBounds(10+imageheight-40,10+imageheight+10,80,20);
		xmaxval.addActionListener(this);
		yminval=new TextField(""+ymin,10);
		yminval.addActionListener(this);
		yminval.setBounds(10+imageheight+10,10,80,20);
		ymaxval=new TextField(""+ymax,10);
		ymaxval.setBounds(10+imageheight+10,imageheight-10,80,20);
		ymaxval.addActionListener(this);

		c2minval=new TextField(""+c2min,10);
		c2minval.setBounds(10+imageheight+100,10,80,20);
		c2minval.addActionListener(this);
		c2minlabel=new Label("chi^2 min");
		c2minlabel.setBounds(10+imageheight+100+100,10,80,20);
		add(c2minlabel);
		c2maxval=new TextField(""+c2max,10);
		c2maxval.setBounds(10+imageheight+100,40,80,20);
		c2maxval.addActionListener(this);
		c2maxlabel=new Label("chi^2 max");
		c2maxlabel.setBounds(10+imageheight+100+100,40,80,20);
		add(c2maxlabel);
		ampminval=new TextField(""+ampmin,10);
		ampminval.setBounds(10+imageheight+100,70,80,20);
		ampminval.addActionListener(this);
		ampminlabel=new Label("bright min");
		ampminlabel.setBounds(10+imageheight+100+100,70,80,20);
		add(ampminlabel);
		ampmaxval=new TextField(""+ampmax,10);
		ampmaxval.setBounds(10+imageheight+100,100,80,20);
		ampmaxval.addActionListener(this);
		ampmaxlabel=new Label("bright max");
		ampmaxlabel.setBounds(10+imageheight+100+100,100,80,20);
		add(ampmaxlabel);

		activitylabel=new Label("");
		activitylabel.setBounds(10+imageheight+100,220,80,20);

		generate_fimage();
		dispmin=jstatistics.getstatistic("Min",fdispimg,null);
		dispmax=jstatistics.getstatistic("Max",fdispimg,null);

		dispminval=new TextField(""+dispmin,10);
		dispminval.setBounds(10+imageheight+100,130,80,20);
		dispminval.addActionListener(this);
		dispminlabel=new Label("intensity min");
		dispminlabel.setBounds(10+imageheight+100+100,130,80,20);
		add(dispminlabel);
		dispmaxval=new TextField(""+dispmax,10);
		dispmaxval.setBounds(10+imageheight+100,160,80,20);
		dispmaxval.addActionListener(this);
		dispmaxlabel=new Label("intensity max");
		dispmaxlabel.setBounds(10+imageheight+100+100,160,80,20);
		add(dispmaxlabel);

		showthreshcheck=new Checkbox("Show Threshhold",showthresh);
		showthreshcheck.setBounds(10+imageheight+100,190,150,20);
		showthreshcheck.addItemListener(this);

		// activity label is at y position 220

		rescalezoombutton=new Button("Rescale Zoom");
		rescalezoombutton.setBounds(10+imageheight+100,250,120,40);
		rescalezoombutton.addActionListener(this);
		add(rescalezoombutton);

		rescaleintensitybutton=new Button("Rescale Intensity");
		rescaleintensitybutton.setBounds(10+imageheight+100,300,120,40);
		rescaleintensitybutton.addActionListener(this);
		add(rescaleintensitybutton);

		addnoisecheck=new Checkbox("Add Noise",addnoise);
		addnoisecheck.setBounds(10+imageheight+100,350,150,20);
		addnoisecheck.addItemListener(this);

		saveimagebutton=new Button("Save Image");
		saveimagebutton.setBounds(10+imageheight+100,380,120,40);
		saveimagebutton.addActionListener(this);
		add(saveimagebutton);

		gainlabel=new Label("Gain");
		gainlabel.setBounds(10+imageheight+100,430,80,20);
		add(gainlabel);
		gainval=new TextField(""+gain,10);
		gainval.setBounds(10+imageheight+100+100,430,60,20);
		gainval.addActionListener(this);
		add(gainval);

		showorig=false;
		showorigcheck=new Checkbox("Show Original");
		showorigcheck.setBounds(10+imageheight+100,460,100,20);
		showorigcheck.addItemListener(this);
		add(showorigcheck);

		resmultlabel=new Label("Res Multiplier");
		resmultlabel.setBounds(10+imageheight+100,490,80,20);
		add(resmultlabel);
		resmultval=new TextField(""+resmult,10);
		resmultval.setBounds(10+imageheight+100+100,490,60,20);
		resmultval.addActionListener(this);
		add(resmultval);

		nudgerightbutton=new Button(">");
		nudgerightbutton.setBounds(10+imageheight+5,250,40,40);
		nudgerightbutton.addActionListener(this);
		add(nudgerightbutton);

		nudgeleftbutton=new Button("<");
		nudgeleftbutton.setBounds(10+imageheight+5,200,40,40);
		nudgeleftbutton.addActionListener(this);
		add(nudgeleftbutton);

		nudgeupbutton=new Button("v");
		nudgeupbutton.setBounds(250,10+imageheight+5,40,40);
		nudgeupbutton.addActionListener(this);
		add(nudgeupbutton);

		nudgedownbutton=new Button("^");
		nudgedownbutton.setBounds(200,10+imageheight+5,40,40);
		nudgedownbutton.addActionListener(this);
		add(nudgedownbutton);

		scaleupbutton=new Button("ZIn");
		scaleupbutton.setBounds(200,10+imageheight+50,40,40);
		scaleupbutton.addActionListener(this);
		add(scaleupbutton);

		scaledownbutton=new Button("ZOut");
		scaledownbutton.setBounds(250,10+imageheight+50,40,40);
		scaledownbutton.addActionListener(this);
		add(scaledownbutton);

		add(xminval);
		add(xmaxval);
		add(yminval);
		add(ymaxval);
		add(c2minval);
		add(c2maxval);
		add(ampminval);
		add(ampmaxval);
		add(dispminval);
		add(dispmaxval);
		add(showthreshcheck);
		add(activitylabel);
		add(addnoisecheck);

		generate_image();
		repaint();
		// set_options();
	}

	public void paint(Graphics g){
		g.drawImage(dispimg,10,10,this);
		if(inroi){
			g.setColor(Color.blue);
			g.drawRect(startx,starty,roiwidth,roiwidth);
		}
	}

	public void update(Graphics g){
		paint(g);
	}

	public void itemStateChanged(ItemEvent e){
		showthresh=showthreshcheck.getState();
		addnoise=addnoisecheck.getState();
		showorig=showorigcheck.getState();
		generate_fimage();
		generate_image();
		repaint();
	}

	public void actionPerformed(ActionEvent e){
		float xminold=xmin;
		float xmaxold=xmax;
		float yminold=ymin;
		float ymaxold=ymax;
		float c2minold=c2min;
		float c2maxold=c2max;
		float ampminold=ampmin;
		float ampmaxold=ampmax;
		float gainold=gain;
		float resmultold=resmult;
		boolean updateimg=false;
		boolean updatescale=false;
		xmin=Float.parseFloat(xminval.getText());
		if(xmin!=xminold){
			updateimg=true;
		}
		xmax=Float.parseFloat(xmaxval.getText());
		if(xmax!=xmaxold){
			updateimg=true;
		}
		ymin=Float.parseFloat(yminval.getText());
		if(ymin!=yminold){
			updateimg=true;
		}
		ymax=Float.parseFloat(ymaxval.getText());
		if(ymax!=ymaxold){
			updateimg=true;
		}
		c2min=Float.parseFloat(c2minval.getText());
		if(c2min!=c2minold){
			updateimg=true;
			updatescale=true;
		}
		c2max=Float.parseFloat(c2maxval.getText());
		if(c2max!=c2maxold){
			updateimg=true;
			updatescale=true;
		}
		ampmin=Float.parseFloat(ampminval.getText());
		if(ampmin!=ampminold){
			updateimg=true;
			updatescale=true;
		}
		ampmax=Float.parseFloat(ampmaxval.getText());
		if(ampmax!=ampmaxold){
			updateimg=true;
			updatescale=true;
		}
		gain=Float.parseFloat(gainval.getText());
		if(gain!=gainold){
			updateimg=true;
			updatescale=true;
		}
		resmult=Float.parseFloat(resmultval.getText());
		if(resmult!=resmultold){
			updateimg=true;
			updatescale=true;
		}
		if(e.getSource()==rescalezoombutton){
			xmin=0.0f;
			xmax=width-1;
			ymin=0.0f;
			ymax=height-1;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			updateimg=true;
			// updatescale=true;
		}
		if(e.getSource()==saveimagebutton){
			GenericDialog gd=new GenericDialog("Options");
			gd.addCheckbox("Color_Output?",false);
			gd.showDialog();
			if(gd.wasCanceled()){
				return;
			}
			boolean color=gd.getNextBoolean();
			ImagePlus imp;
			if(color)
				imp=new ImagePlus("PALM Image",new ColorProcessor(dispimg));
			else
				imp=new ImagePlus("PALM Image",new FloatProcessor(512,512,fdispimg,null));
			jutils.set_psize(imp,(xmax-xmin)/512.0,"units");
			imp.show();
		}
		if(e.getSource()==nudgedownbutton){
			scaling=512.0f/(xmax-xmin);
			ymax+=32.0f/scaling;
			ymin+=32.0f/scaling;
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			updateimg=true;
		}
		if(e.getSource()==nudgeupbutton){
			scaling=512.0f/(xmax-xmin);
			ymax-=32.0f/scaling;
			ymin-=32.0f/scaling;
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			updateimg=true;
		}
		if(e.getSource()==nudgeleftbutton){
			scaling=512.0f/(xmax-xmin);
			xmax+=32.0f/scaling;
			xmin+=32.0f/scaling;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			updateimg=true;
		}
		if(e.getSource()==nudgerightbutton){
			scaling=512.0f/(xmax-xmin);
			xmax-=32.0f/scaling;
			xmin-=32.0f/scaling;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			updateimg=true;
		}
		if(e.getSource()==scaleupbutton){
			scaling=512.0f/(xmax-xmin);
			float centerx=0.5f*(xmax+xmin);
			float centery=0.5f*(ymax+ymin);
			scaling*=1.5f;
			xmax=centerx+256.0f/scaling;
			xmin=xmax-512.0f/scaling;
			ymax=centery+256.0f/scaling;
			ymin=ymax-512.0f/scaling;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			updateimg=true;
		}
		if(e.getSource()==scaledownbutton){
			scaling=512.0f/(xmax-xmin);
			float centerx=0.5f*(xmax+xmin);
			float centery=0.5f*(ymax+ymin);
			scaling/=1.5f;
			xmax=centerx+256.0f/scaling;
			xmin=xmax-512.0f/scaling;
			ymax=centery+256.0f/scaling;
			ymin=ymax-512.0f/scaling;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			updateimg=true;
		}
		if((xmax-xmin)>(ymax-ymin)){
			scaling=512.0f/(xmax-xmin);
			ymax=ymin+512.0f/scaling;
			ymaxval.setText(""+ymax);
		}else{
			scaling=512.0f/(ymax-ymin);
			xmax=xmin+512.0f/scaling;
			xmaxval.setText(""+xmax);
		}
		dispmin=Float.parseFloat(dispminval.getText());
		dispmax=Float.parseFloat(dispmaxval.getText());
		if(e.getSource()==rescaleintensitybutton){
			updatescale=true;
		}
		if(updateimg){
			generate_fimage();
		}
		if(updatescale){
			dispmin=jstatistics.getstatistic("Min",fdispimg,null);
			dispmax=jstatistics.getstatistic("Max",fdispimg,null);
			dispminval.setText(""+dispmin);
			dispmaxval.setText(""+dispmax);
		}
		generate_image();
		repaint();
		// set_options();
	}

	public void adjustmentValueChanged(AdjustmentEvent e){
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
		int x=e.getX();
		int y=e.getY();
		if((x>10&&x<522)&&(y>10&&y<522)){
			inroi=true;
			startx=x;
			starty=y;
		}
	}

	public void mouseReleased(MouseEvent e){
		if(inroi){
			inroi=false;
			float xinc=((startx-10.0f)/512.0f)*(xmax-xmin);
			float yinc=((starty-10.0f)/512.0f)*(ymax-ymin);
			float xinc2=(((float)startx+(float)roiwidth-10.0f)/512.0f)*(xmax-xmin);
			float newwidth=xinc2-xinc;
			xmin+=xinc;
			ymin+=yinc;
			xmax=xmin+newwidth;
			ymax=ymin+newwidth;
			xminval.setText(""+xmin);
			xmaxval.setText(""+xmax);
			yminval.setText(""+ymin);
			ymaxval.setText(""+ymax);
			scaling=512.0f/newwidth;
			generate_fimage();
			generate_image();
			repaint();
		}
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		if(inroi){
			int x=e.getX();
			int y=e.getY();
			int xwidth=x-startx;
			int ywidth=y-starty;
			if(ywidth>xwidth){
				roiwidth=ywidth;
			}else{
				roiwidth=xwidth;
			}
			if(roiwidth<0){
				roiwidth=0;
			}
			repaint();
		}
	}

	private void generate_fimage(){
		PALMutils pu=new PALMutils(stdev*resmult,baseline,gain,512,512);
		activitylabel.setText("Calculating Image");
		float[] limits={xmin,xmax,ymin,ymax};
		float mult=100.0f;
		if(showorig){
			fdispimg=pu.render_const_size(molecules,limits,stdev*resmult,1.0f);
			fthreshimg=new float[512*512];
		}else{
			float tampmin=get_gaus_amp(ampmin*gain,stdev);
			float tampmax=get_gaus_amp(ampmax*gain,stdev);
			float[][] threshvals={null,null,{tampmin,tampmax},{c2min,c2max},null,null};
			float[][][] filtmol=pu.get_filtered_molecules(molecules,threshvals);
			if(showthresh){
				if(addnoise){
					fdispimg=pu.render_err_size(filtmol[0],limits,mult);
					fthreshimg=pu.render_err_size(filtmol[1],limits,mult);
				}else{
					fdispimg=pu.render_pix_size(filtmol[0],limits,mult);
					fthreshimg=pu.render_pix_size(filtmol[1],limits,mult);
				}
			}else{
				if(addnoise){
					fdispimg=pu.render_err_size(filtmol[0],limits,mult);
				}else{
					fdispimg=pu.render_pix_size(filtmol[0],limits,mult);
				}
				fthreshimg=new float[512*512];
			}
		}
		activitylabel.setText("Ready");
		// new ImagePlus("Test",new
		// FloatProcessor(512,512,fdispimg,null)).show();
	}

	private void generate_image(){
		int[] pixels=new int[512*512];
		for(int i=0;i<512*512;i++){
			int dumint=(int)(255.0f*((fdispimg[i]-dispmin)/(dispmax-dispmin)));
			int red=dumint;
			int bg=dumint;
			if(showthresh&&!showorig){
				red=(int)(255.0f*((fthreshimg[i]-dispmin)/(dispmax-dispmin)));
			}
			if(red>255){
				red=255;
			}
			if(red<0){
				red=0;
			}
			if(bg>255){
				bg=255;
			}
			if(bg<0){
				bg=0;
			}
			pixels[i]=0xff000000+(red<<16)+(bg<<8)+bg;
		}
		ColorProcessor cp=new ColorProcessor(512,512,pixels);
		dispimg=cp.createImage();
	}

	public float get_gaus_integral(float amp,float stdev){
		return amp*stdev*stdev*2.0f*(float)Math.PI;
	}

	public float get_gaus_amp(float integral,float stdev){
		return integral/(stdev*stdev*2.0f*(float)Math.PI);
	}

}
