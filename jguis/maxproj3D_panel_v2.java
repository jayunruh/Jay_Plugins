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
import ij.ImageListener;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
import j3D.element3D;
import j3D.j3Dutils;
import j3D.line3D;
import j3D.point3D;
import j3D.rect3D;
import j3D.renderer;
import j3D.spot3D;
import jalgs.matrixsolve;
import jalgs.profiler;
import jalgs.jseg.jsobel;
import jalgs.jsim.rngs;

import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Label;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JButton;
import javax.swing.JPanel;

public class maxproj3D_panel_v2 extends JPanel implements ActionListener,ItemListener,ImageListener{
	//in this version the data cube is static while the viewer rotates

	private JButton rotright,rotleft,rotup,rotdown,rotclock,rotcounter,reset,makemovie,droppoint,delpoint,exportbutton;
	private JButton profilebutton,importbutton;
	private Choice methodchoice,projchoice,roichoice;
	private Label threshlabel,methodlabel,xysteplabel,projlabel,limitslabel,roilabel;
	private TextField threshval,xystepval,xstartval,xendval,ystartval,yendval,zstartval,zendval;
	public ImagePlus imp3D,srcimp;
	public ImageStack stack,stack3D;
	public int method,width,height,slices,channels,frames,newzslices,maxsize,imagetype;
	public int currframe,currchan,xystep,nthreads,xstart1,xend1,ystart1,yend1,zstart1,zend1;
	public int xrem,yrem,zrem,zremcorr,projmeth,centerindex,topindex,selroi;
	public float zratio,xrot,yrot,zrot,invzr,halfzr;
	private float[] distances,thresh;
	private float[][] images;
	private rngs random;
	private element3D[] elements;
	private line3D[] lines;
	private renderer r,r2;
	private long starttime;
	public long timeout=5L;

	public static Frame launch_frame(maxproj3D_panel_v2 panel){
		final Frame f=new Frame("3D Viewer");
		f.setLocation(300,50);
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
		panel.setBounds(10,40,440,400);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(450,450));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(ImagePlus imp,float zratio,float thresh){
		init(imp,zratio,thresh,1);
	}

	public void init(ImagePlus imp,float zratio,float thresh,int nthreads){
		setLayout(null);
		this.zratio=zratio;
		this.nthreads=nthreads;
		xystep=1;
		srcimp=imp;
		int bsize=55;
		// 55
		int bheight=30;
		int theight=20;
		Font bfont=new Font("Lucida Sans Unicode",Font.BOLD,12);
		rotup=new JButton("\u2191");
		rotup.setFont(bfont);
		rotup.setBounds(10,10,bsize,bheight);
		rotup.addActionListener(this);
		add(rotup);
		rotdown=new JButton("\u2193");
		rotdown.setFont(bfont);
		rotdown.setBounds(10+bsize,10,bsize,bheight);
		rotdown.addActionListener(this);
		add(rotdown);
		rotleft=new JButton("\u2190");
		rotleft.setFont(bfont);
		rotleft.setBounds(10+2*bsize,10,bsize,bheight);
		rotleft.addActionListener(this);
		add(rotleft);
		rotright=new JButton("\u2192");
		rotright.setFont(bfont);
		rotright.setBounds(10+3*bsize,10,bsize,bheight);
		rotright.addActionListener(this);
		add(rotright);
		rotclock=new JButton("\u21b6");
		rotclock.setFont(bfont);
		rotclock.setBounds(10+4*bsize,10,bsize,bheight);
		rotclock.addActionListener(this);
		add(rotclock);
		rotcounter=new JButton("\u21b7");
		rotcounter.setFont(bfont);
		rotcounter.setBounds(10+5*bsize,10,bsize,bheight);
		rotcounter.addActionListener(this);
		add(rotcounter);
		reset=new JButton("Reset");
		reset.setFont(bfont);
		reset.setBounds(10+6*bsize,10,75,bheight);
		reset.addActionListener(this);
		add(reset);
		methodlabel=new Label("Interpolation");
		methodlabel.setBounds(10,bheight+20,80,20);
		add(methodlabel);
		methodchoice=new Choice();
		methodchoice.add("None");
		methodchoice.add("Random");
		methodchoice.add("Linear");
		methodchoice.setBounds(90,bheight+20,100,20);
		methodchoice.addItemListener(this);
		add(methodchoice);
		projlabel=new Label("Projection");
		projlabel.setBounds(200,bheight+20,70,20);
		add(projlabel);
		projchoice=new Choice();
		projchoice.add("Max");
		projchoice.add("First");
		projchoice.add("Avg");
		projchoice.add("Surface");
		projchoice.setBounds(280,bheight+20,100,20);
		projchoice.addItemListener(this);
		add(projchoice);
		threshlabel=new Label("Threshold");
		threshlabel.setBounds(10,bheight+50,80,20);
		add(threshlabel);
		threshval=new TextField(""+thresh);
		threshval.setBounds(90,bheight+50,60,20);
		threshval.addActionListener(this);
		add(threshval);
		xysteplabel=new Label("Resample");
		xysteplabel.setBounds(160,bheight+50,80,20);
		add(xysteplabel);
		xystepval=new TextField(""+xystep);
		xystepval.setBounds(250,bheight+50,60,20);
		xystepval.addActionListener(this);
		add(xystepval);
		
		width=imp.getWidth();
		height=imp.getHeight();
		stack=imp.getStack();
		slices=imp.getNSlices();
		channels=imp.getNChannels();
		frames=imp.getNFrames();
		currframe=imp.getFrame()-1;
		currchan=imp.getChannel()-1;
		xstart1=0;
		xend1=width;
		ystart1=0;
		yend1=height;
		zstart1=0;
		zend1=slices;
		
		limitslabel=new Label("Limits:");
		limitslabel.setBounds(10,bheight+80,80,20);
		add(limitslabel);
		xstartval=new TextField(""+xstart1);
		xstartval.setBounds(10,bheight+110,60,20);
		xstartval.addActionListener(this);
		add(xstartval);
		xendval=new TextField(""+xend1);
		xendval.setBounds(75,bheight+110,60,20);
		xendval.addActionListener(this);
		add(xendval);
		ystartval=new TextField(""+ystart1);
		ystartval.setBounds(140,bheight+110,60,20);
		ystartval.addActionListener(this);
		add(ystartval);
		yendval=new TextField(""+yend1);
		yendval.setBounds(205,bheight+110,60,20);
		yendval.addActionListener(this);
		add(yendval);
		zstartval=new TextField(""+zstart1);
		zstartval.setBounds(270,bheight+110,60,20);
		zstartval.addActionListener(this);
		add(zstartval);
		zendval=new TextField(""+zend1);
		zendval.setBounds(335,bheight+110,60,20);
		zendval.addActionListener(this);
		add(zendval);
		
		makemovie=new JButton("Make Movie");
		makemovie.setFont(bfont);
		makemovie.setBounds(10,bheight+140,140,bheight);
		makemovie.addActionListener(this);
		add(makemovie);
		droppoint=new JButton("Drop Point");
		droppoint.setFont(bfont);
		droppoint.setBounds(10+140,bheight+140,120,bheight);
		droppoint.addActionListener(this);
		add(droppoint);
		delpoint=new JButton("Del Point");
		delpoint.setFont(bfont);
		delpoint.setBounds(10+140+120,bheight+140,120,bheight);
		delpoint.addActionListener(this);
		add(delpoint);
		
		roilabel=new Label("Sel Points");
		roilabel.setBounds(10,bheight+140+bheight+15,80,20);
		add(roilabel);
		roichoice=new Choice();
		roichoice.add("None");
		roichoice.add("All");
		roichoice.setBounds(100,bheight+140+bheight+15,70,20);
		roichoice.addItemListener(this);
		add(roichoice);
		selroi=0;
		exportbutton=new JButton("Export Pts");
		exportbutton.setFont(bfont);
		exportbutton.setBounds(180,bheight+140+bheight+10,120,30);
		exportbutton.addActionListener(this);
		add(exportbutton);
		profilebutton=new JButton("Straighten Profile");
		profilebutton.setFont(bfont);
		profilebutton.setBounds(10,bheight+140+bheight+10+bheight+10,170,30);
		profilebutton.addActionListener(this);
		add(profilebutton);
		importbutton=new JButton("Import Pts");
		importbutton.setFont(bfont);
		importbutton.setBounds(180,bheight+140+bheight+10+bheight+10,120,30);
		importbutton.addActionListener(this);
		add(importbutton);
		
		this.thresh=new float[channels];
		for(int i=0;i<channels;i++)
			this.thresh[i]=thresh;
		newzslices=(int)(zratio*slices);
		maxsize=(int)Math.sqrt(width*width+height*height+newzslices*newzslices);
		xrem=(maxsize-width)/2;
		yrem=(maxsize-height)/2;
		zrem=(maxsize-newzslices)/2;
		zremcorr=(int)(zrem/zratio);
		lines=new line3D[maxsize*maxsize];
		elements=new element3D[maxsize*maxsize];
		for(int i=0;i<maxsize;i++){
			for(int j=0;j<maxsize;j++){
				lines[j+i*maxsize]=new line3D(j,i,0,j,i,maxsize,Color.black);
				elements[j+i*maxsize]=lines[j+i*maxsize];
			}
		}
		r=new renderer(maxsize,maxsize);
		r.centerz=maxsize/2;
		r2=new renderer(maxsize,maxsize);
		r2.centerz=r.centerz;
		element3D[] cube=new element3D[1];
		cube[0]=new rect3D(r.centerx,r.centery,r.centerz,width,height,newzslices,Color.white);
		//cube[0]=new cube3D(xrem,yrem,zrem,5,Color.white);
		r2.setelementarray(cube);
		r2.background=Color.black;
		r.setelementarray(elements);
		stack3D=new ImageStack(maxsize,maxsize);
		images=new float[channels+1][maxsize*maxsize];
		for(int i=0;i<channels;i++)
			stack3D.addSlice("",images[i]);
		imp3D=new ImagePlus("3D Rendering",stack3D);
		imp3D.setOpenAsHyperStack(true);
		imp3D.setDimensions(channels,1,1);
		imp3D=new CompositeImage(imp3D,CompositeImage.COMPOSITE);
		imp3D.copyScale(imp);
		//if(imp.isComposite()){
			/*((CompositeImage)imp3D).copyLuts(imp);
			LUT[] luts=((CompositeImage)imp).getLuts();
			LUT[] luts2=new LUT[channels];
			System.arraycopy(luts,0,luts2,0,channels);
			//luts2[channels]=LUT.createLutFromColor(Color.white);
			((CompositeImage)imp3D).setLuts(luts2);*/
		//}else{
			//((CompositeImage)imp3D).setChannelLut(LUT.createLutFromColor(Color.white),channels+1);
		//}
		imp3D.show();
		//r.set_perspective(maxsize*2.0);
		xrot=0.0f;
		yrot=0.0f;
		zrot=0.0f;
		init_distances();
		method=0;
		projmeth=0;
		random=new rngs();
		imagetype=0;
		if(stack.getPixels(1) instanceof short[])
			imagetype=1;
		if(stack.getPixels(1) instanceof float[])
			imagetype=2;
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		r2.setrotation((int)xrot,(int)yrot,(int)zrot);
		topindex=(int)(0.5f*maxsize);
		centerindex=topindex*maxsize+topindex;
		ImagePlus.addImageListener(this);
		update_image();
		((CompositeImage)imp3D).resetDisplayRanges();
	}

	public void setVisible(boolean b){
		ImagePlus.removeImageListener(this);
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==makemovie){
			make_movie();
			return;
		}
		if(e.getSource()==exportbutton){
			GenericDialog gd=new GenericDialog("Options");
			gd.addCheckbox("To RoiManager",false);
			gd.showDialog(); if(gd.wasCanceled()) return;
			boolean torman=gd.getNextBoolean();
			float[][] coords=new float[3][r2.elements.length-1];
			for(int i=0;i<r2.elements.length-1;i++){
				coords[0][i]=((spot3D)r2.elements[i+1]).point.x-xrem;
				coords[1][i]=((spot3D)r2.elements[i+1]).point.y-yrem;
				coords[2][i]=((spot3D)r2.elements[i+1]).point.z-zrem;
				if(torman) coords[2][i]/=zratio;
			}
			if(torman){
				RoiManager rman=RoiManager.getInstance();
				if(rman==null) rman=new RoiManager();
				for(int i=0;i<coords[0].length;i++){
					PointRoi roi=new PointRoi(coords[0][i],coords[1][i]);
					int slice=(int)coords[2][i];
					int position=currframe*slices*channels+slice*channels+currchan+1;
					roi.setPosition(position);
					rman.addRoi(roi);
				}
			} else {
				new PlotWindow3D("Plot 3D",new Traj3D("x","y","z",coords[0],coords[1],coords[2])).draw();
			}
			return;
		}
		if(e.getSource()==profilebutton){
			GenericDialog gd=new GenericDialog("Options");
			gd.addNumericField("Line Thickness",5,2);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			int linethickness=(int)gd.getNextNumber();
			float[][] coords=new float[3][r2.elements.length-1];
			for(int i=0;i<r2.elements.length-1;i++){
				coords[0][i]=((spot3D)r2.elements[i+1]).point.x-xrem;
				coords[1][i]=((spot3D)r2.elements[i+1]).point.y-yrem;
				coords[2][i]=((spot3D)r2.elements[i+1]).point.z-zrem;
			}
			float[][] profile=new float[channels][];
			float[][][] straightened=new float[channels][][];
			for(int j=0;j<channels;j++){
				Object[] stack2=jutils.get3DZSeries(stack,j,currframe,frames,slices,channels);
				profile[j]=profiler.get3DThickProfile(stack2,width,height,coords[0],coords[1],coords[2],false,linethickness,0,zratio);
				straightened[j]=profiler.get3DThickStraightened(stack2,width,height,coords[0],coords[1],coords[2],false,linethickness,0,zratio);
			}
			int length=profile[0].length;
			new PlotWindow4("3D Profile","length","intensity",profile,null).draw();
			int straightlength=length;
			ImageStack sstack=new ImageStack(linethickness,straightlength);
			for(int i=0;i<straightened[0].length;i++){
				for(int k=0;k<channels;k++){
					sstack.addSlice("",straightened[k][i]);
				}
			}
			ImagePlus simp=new ImagePlus("Straightened Profile",sstack);
			simp.setOpenAsHyperStack(true);
			simp.setDimensions(channels,linethickness,1);
			simp.copyScale(srcimp);
			simp.show();
			return;
		}
		if(e.getSource()==droppoint){
			Roi roi=imp3D.getRoi();
			Rectangle r=roi.getBounds();
			int index=r.x+r.y*maxsize;
			float[] temp=get_line_max(lines[index],distances[index],currchan);
			spot3D spot=new spot3D((int)temp[2],(int)temp[3],(int)temp[4],0,Color.white);
			int opt=2;
			if(selroi>1){
				GenericDialog gd=new GenericDialog("Options");
				String[] options={"Insert After","Replace","Append"};
				gd.addChoice("Options",options,options[0]);
				gd.showDialog(); if(gd.wasCanceled()) return;
				opt=gd.getNextChoiceIndex();
			}
			if(opt==2){
				r2.addelement(spot);
				roichoice.add(""+(r2.elements.length-1));
			} else if(opt==1){
				r2.replaceelement(spot,selroi-1);
				r2.elements[selroi-1].color=Color.yellow;
			} else {
				r2.insertelement(spot,selroi-1);
				roichoice.insert(""+(selroi),selroi+1);
				String[] list=get_dropdown_list(roichoice);
				for(int i=selroi+2;i<=r2.elements.length;i++){
					list[i]=""+(i-1);
				}
				set_dropdown_list(roichoice,list);
			}
		}
		if(e.getSource()==importbutton){
			ImageWindow[] iw=jutils.selectPlotFamily(true,1);
			if(iw==null) return;
			if(iw[0]!=null){
				GenericDialog gd=new GenericDialog("Options");
				gd.addCheckbox("Scale Profile",true);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				boolean scaleprof=gd.getNextBoolean();
				float psize=(float)jutils.get_psize(srcimp);
				float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getXValues");
				float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
				float[][][] zvals=(float[][][])jutils.runPW4VoidMethod(iw[0],"getZValues");
    			for(int i=0;i<xvals[0].length;i++){
    				float x=xvals[0][i];
    				float y=yvals[0][i];
    				float z=zvals[0][0][i];
    				if(scaleprof){x/=psize; y/=psize; z/=psize;}
    				spot3D spot=new spot3D((int)x+xrem,(int)y+yrem,(int)z+zrem,0,Color.white);
    				r2.addelement(spot);
    				roichoice.add(""+(r2.elements.length-1));
    			}
			} else {
				RoiManager rman=RoiManager.getInstance();
				if(rman==null) return;
				Roi[] rois=rman.getRoisAsArray();
    			for(int i=0;i<rois.length;i++){
    				int x=rois[i].getBounds().x;
    				int y=rois[i].getBounds().y;
    				int z=rois[i].getPosition();
    				z=(int)(z*zratio);
    				spot3D spot=new spot3D(x+xrem,y+yrem,z+zrem,0,Color.white);
    				r2.addelement(spot);
    				roichoice.add(""+(r2.elements.length-1));
    			}
			}
		}
		if(e.getSource()==delpoint){
			/*Roi roi=imp3D.getRoi();
			if(roi==null) return;
			Rectangle r=roi.getBounds();
			int index=r.x+r.y*maxsize;
			float[] temp=get_line_max(lines[index],distances[index],currchan);*/
			//find the closest point
			element3D[] spots=r2.elements;
			if(spots.length<2) return;
			/*if(spots.length<3){
				r2.removeelement(1);
				roichoice.remove(1);
				roichoice.select(0);
				IJ.showStatus("deleting"+1);
			} else {
    			float mindist=get_2D_dist(((spot3D)spots[1]).point,new point3D(r.x,r.y,(int)temp[1]));
    			int minpt=0;
    			for(int i=1;i<(spots.length-1);i++){
    				float temp2=get_2D_dist(((spot3D)spots[i+1]).point,new point3D(r.x,r.y,(int)temp[1]));
    				if(temp2<mindist){
    					mindist=temp2;
    					minpt=i;
    				}
    			}
				
    			IJ.showStatus("deleting"+(minpt+1));
    			r2.removeelement(minpt+1);
			}*/
			int selroi=roichoice.getSelectedIndex();
			if(selroi==0) return;
			if(selroi==1){
				int nrois=r2.elements.length-1;
				for(int i=0;i<nrois;i++){
					r2.removeelement(1);
					roichoice.remove(2);
				}
				roichoice.select(0);
			} else {
    			IJ.showStatus("deleting"+selroi);
    			r2.removeelement(selroi-1);
    			roichoice.remove(selroi);
    			roichoice.select(0);
			}
		}
		if(e.getSource()==rotup){
			xrot+=10.0f;
		} if(e.getSource()==rotdown){
			xrot-=10.0f;
		} if(e.getSource()==rotright){
			yrot+=10.0f;
		}if(e.getSource()==rotleft){
			yrot-=10.0f;
		}if(e.getSource()==rotclock){
			zrot-=10.0f;
		}if(e.getSource()==rotcounter){
			zrot+=10.0f;
		}if(e.getSource()==reset){
			xrot=0.0f;
			yrot=0.0f;
			zrot=0.0f;
		}
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		update_cube();
		float temp=Float.parseFloat(threshval.getText());
		thresh[currchan]=temp;
		xystep=(int)Float.parseFloat(xystepval.getText());
		xstart1=(int)Float.parseFloat(xstartval.getText());
		xend1=(int)Float.parseFloat(xendval.getText());
		ystart1=(int)Float.parseFloat(ystartval.getText());
		yend1=(int)Float.parseFloat(yendval.getText());
		zstart1=(int)Float.parseFloat(zstartval.getText());
		zend1=(int)Float.parseFloat(zendval.getText());
		update_image();
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==methodchoice){
			method=methodchoice.getSelectedIndex();
			update_image();
		}
		if(e.getSource()==projchoice){
			projmeth=projchoice.getSelectedIndex();
			update_image();
		}
		if(e.getSource()==roichoice){
			for(int i=0;i<r2.elements.length;i++) r2.elements[i].color=Color.white;
			selroi=roichoice.getSelectedIndex();
			if(selroi==1){
				for(int i=1;i<r2.elements.length;i++) r2.elements[i].color=Color.yellow;
			}
			if(selroi>1){
				r2.elements[selroi-1].color=Color.yellow;
			}
			int[] overpix=r2.renderrgb();
			ColorProcessor cp=new ColorProcessor(maxsize,maxsize,overpix);
			ImageRoi roi=new ImageRoi(0,0,cp);
			roi.setZeroTransparent(true);
			imp3D.setOverlay(new Overlay(roi));
			//update_image();
		}
	}

	public void imageClosed(ImagePlus imp){
	}

	public void imageOpened(ImagePlus imp){
	}

	public void imageUpdated(ImagePlus imp){
		if(imp.equals(imp3D)){
			int tempcurrchan=imp3D.getChannel()-1;
			if(tempcurrchan!=currchan&&tempcurrchan<channels){
				currchan=tempcurrchan;
				srcimp.setPosition(currchan+1,srcimp.getSlice(),currframe+1);
				threshval.setText(""+thresh[currchan]);
			}
		}else{
			if(imp.equals(srcimp)){
				int tempcurrchan=srcimp.getChannel()-1;
				if(tempcurrchan!=currchan){
					currchan=tempcurrchan;
					imp3D.setPosition(currchan+1,1,1);
					threshval.setText(""+thresh[currchan]);
				}
				int tempcurrframe=srcimp.getFrame()-1;
				if(tempcurrframe!=currframe){
					currframe=tempcurrframe;
					update_image();
				}
			}
		}
	}

	public void init_distances(){
		distances=new float[maxsize*maxsize];
		for(int i=0;i<maxsize*maxsize;i++){
			float distx=(lines[i].pt2.rx-lines[i].pt1.rx);
			float disty=(lines[i].pt2.ry-lines[i].pt1.ry);
			float distz=(lines[i].pt2.rz-lines[i].pt1.rz);
			distances[i]=(float)Math.sqrt(distx*distx+disty*disty+distz*distz);
		}
	}

	public void update_image(){
		starttime=System.currentTimeMillis();
		float[][] images2=maxprojimage();
		long temp=System.currentTimeMillis();
		//IJ.log("totprojtime = "+(temp-starttime)/1000.0f);
		for(int i=0;i<channels;i++){
			System.arraycopy(images2[i],0,images[i],0,maxsize*maxsize);
		}
		long temp2=System.currentTimeMillis();
		//IJ.log("totcopytime = "+(temp2-temp)/1000.0f);
		int[] overpix=r2.renderrgb();
		ColorProcessor cp=new ColorProcessor(maxsize,maxsize,overpix);
		ImageRoi roi=new ImageRoi(0,0,cp);
		roi.setZeroTransparent(true);
		imp3D.setOverlay(new Overlay(roi));
		/*float[] grid=r2.renderfloat(255.0f);
		System.arraycopy(grid,0,images[channels],0,maxsize*maxsize);
		long temp3=System.currentTimeMillis();*/
		//IJ.log("totgridtime = "+(temp3-temp2)/1000.0f);
		imp3D.updateAndDraw();
	}
	
	public float get_3D_nearest_value(float x,float y,float z,int channel){
		int x2=Math.round(x)-xrem;
		int y2=Math.round(y)-yrem;
		float scaledz=z-zremcorr;
		int z2=Math.round(scaledz);
		if(x2<0) return 0.0f;
		if(y2<0) return 0.0f;
		if(z2<0) return 0.0f;
		if(x2>=width) return 0.0f;
		if(y2>=height) return 0.0f;
		if(z2>=slices) return 0.0f;
		if(method==0){
			return get_value(x2+y2*width,z2,channel);
		} else if(method==1) {
			return randinterpz(x2+y2*width,scaledz,channel);
		} else {
			return interpolatez(x2+y2*width,scaledz,channel);
		}
	}

	public float interpolatez(int index,float scaledz,int channel){
		int zprev=(int)scaledz;
		if(scaledz==zprev) return get_value(index,zprev,channel);
		if(zprev<0) return 0.0f;
		int znext=zprev+1;
		if(znext>=slices) return 0.0f;
		float zrem2=scaledz-zprev;
		float z1=get_value(index,zprev,channel);
		float z2=get_value(index,znext,channel);
		return z1+zrem2*(z2-z1);
	}
	
	public float randinterpz(int index,float scaledz,int channel){
		int zprev=(int)scaledz;
		if(zprev<0) return 0.0f;
		if(scaledz==zprev) return get_value(index,zprev,channel);
		if(zprev>=(slices-1)) return 0.0f;
		float zrem2=scaledz-zprev;
		if((float)random.unidev(1.0,0.0)<zrem2) return get_value(index,zprev,channel);
		else return get_value(index,zprev+1,channel);
	}
	
	private String[] get_dropdown_list(Choice choice){
		int nchoices=choice.getItemCount();
		String[] choices=new String[nchoices];
		for(int i=0;i<choices.length;i++){
			choices[i]=choice.getItem(i);
		}
		return choices;
	}
	
	private void set_dropdown_list(Choice choice,String[] list){
		int sel=choice.getSelectedIndex();
		choice.removeAll();
		for(int i=0;i<list.length;i++) choice.add(list[i]);
		if(sel>=list.length) sel=list.length-1;
		choice.select(sel);
	}
	
	private float get_3D_dist(point3D pt1,point3D pt2){
		return (float)Math.sqrt((pt1.rx-pt2.rx)*(pt1.rx-pt2.rx)+(pt1.ry-pt2.ry)*(pt1.ry-pt2.ry)+(pt1.rz-pt2.rz)*(pt1.rz-pt2.rz));
	}
	
	private float get_2D_dist(point3D pt1,point3D pt2){
		return (float)Math.sqrt((pt1.rx-pt2.rx)*(pt1.rx-pt2.rx)+(pt1.ry-pt2.ry)*(pt1.ry-pt2.ry));
	}
	
	private void update_cube(){
		//we need to transform from our observational coordinate system back to the stack coordinate system
		//firstly get the observational coordinate system
		/*line3D center=lines[centerindex];
		line3D top=lines[topindex];
		line3D right=lines[centerindex+topindex-1];
		double[] obsz=center.unitvec();
		obsz[0]=-obsz[0]; obsz[1]=-obsz[1]; obsz[2]=-obsz[2];
		double[] obsy={top.pt1.rx-center.pt1.rx,top.pt1.ry-center.pt1.ry,top.pt1.rz-center.pt1.rz};
		obsy=j3Dutils.norm_vec(obsy);
		double[] obsx={right.pt1.rx-center.pt1.rx,right.pt1.ry-center.pt1.ry,right.pt1.rz-center.pt1.rz};
		obsx=j3Dutils.norm_vec(obsx);
		double[][] transmat={obsx,obsy,obsz};
		IJ.log(table_tools.print_double_array(transmat));
		//first transform into the cartesian system
		double[][] inv=(new matrixsolve()).gjinv(transmat,3);
		//double[][] inv=transmat;
		//now transform into our data coordinate system (reflect y and z)
		inv[1][0]=-inv[1][0]; inv[1][1]=-inv[1][1]; inv[1][2]=-inv[1][2];
		inv[2][0]=-inv[2][0]; inv[2][1]=-inv[2][1]; inv[2][2]=-inv[2][2];
		//inv[0]=j3Dutils.norm_vec(inv[0]);
		//inv[1]=j3Dutils.norm_vec(inv[1]);
		//inv[2]=j3Dutils.norm_vec(inv[2]);
		//IJ.log(table_tools.print_double_array(inv));
		r2.transform(inv);*/
		double[] origin=find_coordinate(width/2,height/2,newzslices/2);
		double[] right=find_coordinate(width,height/2,newzslices/2);
		double[] bot=find_coordinate(width/2,height,newzslices/2);
		double[] back=find_coordinate(width/2,height/2,newzslices);
		//IJ.log("Origin: "+origin[0]+" , "+origin[1]+" , "+origin[2]);
		//float[] xvals={(float)origin[0],(float)right[0],(float)bot[0],(float)back[0]};
		//float[] yvals={(float)origin[1],(float)right[1],(float)bot[1],(float)back[1]};
		//imp3D.setRoi(new PointRoi(xvals,yvals,4));
		double[] xvec={right[0]-origin[0],right[1]-origin[1],right[2]-origin[2]};
		xvec=j3Dutils.norm_vec(xvec);
		double[] yvec={bot[0]-origin[0],bot[1]-origin[1],bot[2]-origin[2]};
		yvec=j3Dutils.norm_vec(yvec);
		double[] zvec={back[0]-origin[0],back[1]-origin[1],back[2]-origin[2]};
		zvec=j3Dutils.norm_vec(zvec);
		double[][] transmat={xvec,yvec,zvec};
		transmat=matrixsolve.transpose(transmat);
		r2.transform(transmat);
	}
	
	private double[] find_coordinate(int x,int y,int z){
		//this finds the projection coordinates for a point in data space
		//find the distance from each top corner
		point3D pt=new point3D(x+xrem,y+yrem,z+zrem);
		float[] temp1=get_closest_dist(lines[0],pt);
		float[] temp2=get_closest_dist(lines[maxsize-1],pt);
		float d1=temp1[0];
		float d2=temp2[0];
		float d3=maxsize;
		//now triangulate
		double costheta1=(d1*d1-d2*d2+d3*d3)/(2.0*d1*d3);
		double costheta2=(d2*d2-d1*d1+d3*d3)/(2.0*d2*d3);
		double sintheta1=Math.sqrt(1.0-costheta1*costheta1);
		double x2=d3-d2*costheta2;
		double y2=d1*sintheta1;
		//IJ.log("distances: "+d1+" , "+d2);
		return new double[]{x2,y2,0.5*(temp1[1]+temp2[1])};
	}

	public float[] get_closest_dist(line3D line,point3D pt3){
		// get the closest distance between p3 and the line p1-p2
		// first get the distance of the tangent intersection from x1,y1
		// this is the dot product of the line and the line start to the point
		// divided by the length of the line
		float length=get_3D_dist(line.pt1,line.pt2);
		float interdist=((line.pt2.rx-line.pt1.rx)*(pt3.rx-line.pt1.rx)+(line.pt2.ry-line.pt1.ry)*(pt3.ry-line.pt1.ry)+(line.pt2.rz-line.pt1.rz)*(pt3.rz-line.pt1.rz))/length;
		// don't allow the intersection to move outside the line segment
		if(interdist>length)
			interdist=length;
		if(interdist<0.0f)
			interdist=0.0f;
		// now get the tangent intersection
		float xinc=(line.pt2.rx-line.pt1.rx)/length;
		float yinc=(line.pt2.ry-line.pt1.ry)/length;
		float zinc=(line.pt2.rz-line.pt1.rz)/length;
		float interx=xinc*interdist+line.pt1.rx;
		float intery=yinc*interdist+line.pt1.ry;
		float interz=zinc*interdist+line.pt1.rz;
		return new float[]{get_3D_dist(new point3D((int)interx,(int)intery,(int)interz),pt3),interdist};
	}

	public float get_value(int index,int z,int channel){
		int index2=currframe*slices*channels+z*channels+channel;
		return get_value(index,index2);
	}

	public float get_value(int index,int index2){
		if(imagetype==0){
			return ((byte[])stack.getPixels(index2+1))[index]&0xff;
		}else{
			if(imagetype==1){
				return ((short[])stack.getPixels(index2+1))[index]&0xffff;
			}else{
				return ((float[])stack.getPixels(index2+1))[index];
			}
		}
	}

	public float[] get_line_max(float xstart,float ystart,float zstart,float xinc,float yinc,float zinc,int length,int channel){
		float x=xstart-xinc;
		float y=ystart-yinc;
		float scaledzinc=zinc/zratio;
		float z=zstart/zratio-scaledzinc;
		// here xinc, yinc, and zinc are the dimensional steps corresponding to a pixel unit
		float max=0.0f; float maxpos=Float.NaN;
		for(int k=0;k<length;k++){
			x+=xinc;
			y+=yinc;
			z+=scaledzinc;
			if(x<xrem+xstart1) continue;
			if(y<yrem+ystart1) continue;
			if(z<zremcorr+zstart1) continue;
			if(x>=xrem+xend1) continue;
			if(y>=yrem+yend1) continue;
			if(z>=zremcorr+zend1) continue;
			float temp=get_3D_nearest_value(x,y,z,channel);
			if(temp<thresh[channel]) temp=0.0f;
			if(projmeth==0 || projmeth==3){
				if(temp>max){max=temp; maxpos=(float)k;}
			} else if(projmeth==1){
				max=temp;
				if(temp>0.0f) return new float[]{max,(float)k};
			} else {
				max+=temp;
			}
		}
		if(projmeth==2) max/=length;
		return new float[]{max,maxpos};
	}
	
	public float[] get_line_max(line3D line,float dist,int channel){
		float distx=(line.pt2.rx-line.pt1.rx);
		float disty=(line.pt2.ry-line.pt1.ry);
		float distz=(line.pt2.rz-line.pt1.rz);
		float xinc=distx/dist;
		float yinc=disty/dist;
		float zinc=distz/dist;
		float xstart=line.pt1.rx;
		float ystart=line.pt1.ry;
		float zstart=line.pt1.rz;
		float x=xstart-xinc;
		float y=ystart-yinc;
		float scaledzinc=zinc/zratio;
		float z=zstart/zratio-scaledzinc;
		// here xinc, yinc, and zinc are the dimensional steps corresponding to a pixel unit
		float max=0.0f;
		float maxk=0.0f;
		float maxx=0.0f,maxy=0.0f,maxz=0.0f;
		for(int k=0;k<maxsize;k++){
			x+=xinc;
			y+=yinc;
			z+=scaledzinc;
			if(x<xrem+xstart1) continue;
			if(y<yrem+ystart1) continue;
			if(z<zremcorr+zstart1) continue;
			if(x>=xrem+xend1) continue;
			if(y>=yrem+yend1) continue;
			if(z>=zremcorr+zend1) continue;
			float temp=get_3D_nearest_value(x,y,z,channel);
			if(temp<thresh[channel]) temp=0.0f;
			if(projmeth==0){
				if(temp>max){max=temp; maxk=k; maxx=x; maxy=y; maxz=z*zratio;}
			} else if(projmeth==1){
				max=temp;
				if(temp>0.0f) return new float[]{max,k, maxx=x, maxy=y, maxz=z*zratio};
			} else {
				max+=temp;
				maxk+=temp*k;
			}
		}
		if(projmeth==2){
			maxk/=max;
			max/=maxsize;
			maxx=xstart+xinc*maxk;
			maxy=ystart+yinc*maxk;
			maxz=zstart+zinc*maxk;
		}
		return new float[]{max,maxk,maxx,maxy,maxz};
	}

	public void make_movie(){
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_steps",36,0);
		gd.addNumericField("X_Angle_Step_(deg)",10.0,5,15,null);
		gd.addNumericField("Y_Angle_Step_(deg)",10.0,5,15,null);
		gd.addNumericField("Z_Angle_Step_(deg)",0.0,5,15,null);
		gd.addNumericField("Frame_Step",0,0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		int nsteps=(int)gd.getNextNumber();
		float dy=(float)gd.getNextNumber();
		float dx=(float)gd.getNextNumber();
		float dz=(float)gd.getNextNumber();
		int dt=(int)gd.getNextNumber();
		ImageStack movstack=new ImageStack(maxsize,maxsize);
		Overlay overlay=new Overlay();
		float oldxrot=xrot;
		float oldyrot=yrot;
		float oldzrot=zrot;
		int oldt=currframe;
		for(int i=0;i<nsteps;i++){
			float[][] temp=maxprojimage();
			for(int j=0;j<temp.length;j++){
				movstack.addSlice("",temp[j]);
			}
			int[] overpix=r2.renderrgb();
			ColorProcessor cp=new ColorProcessor(maxsize,maxsize,overpix);
			ImageRoi roi=new ImageRoi(0,0,cp);
			roi.setZeroTransparent(true);
			roi.setPosition(i+1);
			overlay.add(roi);
			//movstack.addSlice("",r2.renderfloat(255.0f));
			xrot+=dx;
			yrot+=dy;
			zrot-=dz;
			currframe+=dt;
			if(currframe>=frames)
				currframe=0;
			r.setrotation(xrot,yrot,zrot);
			update_cube();
			IJ.showProgress(i,nsteps);
			IJ.showStatus("frame "+(i+1)+" of "+nsteps+" rendered");
		}
		currframe=oldt;
		xrot=oldxrot;
		yrot=oldyrot;
		zrot=oldzrot;
		r.setrotation(xrot,yrot,zrot);
		update_cube();
		ImagePlus imp=new ImagePlus("Max Proj Movie",movstack);
		imp.setOverlay(overlay);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(channels,1,nsteps);
		imp=new CompositeImage(imp,CompositeImage.COMPOSITE);
		imp.copyScale(imp3D);
		((CompositeImage)imp).copyLuts(imp3D);
		imp.show();
	}
	
	private float[] sobel_point(float[] pixels,int x,int y){
		int offset=x+y*maxsize;
		float[] output=new float[2];
		if(x<=0)
			return output;
		if(x>=(maxsize-1))
			return output;
		if(y<=0)
			return output;
		if(y>=(maxsize-1))
			return output;
		if(Float.isNaN(pixels[offset])){
			return output;
		}
		float val=pixels[offset];
		float[][] n=new float[3][3];
		n[0][0]=pixels[offset-maxsize-1];
		n[0][1]=pixels[offset-maxsize];
		n[0][2]=pixels[offset-maxsize+1];
		n[1][0]=pixels[offset-1];
		n[1][1]=pixels[offset];
		n[1][2]=pixels[offset+1];
		n[2][0]=pixels[offset+maxsize-1];
		n[2][1]=pixels[offset+maxsize];
		n[2][2]=pixels[offset+maxsize+1];
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(Float.isNaN(n[i][j]))
					n[i][j]=val;
			}
		}
		output[0]=n[2][2]+2.0f*n[1][2]+n[0][2]-n[2][0]-2.0f*n[1][0]-n[0][0];
		output[1]=n[2][0]+2.0f*n[2][1]+n[2][2]-n[0][0]-2.0f*n[0][1]-n[0][2];
		return output;
	}

	public float[][] maxprojimage(){
		float[][] outimage=new float[channels][maxsize*maxsize];
		ExecutorService executor=null;
		//these are 0maxsize,1channels,2slices,3frames,4xystep,5xrem,6yrem,7zremcorr,8projmeth,9currframe,10imagetype,11method
		int[] intparams={maxsize,channels,slices,frames,xystep,xrem,yrem,zremcorr,projmeth,currframe,imagetype,method,xstart1,xend1,ystart1,yend1,zstart1,zend1};
		//Object[] stackarray=stack.getImageArray();
		Object[] stackarray=jutils.stack2array(stack);
		//long tstarttime=starttime;
		//float[][] depth=new float[channels][maxsize*maxsize];
		for(int ch=0;ch<channels;ch++){
			float[] depth=new float[maxsize*maxsize];
			//run a new thread pool for each channel
			if(nthreads>1) executor=Executors.newFixedThreadPool(nthreads);
			for(int i=0;i<maxsize;i+=xystep){
				if(nthreads<2){
					for(int l=0;l<maxsize;l+=xystep){	
						int index=i*maxsize+l;
						float distx=(lines[index].pt2.rx-lines[index].pt1.rx);
						float disty=(lines[index].pt2.ry-lines[index].pt1.ry);
						float distz=(lines[index].pt2.rz-lines[index].pt1.rz);
						float dist=distances[index];
						float xunit2=distx/dist;
						float yunit2=disty/dist;
						float zunit2=distz/dist;
						// float xunit2=xunit*zratio;
						// float yunit2=yunit*zratio;
						float xstart=lines[index].pt1.rx;
						float ystart=lines[index].pt1.ry;
						float zstart=lines[index].pt1.rz;
						float[] putval=get_line_max(xstart,ystart,zstart,xunit2,yunit2,zunit2,maxsize,ch);
    					//float putval=jstatistics.getstatistic("Max",line,null);
    					for(int j=i;j<(i+xystep);j++){
    						for(int k=l;k<(l+xystep);k++){
    							outimage[ch][k+j*maxsize]=putval[0];
    							if(projmeth==3) depth[k+j*maxsize]=putval[1];
    						}
    					}
					}
				} else {
					if(IJ.escapePressed()) break;
					float[] fltparams={zratio,thresh[ch]};
					//int[] intparams,float[] fltparams,ImageStack stack,float[] outimage,int width,int height,int xc,int yc,int ch
					Runnable worker=new interp_line_v2(intparams,fltparams,lines,distances,stackarray,outimage[ch],width,height,i,ch,depth);
					executor.execute(worker);
				}
			}
			if(nthreads>1){
				if(IJ.escapePressed()) executor.shutdownNow();
				else executor.shutdown();
				// Wait until all threads are finished
				try{
					boolean success=executor.awaitTermination(timeout,TimeUnit.SECONDS);
					if(!success) IJ.showStatus("3D Rendering Timeout");
				}catch(InterruptedException e){}
			}
			//if we have a surface, apply the normals--use sobel for the derivative
			if(projmeth==3){
				for(int i=1;i<(maxsize-1);i++){
					for(int l=1;l<(maxsize-1);l++){
						float[] der=sobel_point(depth,l,i);
						float mag=(float)Math.sqrt(der[0]*der[0]+der[1]*der[1]);
						//the normal is the cross product of the derivative vectors
						float[] norm={-der[0]/mag,-der[1]/mag,1.0f/mag};
						//now the intensity is calculated as the dot product of the norm and the light vector
						//the light will come in at a 30 deg angle to the y axis
						float dotprod=(float)Math.abs(norm[0]*0.5f+norm[2]*0.866f);
						outimage[ch][l+i*maxsize]*=dotprod;
    				}
    			}
			}
			if(IJ.escapePressed()) return outimage;
			//long temp=System.currentTimeMillis();
			//IJ.log("line "+i+" time = "+(temp-tstarttime)/1000.0f);
			//tstarttime=temp;
		}

		return outimage;
	}

}

class interp_line_v2 implements Runnable{
	private float[] fltparams; //these are 0xstart,1ystart,2zstart,3xinc,4yinc,5zinc,6zratio,7thresh
	private final int[] intparams; //these are 0maxsize,1channels,2slices,3frames,4xystep,5xrem,6yrem,7zremcorr,8method,9currframe,10imagetype,11projmeth
	private final int width,height;
	private int xc,yc,ch;
	//public ImageStack stack;
	public Object[] stack;
	private float[] outimage;
	private float[] depth;
	public line3D[] lines;
	public float[] distances;
	rngs random;
	
	public interp_line_v2(int[] intparams,float[] fltparams,Object[] stack,float[] outimage,int width,int height,int xc,int yc,int ch){
		this.intparams=intparams;
		this.fltparams=fltparams;
		this.outimage=outimage;
		this.width=width;
		this.height=height;
		this.stack=stack;
		this.xc=xc;
		this.yc=yc;
		this.ch=ch;
		random=new rngs();
	}
	
	public interp_line_v2(int[] intparams,float[] fltparams2,line3D[] lines,float[] distances,Object[] stack,float[] outimage,int width,int height,int yc,int ch,float[] depth){
		this.intparams=intparams;
		this.outimage=outimage;
		this.depth=depth;
		this.width=width;
		this.height=height;
		this.stack=stack;
		//this.xc=xc;
		this.yc=yc;
		this.ch=ch;
		this.lines=lines;
		this.distances=distances;
		//this.fltparams=fltparams;
		fltparams=new float[8];
		fltparams[6]=fltparams2[0];
		fltparams[7]=fltparams2[1];
		random=new rngs();
	}
	
	public void run(){
		int index2=yc*intparams[0];
		for(int i=0;i<intparams[0];i+=intparams[4]){
			xc=i;
			int index=index2+xc;
			fltparams[0]=lines[index].pt1.rx;
			fltparams[1]=lines[index].pt1.ry;
			fltparams[2]=lines[index].pt1.rz;
			fltparams[3]=(lines[index].pt2.rx-lines[index].pt1.rx)/distances[index];
			fltparams[4]=(lines[index].pt2.ry-lines[index].pt1.ry)/distances[index];
			fltparams[5]=(lines[index].pt2.rz-lines[index].pt1.rz)/distances[index];
    		float[] putval=get_line_max2(fltparams[0],fltparams[1],fltparams[2],fltparams[3],fltparams[4],fltparams[5],intparams[0],ch);
    		for(int j=yc;j<(yc+intparams[4]);j++){
    			for(int k=xc;k<(xc+intparams[4]);k++){
    				outimage[k+j*intparams[0]]=putval[0];
    				if(intparams[8]==3) depth[k+j*intparams[0]]=putval[1];
    			}
    		}
		}
	}
	
	public float get_3D_nearest_value2(float x,float y,float z,int channel){
		int x2=Math.round(x)-intparams[5];
		int y2=Math.round(y)-intparams[6];
		float scaledz=z-intparams[7];
		int z2=Math.round(scaledz);
		if(x2<0) return 0.0f;
		if(y2<0) return 0.0f;
		if(z2<0) return 0.0f;
		if(x2>=width) return 0.0f;
		if(y2>=height) return 0.0f;
		if(z2>=intparams[2]) return 0.0f;
		if(intparams[11]==0){
			return get_value2(x2+y2*width,z2,channel);
			//return 100.0f;
		} else if(intparams[11]==1) {
			return randinterpz2(x2+y2*width,scaledz,channel);
		} else {
			return interpolatez2(x2+y2*width,scaledz,channel);
		}
	}
	//these are 0maxsize,1channels,2slices,3frames,4xystep,5xrem,6yrem,7zremcorr,8method,9currframe,10imagetype,11projmeth
	public float interpolatez2(int index,float scaledz,int channel){
		int zprev=(int)scaledz;
		if(scaledz==zprev) return get_value2(index,zprev,channel);
		if(zprev<0) return 0.0f;
		int znext=zprev+1;
		if(znext>=intparams[2]) return 0.0f;
		float zrem2=scaledz-zprev;
		float z1=get_value2(index,zprev,channel);
		float z2=get_value2(index,znext,channel);
		return z1+zrem2*(z2-z1);
	}
	
	public float randinterpz2(int index,float scaledz,int channel){
		int zprev=(int)scaledz;
		if(zprev<0) return 0.0f;
		if(scaledz==zprev) return get_value2(index,zprev,channel);
		if(zprev>=(intparams[2]-1)) return 0.0f;
		float zrem2=scaledz-zprev;
		if((float)random.unidev(1.0,0.0)<zrem2) return get_value2(index,zprev,channel);
		else return get_value2(index,zprev+1,channel);
	}
	
	public float get_value2(int index,int z,int channel){
		int index2=intparams[9]*intparams[2]*intparams[1]+z*intparams[1]+channel;
		return get_value2(index,index2);
		//return index;
	}

	public float get_value2(int index,int index2){
		if(intparams[10]==0){
			//return (float)(((byte[])stack.getPixels(index2+1))[index]&0xff);
			return ((byte[])stack[index2])[index]&0xff;
		}else{
			if(intparams[10]==1){
				//return (float)(((short[])stack.getPixels(index2+1))[index]&0xffff);
				return ((short[])stack[index2])[index]&0xffff;
			}else{
				return ((float[])stack[index2])[index];
				//return ((float[])stack.getPixels(index2+1))[index];
			}
		}
	}
	//these are 0maxsize,1channels,2slices,3frames,4xystep,5xrem,6yrem,7zremcorr,8method,9currframe,10imagetype,11projmeth
	public float[] get_line_max2(float xstart,float ystart,float zstart,float xinc,float yinc,float zinc,int length,int channel){
		float x=xstart-xinc;
		float y=ystart-yinc;
		float scaledzinc=zinc/fltparams[6];
		float z=zstart/fltparams[6]-scaledzinc;
		// here xinc, yinc, and zinc are the dimensional steps corresponding to a pixel unit
		float max=0.0f; float maxpos=Float.NaN;
		for(int k=0;k<length;k++){
			x+=xinc;
			y+=yinc;
			z+=scaledzinc;
			if(x<intparams[5]+intparams[12]) continue;
			if(y<intparams[6]+intparams[14]) continue;
			if(z<intparams[7]+intparams[16]) continue;
			if(x>=intparams[5]+intparams[13]) continue;
			if(y>=intparams[6]+intparams[15]) continue;
			if(z>=intparams[7]+intparams[17]) continue;
			float temp=get_3D_nearest_value2(x,y,z,channel);
			if(temp<fltparams[7]) temp=0.0f;
			if(intparams[8]==0){ //max proj
				if(temp>max) {max=temp; maxpos=k;}
			} else if(intparams[8]==1 || intparams[8]==3){ //first proj and surface
				max=temp;
				if(temp>0.0f) return new float[]{max,(float)k};
			} else { //avg proj
				max+=temp;
			}
		}
		if(intparams[8]==2) max/=length;
		return new float[]{max,maxpos};
	}
	
}