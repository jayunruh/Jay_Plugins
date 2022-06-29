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
import jalgs.algutils;
import jalgs.matrixsolve;
import jalgs.profiler;

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

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JPanel;

public class maxproj3D_panel_v3 extends JPanel implements ActionListener,ItemListener,ImageListener{
	//in this version the data cube is static while the viewer rotates
	//in v3, opencl is used to do the rendering

	private JButton rotright,rotleft,rotup,rotdown,rotclock,rotcounter,reset,makemovie,droppoint,delpoint,exportbutton;
	private JButton profilebutton,importbutton;
	private JCheckBox stereocheck;
	private Choice projchoice,roichoice;
	private Label threshlabel,xysteplabel,projlabel,limitslabel,roilabel;
	private TextField threshval,xystepval,xstartval,xendval,ystartval,yendval,zstartval,zendval;
	public ImagePlus imp3D,srcimp;
	public ImageStack stack,stack3D;
	public int width,height,slices,channels,frames,newzslices,maxsize,imagetype;
	public int currframe,currchan,xystep,nthreads,xstart1,xend1,ystart1,yend1,zstart1,zend1;
	public int xrem,yrem,zrem,zremcorr,projmeth,centerindex,topindex,selroi;
	public float zratio,xrot,yrot,zrot,invzr,halfzr;
	private boolean stereo;
	private float[] thresh;
	private float[][] images;
	private element3D[] elements;
	private line3D[] lines;
	private renderer r,r2;
	private opencl_maxproj omp;

	public static Frame launch_frame(maxproj3D_panel_v3 panel){
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
		setLayout(null);
		this.zratio=zratio;
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
		stereo=false;
		stereocheck=new JCheckBox("Stereo",false);
		stereocheck.setBounds(10,bheight+20,100,20);
		stereocheck.addActionListener(this);
		add(stereocheck);
		projlabel=new Label("Projection");
		projlabel.setBounds(200,bheight+20,70,20);
		add(projlabel);
		projchoice=new Choice();
		projchoice.add("Max");
		projchoice.add("Avg");
		projchoice.add("First");
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
		xend1=width-1;
		ystart1=0;
		yend1=height-1;
		zstart1=0;
		zend1=slices-1;
		
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
				lines[j+i*maxsize]=new line3D(j-xrem,i-yrem,-zrem,j-xrem,i-yrem,maxsize-zrem,Color.black);
				elements[j+i*maxsize]=lines[j+i*maxsize];
			}
		}
		r=new renderer(maxsize,maxsize);
		r.centerz=maxsize/2-zrem;
		r.centerx=maxsize/2-xrem;
		r.centery=maxsize/2-yrem;
		//r.centerz=maxsize/2;
		r2=new renderer(maxsize,maxsize);
		r2.centerz=maxsize/2-zrem;
		r2.centerx=maxsize/2-xrem;
		r2.centery=maxsize/2-yrem;
		element3D[] cube=new element3D[1];
		cube[0]=new rect3D(r2.centerx,r2.centery,r2.centerz,width,height,newzslices,Color.white);
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
		projmeth=0;
		imagetype=0;
		if(stack.getPixels(1) instanceof short[])
			imagetype=1;
		if(stack.getPixels(1) instanceof float[])
			imagetype=2;
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		r2.setrotation((int)xrot,(int)yrot,(int)zrot);
		topindex=(int)(0.5f*maxsize);
		centerindex=topindex*maxsize+topindex;
		omp=new opencl_maxproj(jutils.stack2array(stack),width,height,zratio,maxsize,channels,currchan,this.thresh[currchan],projmeth,lines,elements,r);
		ImagePlus.addImageListener(this);
		update_image();
		((CompositeImage)imp3D).resetDisplayRanges();
	}

	public void setVisible(boolean b){
		if(!b){
			ImagePlus.removeImageListener(this);
			omp.dispose();
		}
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
				coords[0][i]=((spot3D)r2.elements[i+1]).point.x;
				coords[1][i]=((spot3D)r2.elements[i+1]).point.y;
				coords[2][i]=((spot3D)r2.elements[i+1]).point.z;
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
			gd.addNumericField("End Extension (frac of len)",0.0,5,15,null);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			int linethickness=(int)gd.getNextNumber();
			float fracext=(float)gd.getNextNumber();
			float[][] coords=new float[3][r2.elements.length-1];
			for(int i=0;i<r2.elements.length-1;i++){
				coords[0][i]=((spot3D)r2.elements[i+1]).point.x;
				coords[1][i]=((spot3D)r2.elements[i+1]).point.y;
				coords[2][i]=((spot3D)r2.elements[i+1]).point.z;
			}
			if(fracext>0.0f) coords=extendEndSegments(coords,fracext);
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
			//float[] temp=get_line_max(lines[index],distances[index],currchan);
			omp.setChannel(currchan);
			omp.setThresh(thresh[currchan]);
			omp.createProjection();
			float[] depth=omp.getDepth();
			float temp=depth[index];
			float length=get_3D_dist(lines[index].pt1,lines[index].pt2);
			float xinc=(float)(lines[index].pt2.x-lines[index].pt1.x)/length;
			float yinc=(float)(lines[index].pt2.y-lines[index].pt1.y)/length;
			float zinc=(float)(lines[index].pt2.z-lines[index].pt1.z)/length;
			float[] temp2={xinc*temp+lines[index].pt1.x,yinc*temp+lines[index].pt1.y,zinc*temp+lines[index].pt1.z};
			spot3D spot=new spot3D((int)temp2[0],(int)temp2[1],(int)temp2[2],0,Color.white);
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
    				spot3D spot=new spot3D((int)x,(int)y,(int)z,0,Color.white);
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
    				spot3D spot=new spot3D(x,y,z,0,Color.white);
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
		stereo=stereocheck.isSelected();
		update_image();
	}
	
	public float[][] extendEndSegments(float[][] coords1,float fracext){
		float[][] coords=algutils.clone_multidim_array(coords1);
		int length=profiler.get3DPolygonLength(coords[0],coords[1],coords[2],false);
		float extlen=fracext*(float)length;
		if(fracext<=0.0f) return coords;
		float len1=profiler.get3DLength(new float[]{coords[0][0],coords[1][0],coords[2][0],coords[0][1],coords[1][1],coords[2][1]});
		float xinc=(coords[0][1]-coords[0][0])/len1;
		float yinc=(coords[1][1]-coords[1][0])/len1;
		float zinc=(coords[2][1]-coords[2][0])/len1;
		coords[0][0]-=xinc*extlen;
		coords[1][0]-=yinc*extlen;
		coords[2][0]-=zinc*extlen;
		int e=coords[0].length-1;
		float len2=profiler.get3DLength(new float[]{coords[0][e],coords[1][e],coords[2][e],coords[0][e-1],coords[1][e-1],coords[2][e-1]});
		xinc=(coords[0][e]-coords[0][e-1])/len2;
		yinc=(coords[1][e]-coords[1][e-1])/len2;
		zinc=(coords[2][e]-coords[2][e-1])/len2;
		coords[0][e]+=xinc*extlen;
		coords[1][e]+=yinc*extlen;
		coords[2][e]+=zinc*extlen;
		return coords;
	}

	public void itemStateChanged(ItemEvent e){
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
			r2.translate(xrem,yrem,0);
			int[] overpix=r2.renderrgb();
			r2.translate(-xrem,-yrem,0);
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

	public void update_image(){
		int[] overpix=null;
		if(!stereo){
			float[][] images2=maxprojimage();
			for(int i=0;i<channels;i++){
				System.arraycopy(images2[i],0,images[i],0,maxsize*maxsize);
			}
			r2.translate(xrem,yrem,0);
			overpix=r2.renderrgb();
			r2.translate(-xrem,-yrem,0);
		} else {
			//first the left eye
			r.setrotation((int)xrot,(int)yrot-6,(int)zrot);
			float[][] images2=maxprojimage();
			for(int i=0;i<channels;i++){
				for(int j=0;j<maxsize;j+=2){
					System.arraycopy(images2[i],j*maxsize,images[i],j*maxsize,maxsize);
				}
			}
			update_cube();
			r2.translate(xrem,yrem,0);
			overpix=r2.renderrgb();
			r2.translate(-xrem,-yrem,0);
			//now the right eye
			r.setrotation((int)xrot,(int)yrot+6,(int)zrot);
			images2=maxprojimage();
			for(int i=0;i<channels;i++){
				for(int j=1;j<maxsize;j+=2){
					System.arraycopy(images2[i],j*maxsize,images[i],j*maxsize,maxsize);
				}
			}
			update_cube();
			r2.translate(xrem,yrem,0);
			int[] overpix2=r2.renderrgb();
			r2.translate(-xrem,-yrem,0);
			for(int i=1;i<maxsize;i+=2){
				System.arraycopy(overpix2,i*maxsize,overpix,i*maxsize,maxsize);
			}
			//now reset everything
			r.setrotation((int)xrot,(int)yrot,(int)zrot);
			update_cube();
		}
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
		point3D pt=new point3D(x,y,z);
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
			/*r2.translate(xrem,yrem,0);
			int[] overpix=r2.renderrgb();
			r2.translate(-xrem,-yrem,0);
			ColorProcessor cp=new ColorProcessor(maxsize,maxsize,overpix);
			ImageRoi roi=new ImageRoi(0,0,cp);
			roi.setZeroTransparent(true);
			roi.setPosition(i+1);
			overlay.add(roi);*/
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

	public float[][] maxprojimage(){
		float[][] outimage=new float[channels][maxsize*maxsize];
		//these are 0maxsize,1channels,2slices,3frames,4xystep,5xrem,6yrem,7zremcorr,8projmeth,9currframe,10imagetype,11method
		//int[] intparams={maxsize,channels,slices,frames,xystep,xrem,yrem,zremcorr,projmeth,currframe,imagetype,method,xstart1,xend1,ystart1,yend1,zstart1,zend1};
		//Object[] stackarray=stack.getImageArray();
		//Object[] stackarray=jutils.stack2array(stack);
		//long tstarttime=starttime;
		//float[][] depth=new float[channels][maxsize*maxsize];
		omp.setProjType(projmeth);
		float[] limits={xstart1,xend1,ystart1,yend1,zstart1,zend1};
		omp.setLimits(limits);
		omp.setLines(lines);
		for(int ch=0;ch<channels;ch++){
			//float[] depth=new float[maxsize*maxsize];
			omp.setChannel(ch);
			omp.setThresh(thresh[ch]);
			outimage[ch]=omp.createProjection();
			if(IJ.escapePressed()) return outimage;
			//long temp=System.currentTimeMillis();
			//IJ.log("line "+i+" time = "+(temp-tstarttime)/1000.0f);
			//tstarttime=temp;
		}

		return outimage;
	}

}