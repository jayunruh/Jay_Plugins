/*******************************************************************************
 * Copyright (c) 2022 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import jalgs.algutils;
import jalgs.jstatistics;
import jalgs.matrixsolve;
import jalgs.jfit.linleastsquares;
import jalgs.jsim.rngs;

import java.awt.Button;
import java.awt.Choice;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class Hist2DWindow_v2 extends Panel implements ActionListener,ItemListener,MouseMotionListener,ImageListener{

	private Button quit_button,update_button,smooth_button,revert_button,prof_button;
	private Button unmix_button, ungeom_button,ascale_button;
	private TextField threshval;
	private Label threshlabel,smtypelabel;
	private Choice smtypechoice;
	public ImagePlus imp;
	public ImageStack datastack;
	public PlotWindow2DHist pw;
	public float[] smoothxpix,smoothypix,smoothzpix,xpix,ypix,dpix,zpix;
	public float[] threshxpix,threshypix,tdpix;
	public int[] threshxcoords,threshycoords,mask;
	public int width,height;
	public float thresh=0f;
	public float xavg,yavg;
	public String xlab="x";
	public String ylab="y";
	public boolean threshfirst=true;
	public int smoothtype=0;
	public static int PWIDTH=250;
	public static int PHEIGHT=600;

	public static Frame launch_frame(Hist2DWindow_v2 panel){
		final Frame f=new Frame("Dynamic 2D Histogram Overlay");
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
		panel.setBounds(10,40,PWIDTH-20,PHEIGHT-20);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(PWIDTH,PHEIGHT));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(int width,int height,float[] xpix,float[] ypix,float[] dpix,ImagePlus dataimp){
		//the arguments here are x and y pixes of a histogram and corresponding display pixels
		//for FLIM phasor these are G and S components and a projected image
		this.width=width;
		this.height=height;
		this.xpix=xpix;
		this.ypix=ypix;
		this.dpix=dpix;
		smoothxpix=xpix.clone();
		smoothypix=ypix.clone();
		this.datastack=dataimp.getStack();
		update_hist();
		update_mask();
		setLayout(null);
		threshlabel=new Label("Threshold:");
		threshval=new TextField(""+thresh);
		addTextField(threshval,threshlabel,10,30);
		int buttonx=10;
		int buttony=30+25;
		int buttons=30;
		smooth_button=new Button("Smooth");
		addButton(smooth_button,buttonx,buttony);
		revert_button=new Button("Revert");
		addButton(revert_button,buttonx,buttony+buttons);
		prof_button=new Button("Get_Profile");
		addButton(prof_button,buttonx,buttony+2*buttons);
		unmix_button=new Button("Unmix");
		addButton(unmix_button,buttonx,buttony+3*buttons);
		ungeom_button=new Button("Unmix2");
		addButton(ungeom_button,buttonx,buttony+4*buttons);
		smtypelabel=new Label("Sm. Type:");
		smtypechoice=new Choice();
		addChoice(smtypechoice,smtypelabel,new String[]{"Mean","Median","Bin"},buttonx,buttony+5*buttons);
		ascale_button=new Button("Autoscale Hist");
		addButton(ascale_button,buttonx,buttony+6*buttons);
		update_button=new Button("Update");
		addButton(update_button,buttonx,PHEIGHT-100-buttons);
		quit_button=new Button("Quit");
		addButton(quit_button,buttonx,PHEIGHT-100);
		update_mask();
	}
	
	public void addChoice(Choice choice,Label lab,String[] vals,int x,int y) {
		lab.setBounds(x,y,80,20);
		choice.setBounds(x+80,y,100,20);
		for(int i=0;i<vals.length;i++) {
			choice.addItem(vals[i]);
		}
		choice.select(0);
		choice.addItemListener(this);
		add(choice);
		add(lab);
	}
	
	public void addTextField(TextField tf,Label lab,int x,int y) {
		lab.setBounds(x,y,80,20);
		tf.setBounds(x+80,y,100,20);
		tf.addActionListener(this);
		add(tf);
		add(lab);
		return;
	}
	
	public void addButton(Button button,int x,int y) {
		button.setBounds(x,y,125,30);
		button.addActionListener(this);
		add(button);
		return;
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void endprog(){
		/*if(imp!=null)
			imp.getCanvas().removeMouseMotionListener(this);
		ImagePlus.removeImageListener(this);*/
		if(pw!=null) {
			pw.getCanvas().removeMouseMotionListener(this);
		}
		this.getParent().setVisible(false);
		setVisible(false);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==quit_button){
			endprog();
			return;
		}
		if(e.getSource()==update_button){
			update_mask();
		}
		if(e.getSource()==smooth_button) {
			smooth();
		}
		if(e.getSource()==revert_button) {
			revert();
		}
		if(e.getSource()==prof_button) {
			save_profile();
		}
		if(e.getSource()==unmix_button) {
			save_unmixed();
		}
		if(e.getSource()==ungeom_button) {
			unmix_geom();
		}
		if(e.getSource()==ascale_button) {
			pw.intautoscale();
		}
		thresh=Float.parseFloat(threshval.getText());
		update_hist();
		update_mask();
	}
	
	public void smooth() {
		if(smoothtype!=2){
			// here we use a median or gaussian smoothing method
			// the weighting is 2 for the pixel and 1 for its neighbors
			// the pixels around the edge are not smoothed
			int i,j,k,l;
			float dumflt;
			float[] pixelvals=new float[10];
			// first smooth the x image
			for(int m=0;m<1;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothxpix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothxpix[i*width+j+m*width*height];
						}
						if(smoothtype==1){
							smoothxpix[i*width+j+m*width*height]=median(pixelvals);
						}else{
							smoothxpix[i*width+j+m*width*height]=avg(pixelvals);
						}
					}
				}
			}
			// now smooth the y image
			for(int m=0;m<1;m++){
				for(i=1;i<(height-1);i++){
					for(j=1;j<(width-1);j++){
						for(k=0;k<3;k++){
							for(l=0;l<3;l++){
								pixelvals[3*k+l]=smoothypix[(i-(k-1))*width+j-(l-1)+m*width*height];
							}
							pixelvals[9]=smoothypix[i*width+j+m*width*height];
						}
						if(smoothtype==1){
							smoothypix[i*width+j+m*width*height]=median(pixelvals);
						}else{
							smoothypix[i*width+j+m*width*height]=avg(pixelvals);
						}
					}
				}
			}
			if(smoothzpix!=null){
				for(int m=0;m<1;m++){
					for(i=1;i<(height-1);i++){
						for(j=1;j<(width-1);j++){
							for(k=0;k<3;k++){
								for(l=0;l<3;l++){
									pixelvals[3*k+l]=smoothzpix[(i-(k-1))*width+j-(l-1)+m*width*height];
								}
								pixelvals[9]=smoothzpix[i*width+j+m*width*height];
							}
							if(smoothtype==1){
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
				for(int m=0;m<1;m++){
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
										if(dpix[j*binsize+l+(i*binsize+k)*width]<thresh){
											abovethresh=false;
										}
									}else{
										if(dpix[j*binsize+l+(i*binsize+k)*width+m*width*height]<thresh){
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
	}
	
	void revert(){
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
	
	void save_profile(){
		if(datastack!=null){
			int slices=1;
			int currslice=0;
			int dtype=algutils.get_array_type(datastack.getPixels(1));
			int stacklength=datastack.getSize()/slices;
			float[] decay=new float[stacklength];
			if(dtype==2){
				for(int i=0;i<width*height;i++){
					if(mask[i]!=0){
						for(int j=0;j<stacklength;j++){
							decay[j]+=((float[])datastack.getPixels(currslice*stacklength+j+1))[i];
						}
					}
				}
			} else if(dtype==1) {
				for(int i=0;i<width*height;i++){
					if(mask[i]!=0){
						for(int j=0;j<stacklength;j++){
							float temp=(float)(((short[])datastack.getPixels(currslice*stacklength+j+1))[i]&0xffff);
							decay[j]+=temp;
						}
					}
				}
			} else {
				for(int i=0;i<width*height;i++){
					if(mask[i]!=0){
						for(int j=0;j<stacklength;j++){
							float temp=(float)(((byte[])datastack.getPixels(currslice*stacklength+j+1))[i]&0xff);
							decay[j]+=temp;
						}
					}
				}
			}
			(new PlotWindow4("Masked Profile","channel","Intensity",decay)).draw();
		}
		//output the roi stats as well for the lifetime phasor
		//if(calltype==5) {
			GenericDialog taugd=new GenericDialog("Options");
			taugd.addNumericField("Frequency (MHz)",80.0,5,15,null);
			taugd.showDialog(); if(taugd.wasCanceled()) return;
			float megafreq=(float)taugd.getNextNumber();
			IJ.log("G avg= "+xavg);
			IJ.log("S avg= "+yavg);
			float tempfreq=2.0f*megafreq*0.001f*(float)Math.PI;
			float tphi=yavg/(xavg*tempfreq);
			float tmod=(float)Math.sqrt((1.0f/(xavg*xavg+yavg*yavg))-1.0f);
			tmod/=tempfreq;
			IJ.log("Tau Mod= "+tmod);
			IJ.log("Tau Phase= "+tphi);
		//}
	}
	
	void save_unmixed(){
		//here we need to gather and unmix each bin of the histogram
		if(datastack==null) return;
		//start by getting the spectra
		ImageWindow[] iw=jutils.selectPlots(false,1,new String[]{"Reference Spectra"});
		if(iw==null) return;
		float[][] spectra=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
		int slen=spectra[0].length;
		int bs=pw.getBinSize();
		float[] lims=pw.getLimits();
		int hs=Plot2DHist.histSize/bs;
		//need to build up a decay matrix for each point in the 2D histogram
		float[][] fractions=new float[hs*hs][spectra.length];
		float[] sums=new float[width*height];
		float[][] decaymatrix=new float[hs*hs][slen];
		int[] histindices=new int[width*height];
		int[] counts=new int[hs*hs];
		Object[] ds2=jutils.stack2array(datastack);
		int dtype=algutils.get_array_type(ds2[0]);
		//walk through the x and y images, adding them to the decay matrix as we go
		for(int i=0;i<threshxpix.length;i++) {
			//map this pixel in the histogram
			int posx=(int)((threshxpix[i]-lims[0])/(lims[1]-lims[0])*hs);
			int posy=(int)((threshypix[i]-lims[2])/(lims[3]-lims[2])*hs);
			if((posx>=0&&posx<hs)&&(posy>=0&&posy<hs)){
				counts[posx+posy*hs]++;
				Object decay1=algutils.get_stack_col(ds2, width, height, threshxcoords[i], threshycoords[i], slen);
				float[] decay=algutils.convert_arr_float(decay1);
				int histidx=posx+posy*hs;
				for(int j=0;j<slen;j++) decaymatrix[histidx][j]+=decay[j];
				int idx=threshxcoords[i]+threshycoords[i]*width;
				sums[idx]+=jstatistics.getstatistic("Sum", decay, null);
				histindices[idx]=histidx;
			}
		}
		//now unmix the decay matrix to get spectral fractions in each bin
		linleastsquares lls=new linleastsquares(spectra);
		for(int i=0;i<decaymatrix.length;i++){
			if(counts[i]>0){
				fractions[i]=algutils.convert_arr_float(lls.fitdata(decaymatrix[i],null));
				float sum=jstatistics.getstatistic("Sum",fractions[i],null);
				for(int j=0;j<fractions[i].length;j++) fractions[i][j]/=sum;
			}
		}
		float[][] unmixed=new float[spectra.length][width*height];
		//now go back to the image and calculate the contributions
		for(int i=0;i<threshxpix.length;i++){
			int idx=threshxcoords[i]+threshycoords[i]*width;
			if(sums[idx]>=0.0f){
				float[] contr=fractions[histindices[idx]];
				for(int j=0;j<contr.length;j++){
					unmixed[j][idx]=contr[j]*sums[idx];
				}
			}
		}
		ImagePlus unimp=new ImagePlus("Unmixed",jutils.array2stack(unmixed,width,height));
		unimp.setOpenAsHyperStack(true);
		unimp.setDimensions(unmixed.length,1,1);
		new CompositeImage(unimp,CompositeImage.COLOR).show();
	}
	
	void unmix_geom(){
		GenericDialog gd2=new GenericDialog("Number_of_components");
		gd2.addCheckbox("Use_ROI_Manager_Points?", false);
		gd2.addNumericField("Otherwise_How_Many_Components?",3,0);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		boolean userman=gd2.getNextBoolean();
		int ncomp=(int)gd2.getNextNumber();
		float[][] refpos=null;
		if(!userman) {
			GenericDialog gd=new GenericDialog("Options");
			for(int i=0;i<ncomp;i++){
				gd.addNumericField("G"+(i+1),0.0,5,15,null);
				gd.addNumericField("S"+(i+1),0.0,5,15,null);
			}
			gd.showDialog(); if(gd.wasCanceled()) return;
			refpos=new float[ncomp][2];
			for(int i=0;i<ncomp;i++){
				refpos[i][0]=(float)gd.getNextNumber();
				refpos[i][1]=(float)gd.getNextNumber();
			}
		} else {
			Roi[] rois=RoiManager.getInstance2().getRoisAsArray();
			ncomp=rois.length;
			refpos=new float[ncomp][2];
			for(int i=0;i<rois.length;i++) {
				Rectangle r=rois[i].getBounds();
				float[] coords=pw.getPlot().getPlotCoords(r.x,r.y);
				refpos[i][0]=coords[0];
				refpos[i][1]=coords[1];
			}
		}
		//this is the approach suggested by Francesco Cutrale
		//here our basis spectra are three reference phasor (G,S) positions
		//anything outside of the resulting triangle is masked out (or not)
		//inside the triangle, pixels are assumed to be linear combinations of the three spectra
		int bs=pw.getBinSize();
		float[] lims=pw.getLimits();
		int hs=Plot2DHist.histSize/bs;
		//need to build up a decay matrix for each point in the 2D histogram
		float[][] fractions=new float[hs*hs][ncomp];
		float[] sums=new float[width*height];
		int[] histindices=new int[width*height];
		int[] counts=new int[hs*hs];
		Object[] ds2=jutils.stack2array(datastack);
		int slen=ds2.length;
		int dtype=algutils.get_array_type(ds2[0]);
		//walk through the x and y images, adding them to the decay matrix as we go
		for(int i=0;i<threshxpix.length;i++) {
			//map this pixel in the histogram
			int posx=(int)((threshxpix[i]-lims[0])/(lims[1]-lims[0])*hs);
			int posy=(int)((threshypix[i]-lims[2])/(lims[3]-lims[2])*hs);
			if((posx>=0&&posx<hs)&&(posy>=0&&posy<hs)){
				counts[posx+posy*hs]++;
				Object decay1=algutils.get_stack_col(ds2, width, height, threshxcoords[i], threshycoords[i], slen);
				//float[] decay=algutils.convert_arr_float(decay1);
				int histidx=posx+posy*hs;
				//for(int j=0;j<slen;j++) decaymatrix[histidx][j]+=decay[j];
				int idx=threshxcoords[i]+threshycoords[i]*width;
				sums[idx]+=jstatistics.getstatistic("Sum", decay1, null);
				histindices[idx]=histidx;
			}
		}
		//now unmix all of the histogram pixels
		if(ncomp==3){
    		double[][] matrix={{refpos[0][0],refpos[1][0],refpos[2][0]},{refpos[0][1],refpos[1][1],refpos[2][1]},{1.0,1.0,1.0}};
    		double[][] minv=(new matrixsolve()).gjinv2(matrix,3);
    		for(int i=0;i<hs;i++){
    			for(int j=0;j<hs;j++){
    				if(counts[j+i*hs]>0){
    					float g=((lims[1]-lims[0])*(float)j)/hs+lims[0];
    					float s=((lims[3]-lims[2])*(float)i)/hs+lims[2];
    					fractions[j+i*hs][0]=g*(float)minv[0][0]+s*(float)minv[0][1]+(float)minv[0][2];
    					fractions[j+i*hs][1]=g*(float)minv[1][0]+s*(float)minv[1][1]+(float)minv[1][2];
    					fractions[j+i*hs][2]=g*(float)minv[2][0]+s*(float)minv[2][1]+(float)minv[2][2];
    				}
    			}
    		}
		} else {
			if(ncomp==2){
				for(int i=0;i<hs;i++){
	    			for(int j=0;j<hs;j++){
	    				if(counts[j+i*hs]>0){
	    					float g=((lims[1]-lims[0])*(float)j)/hs+lims[0];
	    					float s=((lims[3]-lims[2])*(float)i)/hs+lims[2];
	    					float fx=(g-refpos[1][0])/(refpos[0][0]-refpos[1][0]);
	    					float fy=(s-refpos[1][1])/(refpos[0][1]-refpos[1][1]);
	    					float fa=0.5f*(fx+fy);
	    					if(Math.abs(refpos[0][0]-refpos[1][0])<0.00001f) fa=fy;
	    					if(Math.abs(refpos[0][1]-refpos[1][1])<0.00001f) fa=fx;
	    					fractions[j+i*hs][0]=fa;
	    					fractions[j+i*hs][1]=1.0f-fa;
	    				}
	    			}
	    		}
			} else {
				//here use stochastic simulated fractions
				float[][] temp=sim_fractions(refpos,1000000,new float[]{lims[0],lims[1],lims[2],lims[3]},hs);
				for(int i=0;i<temp.length;i++){
					if(temp[i]!=null) fractions[i]=temp[i];
				}
			}
		}
		float[][] unmixed=new float[refpos.length*1][width*height];
		//now go back to the image and calculate the contributions
		for(int i=0;i<width*height;i++){
			if(sums[i]>=0.0f){
				float[] contr=fractions[histindices[i]];
				for(int j=0;j<contr.length;j++){
					unmixed[j][i]=contr[j]*sums[i];
				}
			}
		}
		ImagePlus unimp=new ImagePlus("Unmixed",jutils.array2stack(unmixed,width,height));
		unimp.setOpenAsHyperStack(true);
		unimp.setDimensions(unmixed.length,1,1);
		new CompositeImage(unimp,CompositeImage.COLOR).show();
		new ImagePlus("Fractions",jutils.array2stack(fractions,256,256)).show();
	}
	
	public static float[][] sim_fractions(float[][] points,int nsims,float[] lims,int hs){
		//here we simulate fractional phasor "images" for random combinations of points
		//the phasor goes from -1 to 1 in x and y (or lims) with hs x hs points
		int[] pix={hs,hs};
		float[][] fractions=new float[pix[0]*pix[1]][];
		//int[] counts=new int[pix[0]*pix[1]];
		rngs random=new rngs();
		for(int i=0;i<nsims;i++){
			float[] frac=random.random_fractions(points.length);
			float[] gs=fractions2gs(points,frac);
			int[] xy=gs2xy(gs,lims,pix);
			int index=xy[0]+xy[1]*pix[0];
			if(fractions[index]==null) fractions[index]=frac;
			else{
				for(int j=0;j<points.length;j++) fractions[index][j]+=frac[j];
			}
			//counts[index]++;
		}
		for(int i=0;i<fractions.length;i++){
			if(fractions[i]!=null){
				float sum=0.0f;
				for(int j=0;j<fractions[i].length;j++){sum+=fractions[i][j];}
				for(int j=0;j<fractions[i].length;j++){fractions[i][j]/=sum;}
			}
		}
		return fractions;
	}
	
	public static int[] gs2xy(float[] gs,float[] lims,int[] pix){
		int x=(int)(((gs[0]-lims[0])/(lims[1]-lims[0]))*(float)pix[0]);
		int y=(int)(((gs[1]-lims[2])/(lims[3]-lims[2]))*(float)pix[1]);
		return new int[]{x,y};
	}
	
	public static float[] fractions2gs(float[][] points,float[] fractions){
		float[] gs=new float[2];
		for(int i=0;i<points.length;i++){
			gs[0]+=fractions[i]*points[i][0];
			gs[1]+=fractions[i]*points[i][1];
		}
		return gs;
	}
	
	public void update_hist() {
		//this should happen whenever the smoothpix or the threshold are changed
		threshxcoords=new int[width*height]; //this holds the spatial coordinates
		threshycoords=new int[width*height];
		threshxpix=new float[width*height]; //this holds the actual phasor values
		threshypix=new float[width*height];
		tdpix=new float[width*height];
		int nthresh=0;
		for(int i=0;i<height;i++) {
			for(int j=0;j<width;j++) {
				if(dpix[j+i*width]>thresh) {
					tdpix[j+i*width]=dpix[j+i*width];
					threshxcoords[nthresh]=j;
					threshycoords[nthresh]=i;
					threshxpix[nthresh]=smoothxpix[j+i*width];
					threshypix[nthresh]=smoothypix[j+i*width];
					nthresh++;
				}
			}
		}
		threshxpix=(float[])algutils.get_subarray(threshxpix, 0, nthresh);
		threshypix=(float[])algutils.get_subarray(threshypix, 0, nthresh);
		if(pw==null) {
			pw=new PlotWindow2DHist("2D Histogram",xlab,ylab,threshxpix,threshypix,null);
			pw.draw();
			pw.getImagePlus().setRoi(new Rectangle(100,100,20,20));
			ImagePlus.addImageListener(this);
			pw.getCanvas().addMouseMotionListener(this);
		} else {
			pw.updateData(threshxpix, threshypix, false);
		}
	}

	public void update_mask(){
		//first get the roi
		//Roi roi=pw.getImagePlus().getRoi();
		//Plot2DHist plot=pw.getPlot();
		int[] indices=pw.getroiindices();
		if(indices==null) return;
		mask=new int[width*height];
		xavg=0.0f;
		yavg=0.0f;
		for(int i=0;i<indices.length;i++) {
			int xpos=threshxcoords[indices[i]];
			int ypos=threshycoords[indices[i]];
			mask[xpos+ypos*width]=0xffff0000;
			xavg+=threshxpix[indices[i]];
			yavg+=threshypix[indices[i]];
		}
		xavg/=(float)indices.length;
		yavg/=(float)indices.length;
		ImageRoi roi=new ImageRoi(0,0,new ColorProcessor(width,height,mask));
		roi.setZeroTransparent(true);
		if(imp==null) {
			imp=new ImagePlus("Display Image",new FloatProcessor(width,height,tdpix,null));
			imp.show();
			imp.setOverlay(new Overlay(roi));
		} else {
			imp.setProcessor(new FloatProcessor(width,height,tdpix,null));
			imp.updateAndDraw();
			imp.setOverlay(new Overlay(roi));
		}
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		update_mask();
	}

	public void mouseClicked(MouseEvent arg0){
		update_mask();
	}

	public void mouseEntered(MouseEvent arg0){
	}

	public void mouseExited(MouseEvent arg0){
	}

	public void mousePressed(MouseEvent arg0){
	}

	public void mouseReleased(MouseEvent arg0){
		update_mask();
	}

	public void imageClosed(ImagePlus arg0){
		if(arg0.equals(pw.getImagePlus())||arg0.equals(imp))
			endprog();
	}

	public void imageOpened(ImagePlus arg0){
	}

	public void imageUpdated(ImagePlus arg0){
		if(arg0.equals(pw.getImagePlus())){
			update_mask();
		}
	}

	public void itemStateChanged(ItemEvent e) {
		if(e.getSource()==smtypechoice) {
			smoothtype=smtypechoice.getSelectedIndex();
		}
	}

}
