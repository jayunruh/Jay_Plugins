/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jalgs.jseg.*;
import jalgs.*;
import jguis.*;

public class remove_dark_pixels_jru_v1 implements PlugIn {
	//this plugin takes saturated spots and their surrounding regions and replaces them with the
	//average of an even wider image or a fill value
	//in a multichannel image, any channel over threshold can disqualify a pixel

	public void run(String arg) {
		float thresh=0.0f;
		int border=2;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Border Width",border,0);
		gd.addNumericField("Saturation Threshold",thresh,5,15,null);
		gd.addCheckbox("Fill_With_Number",false);
		gd.addNumericField("Fill_Value",0.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		border=(int)gd.getNextNumber();
		thresh=(float)gd.getNextNumber();
		boolean fillval=gd.getNextBoolean();
		float fillnumb=(float)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nframes=imp.getNFrames();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		for(int i=0;i<nframes*nslices;i++){
			Object[] chanstack=jutils.get3DCSeries(stack,i,0,1,nframes*nslices,nchans);
			fiximage(chanstack,width,height,border,thresh,fillval,fillnumb);
			IJ.showProgress(i,nframes*nslices);
			if(IJ.escapePressed()){return;}
		}
		imp.updateAndDraw();
	}

	public void fiximage(Object[] chanstack,int width,int height,int border,float thresh,boolean fillval,float fillnumb){
		//start by segmenting the saturated regions
		float[] maxproj=algutils.get_stack_proj_stat("Min",chanstack,width,height,chanstack.length,null);
		findblobs3 fb=new findblobs3(width,height);
		byte[] temp=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(maxproj[i]<=thresh) temp[i]=(byte)255;
			if(Float.isNaN(maxproj[i])) temp[i]=(byte)255;
		}
		float[] objects=fb.dofindblobs(temp);
		//now dilate the border
		for(int i=0;i<border;i++) fb.dilateobjects(objects,false);
		int nobj=fb.nobjects;
		//new ImagePlus("Saturated Objects",new FloatProcessor(width,height,objects,null)).show();
		//now get the circ averages
		float[][] circ=new float[nobj][chanstack.length];
		if(fillval){
			//fill the replacement array with the fillnumb value for every channel and every object
			for(int i=0;i<chanstack.length;i++){
				for(int j=0;j<nobj;j++) circ[j][i]=fillnumb;
			}
		} else {
			float[] circmask=fb.get_circ(objects,1); //get a boundary with width 1 pixel
			for(int i=0;i<chanstack.length;i++){
				for(int j=0;j<nobj;j++) circ[j][i]=fb.get_object_stats(circmask,j+1,chanstack[i],"Avg");
			}
		}
		int type=algutils.get_array_type(chanstack[0]);
		Object[] typecirc=algutils.convert_array((Object[])circ,type);
		for(int i=0;i<width*height;i++){
			if(objects[i]>0.0f){
				algutils.set_stack_col(chanstack,width,height,i,chanstack.length,typecirc[(int)objects[i]-1]);
			}
		}
	}

	/*public void fix8bitimage(byte[] image,int width,int height,int border,float thresh){
		//start by segmenting the saturated regions
		findblobs3 fb=new findblobs3(width,height);
		float[] objects=fb.dofindblobs(image,thresh);
		//now dilate the border
		for(int i=0;i<border;i++) fb.dilateobjects(objects,false);
		//new ImagePlus("Saturated Objects",new FloatProcessor(width,height,objects,null)).show();
		//now get the circ
		float[] circ=fb.get_circ(objects,1);
		//get the averages
		int nobj=fb.nobjects;
		float[] avgs=new float[nobj];
		for(int i=0;i<nobj;i++) avgs[i]=fb.get_object_stats(circ,i+1,image,"Avg");
		for(int i=0;i<width*height;i++){
			if(objects[i]>0.0f){
				float fill=avgs[(int)objects[i]-1];
				if(fill>255.0f) fill=255.0f;
				if(fill<0.0f) fill=0.0f;
				float ifill=(int)fill;
				image[i]=(byte)ifill;
			}
		}
	}*/

}
