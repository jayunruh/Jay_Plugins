/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;
import jalgs.jseg.*;
import jalgs.jfit.*;
import jguis.*;

public class edf_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nchan=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Gradient_Smooth_Radius",2.0,5,15,null);
		gd.addNumericField("Frames To Average",5,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float smstdev=(float)gd.getNextNumber();
		int navg=(int)gd.getNextNumber();
		//first get the max sobel positions
		jsobel js=new jsobel(width,height);
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nchan;j++){
				IJ.showStatus("analyzing frame "+(i+1)+", chan "+(j+1));
				Object[] zseries=jutils.get3DZSeries(stack,j,i,nframes,nslices,nchan);
				float[] edfimg=edf(zseries,width,height,smstdev,navg);
				retstack.addSlice("",edfimg);
			}
		}
		jutils.create_hyperstack("EDF",retstack,imp,nframes,1,nchan).show();
		//new ImagePlus("EDF",retstack).show();
	}

	/*public float[] edf(Object[] stack,int width,int height,float blurstdev){
		float[] maxval=new float[width*height];
		float[] prevmaxval=new float[width*height];
		int[] maxpos=new int[width*height];
		int[] prevmaxpos=new int[width*height];
		float[] maxint=new float[width*height];
		float[] prevmaxint=new float[width*height];
		gausfunc gf=new gausfunc();
		jsobel js=new jsobel(width,height);
		for(int k=stack.length-1;k>=0;k--){
			float[] fimage=algutils.convert_arr_float2(stack[k]);
			float[][] sobel=js.do_sobel(fimage);
			//optionally blur the sobel image
			if(blurstdev>0.0f) jsmooth.blur2D(sobel[0],blurstdev,width,height,gf);
			for(int l=0;l<width*height;l++){
				if(sobel[0][l]>=maxval[l]){
					prevmaxval[l]=maxval[l]; prevmaxpos[l]=maxpos[l]; prevmaxint[l]=maxint[l];
					maxval[l]=sobel[0][l]; maxpos[l]=k; maxint[l]=fimage[l];
				}
			}
			IJ.showProgress(stack.length-k-1,stack.length);
		}
		for(int i=0;i<width*height;i++){
			maxint[i]=(maxval[i]*maxint[i]+prevmaxval[i]*prevmaxint[i])/(maxval[i]+prevmaxval[i]);
		}
		return maxint;
	}*/

	public float[] edf(Object[] stack,int width,int height,float blurstdev,int navg){
		float[][] maxval=new float[navg][width*height];
		int[][] maxpos=new int[navg][width*height];
		float[][] maxint=new float[navg][width*height];
		gausfunc gf=new gausfunc();
		jsobel js=new jsobel(width,height);
		for(int k=stack.length-1;k>=0;k--){
			float[] fimage=algutils.convert_arr_float2(stack[k]);
			float[][] sobel=js.do_sobel(fimage);
			//optionally blur the sobel image
			if(blurstdev>0.0f) jsmooth.blur2D(sobel[0],blurstdev,width,height,gf);
			for(int l=0;l<width*height;l++){
				if(sobel[0][l]>=maxval[0][l]){
					for(int i=(navg-1);i>0;i--){
						maxval[i][l]=maxval[i-1][l]; maxpos[i][l]=maxpos[i-1][l]; maxint[i][l]=maxint[i-1][l];
					}
					maxval[0][l]=sobel[0][l]; maxpos[0][l]=k; maxint[0][l]=fimage[l];
				}
			}
			IJ.showProgress(stack.length-k-1,stack.length);
		}
		for(int i=0;i<width*height;i++){
			float sum=0.0f;
			float norm=0.0f;
			for(int j=0;j<navg;j++){
				sum+=maxval[j][i]*maxint[j][i];
				norm+=maxval[j][i];
			}
			maxint[0][i]=sum/norm;
		}
		return maxint[0];
	}

}
