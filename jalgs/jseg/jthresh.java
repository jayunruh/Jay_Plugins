/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.interpolation;
import jalgs.jstatistics;

public class jthresh{

	public static byte[] autothresh(float[] image,int width,int height,float fraction){
		// here we autothreshold as a percentage of the distance between the
		// smoothed
		// minimum and the smoothed maximum
		float thresh=getautothresh(image,width,height,fraction);
		byte[] retimage=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(image[i]>=thresh){
				retimage[i]=(byte)255;
			}else{
				retimage[0]=(byte)0;
			}
		}
		return retimage;
	}

	public static byte[] autothresh(float[] image,int width,int height){
		return autothresh(image,width,height,0.5f);
	}

	public static float getautothresh(float[] image,int width,int height,float fraction){
		// here we autothreshold as a percentage of the distance between the
		// smoothed
		// minimum and the smoothed maximum
		float[] temp=jsmooth.smooth2D(image,width,height);
		float[] minmax=jstatistics.getminmax(temp);
		return fraction*(minmax[1]-minmax[0])+minmax[0];
	}

	public static float getautothresh2(float[] image,int width,int height,float fraction){
		// here we autothreshold as a percentage of the distance between the
		// smoothed
		// minimum and the smoothed maximum
		float[] minmax=jstatistics.getminmax(image);
		return fraction*(minmax[1]-minmax[0])+minmax[0];
	}

	public static byte[] adaptivethresh(float[] image,int width,int height,threshinterface tf,int size,int skip){
		float[] thresh2=getadaptivethresh(image,width,height,tf,size,skip);
		byte[] temp=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(image[i]>=thresh2[i])
				temp[i]=(byte)255;
		}
		return temp;
	}

	public static float[] getadaptivethresh(float[] image,int width,int height,threshinterface tf,int size,int skip){
		int newwidth=(int)((float)width/(float)skip);
		int newheight=(int)((float)height/(float)skip);
		float[] thresh=new float[newwidth*newheight];
		for(int i=0;i<newheight;i++){
			for(int j=0;j<newwidth;j++){
				float[] subregion=algutils.get_region_pad(image,j*skip+skip/2,i*skip+skip/2,size,size,width,height);
				thresh[j+i*newwidth]=tf.calcthresh(subregion);
			}
		}
		float[] temp=new float[width*height];
		for(int i=0;i<height;i++){
			float ypos=(float)i/(float)skip-0.5f;
			for(int j=0;j<width;j++){
				float xpos=(float)j/(float)skip-0.5f;
				temp[j+i*width]=interpolation.interp2D_pad(thresh,newwidth,newheight,xpos,ypos);
			}
		}
		return temp;
	}

}
