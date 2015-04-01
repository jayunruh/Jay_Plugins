/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.jstatistics;

public class gray_processing{

	public static float[] dilate(float[] image,int width,int height){
		float[] image2=new float[image.length];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float[] neighborhood=algutils.getNeighbors2(image,j,i,width,height);
				image2[j+i*width]=jstatistics.getstatistic("Max",neighborhood,null);
			}
		}
		return image2;
	}

	public static float[] erode(float[] image,int width,int height){
		float[] image2=new float[image.length];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float[] neighborhood=algutils.getNeighbors2(image,j,i,width,height);
				image2[j+i*width]=jstatistics.getstatistic("Min",neighborhood,null);
			}
		}
		return image2;
	}

	public static float[] open(float[] image,int iterations,int width,int height){
		float[] temp=image.clone();
		for(int i=0;i<iterations;i++){
			temp=erode(temp,width,height);
		}
		for(int i=0;i<iterations;i++){
			temp=dilate(temp,width,height);
		}
		return temp;
	}

	public static float[] close(float[] image,int iterations,int width,int height){
		float[] temp=image.clone();
		for(int i=0;i<iterations;i++){
			temp=dilate(temp,width,height);
		}
		for(int i=0;i<iterations;i++){
			temp=erode(temp,width,height);
		}
		return temp;
	}

	public static float[] tophat_black(float[] image,int iterations,int width,int height){
		float[] closed=close(image,iterations,width,height);
		float[] output=new float[image.length];
		for(int i=0;i<image.length;i++){
			output[i]=closed[i]-image[i];
		}
		return output;
	}

	public static float[] tophat_white(float[] image,int iterations,int width,int height){
		float[] opened=open(image,iterations,width,height);
		float[] output=new float[image.length];
		for(int i=0;i<image.length;i++){
			output[i]=image[i]-opened[i];
		}
		return output;
	}

}
