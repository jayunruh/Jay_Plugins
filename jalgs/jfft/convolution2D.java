/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class convolution2D{
	public po4realfft2D fft; // this is the 2D fft class used for the
	// calculations
	public int width,height;

	/*
	 * this class does a 2D spatial crosscorrelation of two images the images
	 * must be a power of 2 size in width and height the image must be real once
	 * the class is constructed, it can be used to do any 2D crosscorrelations
	 * on images of the same size Copyright Jay Unruh Stowers Institute for
	 * Medical Research 4/25/08
	 */

	public convolution2D(int width1,int height1){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height);
	}
	
	public convolution2D(int width1,int height1,int fftindex1,int fftindex2){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height,fftindex1,fftindex2);
	}

	public float[] convolve2D(float[] data1,float[] data2){
		float[] real1=new float[width*height];
		float[] im1=new float[width*height];
		float[] real2=new float[width*height];
		float[] im2=new float[width*height];
		for(int j=0;j<width*height;j++){
			real1[j]=data1[j];
			real2[j]=data2[j];
		}
		fft.dorealfft2D(real1,im1,false);
		fft.dorealfft2D(real2,im2,false);
		for(int j=0;j<width*height;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c-b*d;
			im1[j]=b*c+a*d;
		}
		real2=null;
		im2=null;
		fft.dorealfft2D(real1,im1,true);
		return real1;
	}

	public float[] convolve2D(float[] data1,float[] rdata2,float[] idata2){
		// here we already have the fft of one of the data sets
		float[] real1=new float[width*height];
		float[] im1=new float[width*height];
		System.arraycopy(data1,0,real1,0,width*height);
		fft.dorealfft2D(real1,im1,false);
		for(int j=0;j<width*height;j++){
			float a=real1[j];
			float b=im1[j];
			float c=rdata2[j];
			float d=0.0f;
			if(idata2!=null)
				d=idata2[j];
			real1[j]=a*c-b*d;
			im1[j]=b*c+a*d;
		}
		fft.dorealfft2D(real1,im1,true);
		return real1;
	}

}
