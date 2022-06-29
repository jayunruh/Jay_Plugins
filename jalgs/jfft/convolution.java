/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class convolution{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public convolution(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[] convolve(float[] data1,float[] data2){
		float[] real1=new float[length];
		float[] im1=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		for(int j=0;j<length;j++){
			real1[j]=data1[j];
			real2[j]=data2[j];
		}
		fft.realfft(real1,im1,false);
		fft.realfft(real2,im2,false);
		for(int j=0;j<length;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c-b*d;
			im1[j]=b*c+a*d;
		}
		fft.realfft(real1,im1,true);
		return real1;
	}

}
