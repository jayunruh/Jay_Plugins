/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class complexcc{
	public po4fft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does a crosscorrelation of two complex data sets (real and im
	 * float vectors) the vectors must be a power of 2 length once the class is
	 * constructed, it can be used to do any crosscorrelations on vectors of the
	 * same length Copyright Jay Unruh Stowers Institute for Medical Research
	 * 4/25/08
	 */
	public complexcc(int length1){
		length=length1;
		fft=new po4fft(length);
	}

	public void docrosscorr(float[] real1,float[] im1,float[] real2,float[] im2){
		fft.dopo4fft(real1,im1,false,true);
		fft.dopo4fft(real2,im2,false,true);
		for(int i=0;i<length;i++){
			float a=real1[i];
			float b=im1[i];
			float c=real2[i];
			float d=im2[i];
			real1[i]=a*c+b*d;
			im1[i]=b*c-a*d;
		}
		fft.dopo4fft(real1,im1,true,true);
		for(int j=0;j<length;j++){
			real1[j]/=length;
			im1[j]/=length;
		}
	}
}
