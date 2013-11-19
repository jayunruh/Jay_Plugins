/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class complexac{
	public po4fft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does a complex autocorrelation of a real and imaginary float
	 * vector the vector must be a power of 2 length once the class is
	 * constructed, it can be used to do any autocorrelations on vectors of the
	 * same length Copyright Jay Unruh Stowers Institute for Medical Research
	 * 4/25/08
	 */
	public complexac(int length1){
		length=length1;
		fft=new po4fft(length);
	}

	public void doautocorr(float[] real,float[] im){
		fft.dopo4fft(real,im,false,true);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			real[i]=a*a+b*b;
			im[i]=0.0f;
		}
		fft.dopo4fft(real,im,true,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)length;
			im[j]/=(float)length;
		}
	}
}
