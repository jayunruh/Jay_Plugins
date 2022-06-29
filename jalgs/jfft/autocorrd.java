/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class autocorrd{
	public po4realfftd fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public autocorrd(int length1){
		length=length1;
		fft=new po4realfftd(length);
	}

	public double[][] doautocorr(double[] data,boolean dobrightcorr){
		double[] real=new double[length];
		double[] im=new double[length];
		double avg=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			avg+=data[j]/length;
		}
		fft.realfft(real,im,false);
		for(int i=0;i<length;i++){
			double a=real[i];
			double b=im[i];
			real[i]=a*a+b*b;
			im[i]=0.0;
		}
		fft.realfft(real,im,true);
		for(int j=0;j<length;j++){
			real[j]/=(avg*avg)*length;
			real[j]-=1.0;
			if(dobrightcorr){
				real[j]*=avg;
			}
		}
		double[][] retarray=new double[2][];
		retarray[0]=real;
		retarray[1]=new double[1];
		retarray[1][0]=avg;
		return retarray;
	}
}
