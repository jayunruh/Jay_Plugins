/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.gui_interface;

public class tricrosscorr2{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;
	public gui_interface gi;

	/*
	 * In this version we output a 2D vector with different combinations of
	 * delays for two of the channels Copyright Jay Unruh Stowers Institute for
	 * Medical Research 4/25/08
	 */
	public tricrosscorr2(int length1,gui_interface gi){
		length=length1;
		fft=new po4realfft(length);
		this.gi=gi;
	}

	public float[][] docrosscorr(float[] data,float[] data2,float[] data3){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avg3=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			real2[j]=data2[j];
			real3[j]=data3[j];
			avg+=(double)data[j]/(double)length;
			avg2+=(double)data2[j]/(double)length;
			avg3+=(double)data3[j]/(double)length;
		}
		for(int j=0;j<length;j++){
			real[j]-=avg;
			real2[j]-=avg2;
			real3[j]-=avg3;
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		for(int i=0;i<length;i++){
			float a=real2[i];
			float b=im2[i];
			float c=real3[i];
			float d=im3[i];
			real2[i]=a*c+b*d; // here is the cross correlation of 2 and
			// 3
			im2[i]=b*c-a*d;
		}
		float[][] retarray=new float[length][];
		binmultilog bml=new binmultilog();
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			// correlate shifted versions of the 23 cross corr with 1
			// shifting in fourier space is equivalent to phase rotation
			for(int j=0;j<length;j++){
				float cosval=(float)Math.cos(j*2.0*Math.PI/length);
				float sinval=(float)Math.sin(j*2.0*Math.PI/length);
				float c=cosval*real2[j]-sinval*im2[j];
				float d=cosval*real2[j]+sinval*im2[j];
				real3[j]=a*c+b*d;
				im3[j]=b*c-a*d;
			}
			fft.realfft(real3,im3,true);
			for(int j=0;j<length;j++){
				real3[j]/=avg*avg2*avg3*length;
			}
			retarray[i]=bml.dobinmultilog(real3,length/2);
			gi.showProgress(i,length);
		}
		// now multilog bin along the other direction, transposing the vector
		int binlength=retarray[0].length;
		float[][] retarray2=new float[binlength][];
		for(int i=0;i<binlength;i++){
			float[] temp=new float[length];
			for(int j=0;j<length;j++)
				temp[j]=retarray[j][i];
			retarray2[i]=bml.dobinmultilog(temp,length/2);
		}
		return retarray2;
	}
}
