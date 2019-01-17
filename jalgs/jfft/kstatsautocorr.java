/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.kstats;

public class kstatsautocorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public kstatsautocorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[][] doautocorr(float[] data){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		float[] real4=new float[length];
		float[] im4=new float[length];
		for(int j=0;j<length;j++){
			real[j]=data[j];
			real2[j]=data[j]*data[j];
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			real[i]=a*c+b*d;
			im[i]=b*c-a*d;
			real2[i]=real[i];
			im2[i]=-im[i];
			real3[i]=a*a+b*b;
			im3[i]=0.0f;
			real4[i]=c*c+d*d;
			im4[i]=0.0f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		fft.realfft(real4,im4,true);
		// we will calculate all temporal bivariate k statistics up to order
		// four
		kstats ksclass=new kstats();
		double[] u1=ksclass.rawmoments(data);
		double[][] u=new double[9][9];
		double n=length;
		for(int i=1;i<=8;i++){
			u[i][0]=u1[i];
			u[0][i]=u1[i];
		}
		for(int i=0;i<length;i++){
			u[1][1]=real3[i]/n;
			u[2][1]=real[i]/n;
			u[1][2]=real2[i]/n;
			u[2][2]=real4[i]/n;
			double[][] tempk=ksclass.kstatistics(u,length);
			real[i]=(float)tempk[1][1];
			real2[i]=(float)tempk[2][1];
			real3[i]=(float)tempk[1][2];
			real4[i]=(float)tempk[2][2];
		}
		float[][] retvals=new float[5][];
		retvals[0]=real;
		retvals[1]=real2;
		retvals[2]=real3;
		retvals[3]=real4;
		retvals[4]=new float[4];
		double[] tempk2=ksclass.kstatistics(u1,length);
		float[] kout={(float)tempk2[1],(float)tempk2[2],(float)tempk2[3],(float)tempk2[4]};
		retvals[4]=kout;
		return retvals;
	}

}
