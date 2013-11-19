/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.kstats;

public class tetraautocorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public tetraautocorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[] doautocorr(float[] data){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		float[] real4=new float[length];
		float[] im4=new float[length];
		double avg=0.0;
		double avgsq=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			real2[j]=data[j]*data[j];
			avg+=(double)data[j]/(double)length;
			avgsq+=((double)data[j]*(double)data[j])/(double)length;
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
		float var=(float)(avgsq-avg*avg);
		for(int j=0;j<length;j++){
			float u22=real4[j]/(float)length-(2.0f*real[j]*(float)avg)/(float)length-(2.0f*real2[j]*(float)avg)/(float)length+2.0f*(float)(avgsq*avg*avg)+(4.0f*real3[j]*(float)(avg*avg))
					/(float)length-3.0f*(float)(avg*avg*avg*avg);
			float u21=real2[j]/(float)length-(2.0f*real3[j]*(float)avg)/(float)length-(float)(avgsq*avg)+2.0f*(float)(avg*avg*avg);
			float u12=real[j]/(float)length-(2.0f*real3[j]*(float)avg)/(float)length-(float)(avgsq*avg)+2.0f*(float)(avg*avg*avg);
			float u11=real3[j]/(float)length-(float)(avg*avg);
			real[j]=(u22-var*var-2.0f*u11*u11-u21-u12+u11)/(float)(avg*avg*avg*avg);
		}
		return real;
	}

	public float[][] doautocorr(float[] data,boolean dobrightcorr){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		float[] real4=new float[length];
		float[] im4=new float[length];
		double avg=0.0;
		double avgsq=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			real2[j]=data[j]*data[j];
			avg+=(double)data[j]/(double)length;
			avgsq+=((double)data[j]*(double)data[j])/(double)length;
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
		float var=(float)(avgsq-avg*avg);
		for(int j=0;j<length;j++){
			float u22=real4[j]/(float)length-(2.0f*real[j]*(float)avg)/(float)length-(2.0f*real2[j]*(float)avg)/(float)length+2.0f*(float)(avgsq*avg*avg)+(4.0f*real3[j]*(float)(avg*avg))
					/(float)length-3.0f*(float)(avg*avg*avg*avg);
			float u21=real2[j]/(float)length-(2.0f*real3[j]*(float)avg)/(float)length-(float)(avgsq*avg)+2.0f*(float)(avg*avg*avg);
			float u12=real[j]/(float)length-(2.0f*real3[j]*(float)avg)/(float)length-(float)(avgsq*avg)+2.0f*(float)(avg*avg*avg);
			float u11=real3[j]/(float)length-(float)(avg*avg);
			real[j]=(u22-var*var-2.0f*u11*u11-u21-u12+u11)/(float)(avg*avg*avg*avg);
			if(dobrightcorr){
				real[j]*=(float)(avg*avg*avg);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[1];
		retarray[1][0]=(float)avg;
		return retarray;
	}

	public float[][] doanalogautocorr(float[] data,boolean dobrightcorr,double Gain){
		double W=2.0*Gain;
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
		kstats ksclass=new kstats();
		double[] u1=ksclass.rawmoments(data);
		double[][] u=new double[9][9];
		for(int i=1;i<9;i++){
			u[i][0]=u1[i];
			u[0][i]=u1[i];
		}
		for(int j=0;j<length;j++){
			u[1][1]=(double)real3[j]/(double)length;
			u[2][1]=(double)real[j]/(double)length;
			u[1][2]=(double)real2[j]/(double)length;
			u[2][2]=(double)real4[j]/(double)length;
			double[][] k=ksclass.kstatistics(u,length);
			double afcum11=k[1][1];
			double afcum21=k[2][1]-afcum11*W;
			double afcum12=k[1][2]-afcum11*W;
			real[j]=(float)(k[2][2]-afcum21*W-afcum12*W-afcum11*W*W);
			real[j]/=(float)(k[1][0]*k[1][0]*k[1][0]*k[1][0]);
			if(dobrightcorr){
				real[j]*=(float)(k[1][0]*k[1][0]*k[1][0]);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[1];
		retarray[1][0]=(float)u[1][0];
		return retarray;
	}
}
