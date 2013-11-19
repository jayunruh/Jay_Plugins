/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.kstats;

public class triautocorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public triautocorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[] doautocorr(float[] data){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
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
			real2[i]=a*a+b*b;
			im2[i]=0.0f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg)*(float)length;
			real[j]-=(2.0f*real2[j])/((float)(avg*avg)*(float)length);
			real[j]-=avgsq/(float)(avg*avg);
			real[j]+=2.0f;
			real[j]-=real2[j]/((float)(avg*avg*avg)*(float)length);
			real[j]+=1.0f/(float)(avg);
		}
		return real;
	}

	public float[][] doautocorr(float[] data,boolean dobrightcorr){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
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
			real2[i]=a*a+b*b;
			im2[i]=0.0f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg)*(float)length;
			real[j]-=(2.0f*real2[j])/((float)(avg*avg)*(float)length);
			real[j]-=avgsq/(float)(avg*avg);
			real[j]+=2.0f;
			real[j]-=real2[j]/((float)(avg*avg*avg)*(float)length);
			real[j]+=1.0f/(float)(avg);
			if(dobrightcorr){
				real[j]*=(float)(avg*avg);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[1];
		retarray[1][0]=(float)avg;
		return retarray;
	}

	public float[][] doautocorr_padded(float[] data,boolean dobrightcorr){
		if(data.length>=length){
			return doautocorr(data,dobrightcorr);
		}
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		double avg=0.0;
		double avgsq=0.0;
		for(int j=0;j<data.length;j++){
			real[j]=data[j];
			real2[j]=data[j]*data[j];
			avg+=(double)data[j]/(double)data.length;
			avgsq+=((double)data[j]*(double)data[j])/(double)data.length;
		}
		for(int j=data.length;j<length;j++){
			real[j]=(float)avg;
			real2[j]=(float)(avg*avg);
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
			real2[i]=a*a+b*b;
			im2[i]=0.0f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg)*(float)length;
			real[j]-=(2.0f*real2[j])/((float)(avg*avg)*(float)length);
			real[j]-=avgsq/(float)(avg*avg);
			real[j]+=2.0f;
			real[j]-=real2[j]/((float)(avg*avg*avg)*(float)length);
			real[j]+=1.0f/(float)(avg);
			real[j]*=(float)length/(float)data.length;
			if(dobrightcorr){
				real[j]*=(float)(avg*avg);
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
			real2[i]=a*a+b*b;
			im2[i]=0.0f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		double n=(double)length;
		kstats ksclass=new kstats();
		double[] u1=ksclass.rawmomentsshort(data);
		double[][] u=new double[5][5];
		for(int i=1;i<5;i++){
			u[i][0]=u1[i];
			u[0][i]=u1[i];
		}
		for(int j=0;j<length;j++){
			u[1][1]=real2[j]/n;
			u[2][1]=real[j]/n;
			u[1][2]=real[j]/n;
			double[][] k=ksclass.kstatisticsshort(u,length);
			double afcum11=k[1][1];
			real[j]=(float)(k[2][1]-(float)W*afcum11);
			real[j]/=(float)(u[1][0]*u[1][0]*u[1][0]);
			if(dobrightcorr){
				real[j]*=(float)(u[1][0]*u[1][0]);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[1];
		retarray[1][0]=(float)u[1][0];
		return retarray;
	}
}
