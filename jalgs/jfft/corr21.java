/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class corr21{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public corr21(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[] docrosscorr(float[] data,float[] data2){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avggr=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j]*data2[j];
			real2[j]=data[j];
			real3[j]=data2[j];
			avg+=(double)data[j]/(double)length;
			avg2+=(double)data2[j]/(double)length;
			avggr+=((double)data[j]*(double)data2[j])/length;
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			float e=real3[i];
			float f=im3[i];
			real[i]=a*c+b*d;
			im[i]=b*c-a*d;
			real2[i]=c*c+d*d;
			im2[i]=0.0f;
			real3[i]=c*e+d*f;
			im3[i]=d*e-c*f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg2)*length;
			real[j]-=real2[j]/((float)(avg*avg)*length);
			real[j]-=real3[j]/((float)(avg*avg2)*length);
			real[j]-=avggr/(float)(avg*avg2);
			real[j]+=2.0f;
		}
		return real;
	}

	public float[][] docrosscorr(float[] data,float[] data2,int brightcorrindex){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avggr=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j]*data2[j];
			real2[j]=data[j];
			real3[j]=data2[j];
			avg+=(double)data[j]/(double)length;
			avg2+=(double)data2[j]/(double)length;
			avggr+=((double)data[j]*(double)data2[j])/length;
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			float e=real3[i];
			float f=im3[i];
			real[i]=a*c+b*d;
			im[i]=b*c-a*d;
			real2[i]=c*c+d*d;
			im2[i]=0.0f;
			real3[i]=c*e+d*f;
			im3[i]=d*e-c*f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg2)*length;
			real[j]-=real2[j]/((float)(avg*avg)*length);
			real[j]-=real3[j]/((float)(avg*avg2)*length);
			real[j]-=avggr/(float)(avg*avg2);
			real[j]+=2.0f;
			if(brightcorrindex==1){
				real[j]*=(float)(avg*avg2);
			}
			if(brightcorrindex==2){
				real[j]*=(float)(avg*avg);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[2];
		retarray[1][0]=(float)avg;
		retarray[1][1]=(float)avg2;
		return retarray;
	}

	public float[][] docrosscorr_padded(float[] data,float[] data2,int brightcorrindex){
		int minlength=data.length;
		if(data2.length<minlength){
			minlength=data2.length;
		}
		if(minlength>=length){
			return docrosscorr(data,data2,brightcorrindex);
		}
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avggr=0.0;
		for(int j=0;j<minlength;j++){
			real[j]=data[j]*data2[j];
			real2[j]=data[j];
			real3[j]=data2[j];
			avg+=(double)data[j]/(double)minlength;
			avg2+=(double)data2[j]/(double)minlength;
			avggr+=((double)data[j]*(double)data2[j])/minlength;
		}
		for(int j=minlength;j<length;j++){
			real2[j]=(float)avg;
			real3[j]=(float)avg2;
			real[j]=(float)(avg*avg2);
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			float e=real3[i];
			float f=im3[i];
			real[i]=a*c+b*d;
			im[i]=b*c-a*d;
			real2[i]=c*c+d*d;
			im2[i]=0.0f;
			real3[i]=c*e+d*f;
			im3[i]=d*e-c*f;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg*avg2)*length;
			real[j]-=real2[j]/((float)(avg*avg)*length);
			real[j]-=real3[j]/((float)(avg*avg2)*length);
			real[j]-=avggr/(float)(avg*avg2);
			real[j]+=2.0f;
			if(brightcorrindex==1){
				real[j]*=(float)(avg*avg2);
			}
			if(brightcorrindex==2){
				real[j]*=(float)(avg*avg);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[2];
		retarray[1][0]=(float)avg;
		retarray[1][1]=(float)avg2;
		return retarray;
	}
}
