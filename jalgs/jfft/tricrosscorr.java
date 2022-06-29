/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class tricrosscorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public tricrosscorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public float[] docrosscorr(float[] data,float[] data2,float[] data3){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		float[] real4=new float[length];
		float[] im4=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avg3=0.0;
		double avggy=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j]*data2[j];
			real2[j]=data[j];
			real3[j]=data2[j];
			real4[j]=data3[j];
			avg+=(double)data[j]/(double)length;
			avg2+=(double)data2[j]/(double)length;
			avg3+=(double)data3[j]/(double)length;
			avggy+=((double)data[j]*(double)data2[j])/length;
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		fft.realfft(real4,im4,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			float e=real3[i];
			float f=im3[i];
			float g=real4[i];
			float h=im4[i];
			real[i]=a*g+b*h;
			im[i]=b*g-a*h;
			real2[i]=c*g+d*h;
			im2[i]=d*g-c*h;
			real3[i]=e*g+f*h;
			im3[i]=f*g-e*h;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg2*avg3)*length;
			real[j]-=real2[j]/((float)(avg*avg3)*length);
			real[j]-=real3[j]/((float)(avg2*avg3)*length);
			real[j]-=avggy/(float)(avg*avg2);
			real[j]+=2.0f;
		}
		return real;
	}

	public float[][] docrosscorr(float[] data,float[] data2,float[] data3,int brightcorrindex){
		float[] real=new float[length];
		float[] im=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		float[] real3=new float[length];
		float[] im3=new float[length];
		float[] real4=new float[length];
		float[] im4=new float[length];
		double avg=0.0;
		double avg2=0.0;
		double avg3=0.0;
		double avggy=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j]*data2[j];
			real2[j]=data[j];
			real3[j]=data2[j];
			real4[j]=data3[j];
			avg+=(double)data[j]/(double)length;
			avg2+=(double)data2[j]/(double)length;
			avg3+=(double)data3[j]/(double)length;
			avggy+=((double)data[j]*(double)data2[j])/length;
		}
		fft.realfft(real,im,false);
		fft.realfft(real2,im2,false);
		fft.realfft(real3,im3,false);
		fft.realfft(real4,im4,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			float c=real2[i];
			float d=im2[i];
			float e=real3[i];
			float f=im3[i];
			float g=real4[i];
			float h=im4[i];
			real[i]=a*g+b*h;
			im[i]=b*g-a*h;
			real2[i]=c*g+d*h;
			im2[i]=d*g-c*h;
			real3[i]=e*g+f*h;
			im3[i]=f*g-e*h;
		}
		fft.realfft(real,im,true);
		fft.realfft(real2,im2,true);
		fft.realfft(real3,im3,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg2*avg3)*length;
			real[j]-=real2[j]/((float)(avg*avg3)*length);
			real[j]-=real3[j]/((float)(avg2*avg3)*length);
			real[j]-=avggy/(float)(avg*avg2);
			real[j]+=2.0f;
			if(brightcorrindex==1){
				real[j]*=(float)(avg*avg2);
			}
			if(brightcorrindex==2){
				real[j]*=(float)(avg2*avg3);
			}
			if(brightcorrindex==3){
				real[j]*=(float)(avg*avg3);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[3];
		retarray[1][0]=(float)avg;
		retarray[1][1]=(float)avg2;
		retarray[1][2]=(float)avg3;
		return retarray;
	}
}
