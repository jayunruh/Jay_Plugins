/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class autocorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// autocorrelation
	public int length;

	/*
	 * this class does an autocorrelation of a float vector the vector must be a
	 * power of 2 length once the class is constructed, it can be used to do any
	 * autocorrelations on vectors of the same length Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */
	public autocorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public autocorr(int length1,boolean usefft){
		length=length1;
		if(usefft)
			fft=new po4realfft(length);
	}

	public autocorr(int length1,int fftindex){
		length=length1;
		fft=new po4realfft(length,fftindex);
	}

	public float[] doautocorr(float[] data){
		float[] real=new float[length];
		float[] im=new float[length];
		double avg=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			avg+=(double)data[j]/(double)length;
		}
		fft.realfft(real,im,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			real[i]=a*a+b*b;
			im[i]=0.0f;
		}
		fft.realfft(real,im,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg)*(float)length;
			real[j]-=1.0f;
		}
		return real;
	}

	public float[][] doautocorr(float[] data,boolean dobrightcorr){
		float[] real=new float[length];
		float[] im=new float[length];
		double avg=0.0;
		for(int j=0;j<length;j++){
			real[j]=data[j];
			avg+=(double)data[j]/(double)length;
		}
		fft.realfft(real,im,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			real[i]=a*a+b*b;
			im[i]=0.0f;
		}
		fft.realfft(real,im,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg)*(float)length;
			real[j]-=1.0f;
			if(dobrightcorr){
				real[j]*=(float)avg;
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
		double avg=0.0;
		for(int j=0;j<data.length;j++){
			real[j]=data[j];
			avg+=(double)data[j]/(double)data.length;
		}
		for(int j=data.length;j<length;j++){
			real[j]=(float)avg;
		}
		fft.realfft(real,im,false);
		for(int i=0;i<length;i++){
			float a=real[i];
			float b=im[i];
			real[i]=a*a+b*b;
			im[i]=0.0f;
		}
		fft.realfft(real,im,true);
		for(int j=0;j<length;j++){
			real[j]/=(float)(avg*avg)*(float)length;
			real[j]-=1.0f;
			real[j]*=(float)length/(float)data.length;
			if(dobrightcorr){
				real[j]*=(float)avg;
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real;
		retarray[1]=new float[1];
		retarray[1][0]=(float)avg;
		return retarray;
	}

	public float[][] doautocorrnofft(float[] data1,boolean dobrightcorr){
		// first calculate <I(t)*I(t+tau)> for every possible tau value
		float[] temp=new float[length];
		double avg1=0.0;
		for(int i=0;i<length;i++){ // loop over the tau values
			double temp3=0.0;
			for(int j=0;j<length;j++){
				int temp2=j+i;
				if(temp2>=length)
					temp2-=length;
				temp3+=(double)data1[j]*(double)data1[temp2];
			}
			temp[i]=(float)(temp3/(double)length);
			avg1+=(double)data1[i];
		}
		avg1/=(double)length;
		for(int i=0;i<length;i++){
			temp[i]/=(float)(avg1*avg1);
			temp[i]-=1.0f;
			if(dobrightcorr){
				temp[i]*=(float)Math.sqrt(avg1*avg1);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=temp;
		retarray[1]=new float[1];
		retarray[1][0]=(float)avg1;
		return retarray;
	}
}
