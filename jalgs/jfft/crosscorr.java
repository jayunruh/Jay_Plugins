/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class crosscorr{
	public po4realfft fft; // this is the 1D real fft required for the
	// correlation
	public int length;

	/*
	 * this class does a crosscorrelation of two float vectors the vectors must
	 * be a power of 2 length once the class is constructed, it can be used to
	 * do any crosscorrelations on vectors of the same length Copyright Jay
	 * Unruh Stowers Institute for Medical Research 4/25/08
	 */
	public crosscorr(int length1){
		length=length1;
		fft=new po4realfft(length);
	}

	public crosscorr(int length1,boolean usefft){
		length=length1;
		if(usefft){
			fft=new po4realfft(length);
		}
	}

	public float[] docrosscorr(float[] data1,float[] data2){
		float[] real1=new float[length];
		float[] im1=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		double avg1=0.0;
		double avg2=0.0;
		for(int j=0;j<length;j++){
			real1[j]=data1[j];
			avg1+=data1[j]/(double)length;
			real2[j]=data2[j];
			avg2+=data2[j]/(double)length;
		}
		fft.realfft(real1,im1,false);
		fft.realfft(real2,im2,false);
		for(int j=0;j<length;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c+b*d;
			im1[j]=b*c-a*d;
		}
		fft.realfft(real1,im1,true);
		for(int j=0;j<length;j++){
			real1[j]/=(float)(avg1*avg2)*(float)length;
			real1[j]-=1.0f;
		}
		return real1;
	}

	public float[][] docrosscorr(float[] data1,float[] data2,boolean dobrightcorr){
		float[] real1=new float[length];
		float[] im1=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		double avg1=0.0;
		double avg2=0.0;
		for(int j=0;j<length;j++){
			real1[j]=data1[j];
			avg1+=data1[j]/(double)length;
			real2[j]=data2[j];
			avg2+=data2[j]/(double)length;
		}
		fft.realfft(real1,im1,false);
		fft.realfft(real2,im2,false);
		for(int j=0;j<length;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c+b*d;
			im1[j]=b*c-a*d;
		}
		fft.realfft(real1,im1,true);
		for(int j=0;j<length;j++){
			real1[j]/=(float)(avg1*avg2)*(float)length;
			real1[j]-=1.0f;
			if(dobrightcorr){
				real1[j]*=(float)Math.sqrt(avg1*avg2);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real1;
		retarray[1]=new float[2];
		retarray[1][0]=(float)avg1;
		retarray[1][1]=(float)avg2;
		return retarray;
	}

	public float[][] docrosscorr_padded(float[] data1,float[] data2,boolean dobrightcorr){
		int minlength=data1.length;
		if(data2.length<minlength){
			minlength=data2.length;
		}
		if(minlength>=length){
			return docrosscorr(data1,data2,dobrightcorr);
		}
		float[] real1=new float[length];
		float[] im1=new float[length];
		float[] real2=new float[length];
		float[] im2=new float[length];
		double avg1=0.0;
		double avg2=0.0;
		for(int j=0;j<minlength;j++){
			real1[j]=data1[j];
			avg1+=data1[j]/(double)minlength;
			real2[j]=data2[j];
			avg2+=data2[j]/(double)minlength;
		}
		for(int j=minlength;j<length;j++){
			real1[j]=(float)avg1;
			real2[j]=(float)avg2;
		}
		fft.realfft(real1,im1,false);
		fft.realfft(real2,im2,false);
		for(int j=0;j<length;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c+b*d;
			im1[j]=b*c-a*d;
		}
		fft.realfft(real1,im1,true);
		for(int j=0;j<length;j++){
			real1[j]/=(float)(avg1*avg2)*(float)length;
			real1[j]-=1.0f;
			real1[j]*=(float)length/(float)minlength;
			if(dobrightcorr){
				real1[j]*=(float)Math.sqrt(avg1*avg2);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=real1;
		retarray[1]=new float[2];
		retarray[1][0]=(float)avg1;
		retarray[1][1]=(float)avg2;
		return retarray;
	}

	public float[][] docrosscorrnofft(float[] data1,float[] data2,boolean dobrightcorr){
		// first calculate <I(t)*I(t+tau)> for every possible tau value
		float[] temp=new float[length];
		double avg1=0.0;
		double avg2=0.0;
		for(int i=0;i<length;i++){ // loop over the tau values
			double temp3=0.0;
			for(int j=0;j<length;j++){
				int temp2=j+i;
				if(temp2>=length)
					temp2-=length;
				temp3+=(double)data1[j]*(double)data2[temp2];
			}
			temp[i]=(float)(temp3/(double)length);
			avg1+=(double)data1[i];
			avg2+=(double)data2[i];
		}
		avg1/=(double)length;
		avg2/=(double)length;
		for(int i=0;i<length;i++){
			temp[i]/=(float)(avg1*avg2);
			temp[i]-=1.0f;
			if(dobrightcorr){
				temp[i]*=(float)Math.sqrt(avg1*avg2);
			}
		}
		float[][] retarray=new float[2][];
		retarray[0]=temp;
		retarray[1]=new float[2];
		retarray[1][0]=(float)avg1;
		retarray[1][1]=(float)avg2;
		return retarray;
	}
}
