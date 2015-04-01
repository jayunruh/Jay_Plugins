/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class autocorr2D{
	public po4realfft2D fft; // this is the 2D fft used to calculate the
	// autocorrelation
	public int width,height;

	/*
	 * this class does a 2D spatial autocorrelation of an image the image must
	 * be a power of 2 size in width and height the image must be real once the
	 * class is constructed, it can be used to do any 2D autocorrelations on
	 * images of the same size Copyright Jay Unruh Stowers Institute for Medical
	 * Research 4/25/08
	 */

	public autocorr2D(int width1,int height1){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height);
	}

	public autocorr2D(int width1,int height1,int fftindex1,int fftindex2){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height,fftindex1,fftindex2);
	}

	public float[] doautocorr2D(float[] data,boolean doshiftxcenter,boolean doshiftycenter){
		return doautocorr2D(data,doshiftxcenter,doshiftycenter,true,false);
	}

	public float[] doautocorr2D(float[] data,boolean doshiftxcenter,boolean doshiftycenter,boolean doavgquad,boolean brightcorr){
		// do shiftxcenter shifts the data so that the zero point is in the x
		// center of the image
		float[] real=new float[width*height];
		float[] im=new float[width*height];
		double avg=0.0;
		for(int j=0;j<width*height;j++){
			real[j]=data[j];
			avg+=data[j];
		}
		avg/=((double)width*(double)height);
		fft.dorealfft2D(real,im,false);
		for(int j=0;j<width*height;j++){
			float a=real[j];
			float b=im[j];
			real[j]=a*a+b*b;
			im[j]=0.0f;
		}
		fft.dorealfft2D(real,im,true);
		for(int j=0;j<width*height;j++){
			real[j]/=(float)(avg*avg)*(width*height);
			real[j]-=1.0f;
			if(brightcorr){
				real[j]*=(float)avg;
			}
		}
		manipulate_quads mq=new manipulate_quads();
		if(doavgquad){
			real=mq.avgquadrants(real,width,height);
		}
		if(doshiftxcenter&&doshiftycenter){
			real=mq.shiftxycenter(real,width,height);
		}else{
			if(doshiftxcenter){
				real=mq.shiftxcenter(real,width,height);
			}
			if(doshiftycenter){
				real=mq.shiftycenter(real,width,height);
			}
		}
		return real;
	}

	public float[] doautocorr2D(short[] data,boolean doshiftxcenter,boolean doshiftycenter){
		return doautocorr2D(data,doshiftxcenter,doshiftycenter,true,false);
	}

	public float[] doautocorr2D(short[] data,boolean doshiftxcenter,boolean doshiftycenter,boolean doavgquad,boolean brightcorr){
		// do shiftxcenter shifts the data so that the zero point is in the x
		// center of the image
		float[] real=new float[width*height];
		float[] im=new float[width*height];
		double avg=0.0;
		for(int j=0;j<width*height;j++){
			real[j]=data[j]&0xffff;
			avg+=real[j];
		}
		avg/=((double)width*(double)height);
		fft.dorealfft2D(real,im,false);
		for(int j=0;j<width*height;j++){
			float a=real[j];
			float b=im[j];
			real[j]=a*a+b*b;
			im[j]=0.0f;
		}
		fft.dorealfft2D(real,im,true);
		for(int j=0;j<width*height;j++){
			real[j]/=(float)(avg*avg)*(width*height);
			real[j]-=1.0f;
			if(brightcorr){
				real[j]*=(float)avg;
			}
		}
		manipulate_quads mq=new manipulate_quads();
		if(doavgquad){
			real=mq.avgquadrants(real,width,height);
		}
		if(doshiftxcenter&&doshiftycenter){
			real=mq.shiftxycenter(real,width,height);
		}else{
			if(doshiftxcenter){
				real=mq.shiftxcenter(real,width,height);
			}
			if(doshiftycenter){
				real=mq.shiftycenter(real,width,height);
			}
		}
		return real;
	}

}
