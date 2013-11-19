/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class crosscorr2D{
	public po4realfft2D fft; // this is the 2D fft class used for the
	// calculations
	public int width,height;

	/*
	 * this class does a 2D spatial crosscorrelation of two images the images
	 * must be a power of 2 size in width and height the image must be real once
	 * the class is constructed, it can be used to do any 2D crosscorrelations
	 * on images of the same size Copyright Jay Unruh Stowers Institute for
	 * Medical Research 4/25/08
	 */

	public crosscorr2D(int width1,int height1){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height);
	}

	public float[] docrosscorr2D(float[] data1,float[] data2,boolean doshiftxcenter,boolean doshiftycenter){
		return docrosscorr2D(data1,data2,doshiftxcenter,doshiftycenter,true,false);
	}

	public float[] docrosscorr2D(float[] data1,float[] data2,boolean doshiftxcenter,boolean doshiftycenter,boolean doavgquad,boolean brightcorr){
		float[] real1=new float[width*height];
		float[] im1=new float[width*height];
		float[] real2=new float[width*height];
		float[] im2=new float[width*height];
		double avg1=0.0;
		double avg2=0.0;
		for(int j=0;j<width*height;j++){
			real1[j]=data1[j];
			avg1+=data1[j];
			real2[j]=data2[j];
			avg2+=data2[j];
		}
		avg1/=((double)width*(double)height);
		avg2/=((double)width*(double)height);
		fft.dorealfft2D(real1,im1,false);
		fft.dorealfft2D(real2,im2,false);
		for(int j=0;j<width*height;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c+b*d;
			im1[j]=b*c-a*d;
		}
		real2=null;
		im2=null;
		fft.dorealfft2D(real1,im1,true);
		for(int j=0;j<width*height;j++){
			real1[j]/=(float)(avg1*avg2)*(float)(width*height);
			real1[j]-=1.0f;
			if(brightcorr){
				real1[j]*=(float)Math.sqrt(avg1*avg2);
			}
		}
		manipulate_quads mq=new manipulate_quads();
		if(doavgquad){
			real1=mq.avgquadrants(real1,width,height);
		}
		if(doshiftxcenter&&doshiftycenter){
			real1=mq.shiftxycenter(real1,width,height);
		}else{
			if(doshiftxcenter){
				real1=mq.shiftxcenter(real1,width,height);
			}
			if(doshiftycenter){
				real1=mq.shiftycenter(real1,width,height);
			}
		}
		return real1;
	}

	public float[] docrosscorr2D(short[] data1,short[] data2,boolean doshiftxcenter,boolean doshiftycenter){
		return docrosscorr2D(data1,data2,doshiftxcenter,doshiftycenter,true,false);
	}

	public float[] docrosscorr2D(short[] data1,short[] data2,boolean doshiftxcenter,boolean doshiftycenter,boolean doavgquad,boolean brightcorr){
		float[] real1=new float[width*height];
		float[] im1=new float[width*height];
		float[] real2=new float[width*height];
		float[] im2=new float[width*height];
		double avg1=0.0;
		double avg2=0.0;
		for(int j=0;j<width*height;j++){
			real1[j]=data1[j]&0xffff;
			avg1+=real1[j];
			real2[j]=data2[j]&0xffff;
			avg2+=real2[j];
		}
		avg1/=((double)width*(double)height);
		avg2/=((double)width*(double)height);
		fft.dorealfft2D(real1,im1,false);
		fft.dorealfft2D(real2,im2,false);
		for(int j=0;j<width*height;j++){
			float a=real1[j];
			float b=im1[j];
			float c=real2[j];
			float d=im2[j];
			real1[j]=a*c+b*d;
			im1[j]=b*c-a*d;
		}
		real2=null;
		im2=null;
		fft.dorealfft2D(real1,im1,true);
		for(int j=0;j<width*height;j++){
			real1[j]/=(float)(avg1*avg2)*(float)(width*height);
			real1[j]-=1.0f;
			if(brightcorr){
				real1[j]*=(float)Math.sqrt(avg1*avg2);
			}
		}
		manipulate_quads mq=new manipulate_quads();
		if(doavgquad){
			real1=mq.avgquadrants(real1,width,height);
		}
		if(doshiftxcenter&&doshiftycenter){
			real1=mq.shiftxycenter(real1,width,height);
		}else{
			if(doshiftxcenter){
				real1=mq.shiftxcenter(real1,width,height);
			}
			if(doshiftycenter){
				real1=mq.shiftycenter(real1,width,height);
			}
		}
		return real1;
	}

}
