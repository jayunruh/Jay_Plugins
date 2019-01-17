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
	
	public crosscorr2D(int width1,int height1,int fftindex1,int fftindex2){
		width=width1;
		height=height1;
		fft=new po4realfft2D(width,height,fftindex1,fftindex2);
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
			real1[j]/=(float)(avg1*avg2)*(width*height);
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
	
	public float[][] phaseCorr(float[] pix1,float[] pix2,boolean outcrosscorr){
		float[] real1=pix1.clone();
		float[] real2=pix2.clone();
		float[] im1=new float[width*height];
		float[] im2=new float[width*height];
		fft.dorealfft2D(real1,im1,false);
		fft.dorealfft2D(real2,im2,false);
		
		//now normalize the fft pixels and conjugate multiply
		double avg1=0.0, avg2=0.0, stdev1=0.0, stdev2=0.0;
		for(int i=0;i<width*height;i++){
			avg1+=pix1[i];
			avg2+=pix2[i];
			stdev1+=pix1[i]*pix1[i];
			stdev2+=pix2[i]*pix2[i];
			float length1=(float)Math.sqrt(real1[i]*real1[i]+im1[i]*im1[i]);
			float a1=real1[i], b1=im1[i], c1=real2[i], d1=im2[i]; //these are the raw values
			if(length1<(float)0.00001){
				real1[i]=0.0f; im1[i]=0.0f;
			} else {
				real1[i]/=length1; im1[i]/=length1;
			}
			float length2=(float)Math.sqrt(real2[i]*real2[i]+im2[i]*im2[i]);
			if(length2<(float)0.00001){
				real2[i]=0.0f; im2[i]=0.0f;
			} else {
				real2[i]/=length2; im2[i]/=length2;
			}
			float a=real1[i], b=im1[i], c=real2[i], d=im2[i]; //these are the normalized (and filtered) values
			real1[i]=a*c+b*d;
			im1[i]=a*d-b*c;
			real2[i]=a1*c1+b1*d1;
			im2[i]=a1*d1-b1*c1;
		}
		fft.dorealfft2D(real1,im1,true); //this is the phase correlation
		if(!outcrosscorr) return new float[][]{real1};
		else{ 
			fft.dorealfft2D(real2,im2,true); //this is the cross correlation (unnormalized)
    		double totsize=width*height;
    		avg1/=totsize;
    		avg2/=totsize;
    		stdev1/=totsize;
    		stdev2/=totsize;
    		stdev1-=avg1*avg1; stdev2-=avg2*avg2;
    		stdev1=Math.sqrt(stdev1); stdev2=Math.sqrt(stdev2);
    		float tempa=(float)(avg1*avg2*totsize);
    		float tempb=(float)(stdev1*stdev2*totsize);
    		for(int i=0;i<width*height;i++) real2[i]=(real2[i]-tempa)/tempb;
    		im2=null; im1=null;
    		return new float[][]{real1,real2};
		}
	}
	
	public float[][] phaseCorr(float[] real,float[] im,float[] real2,float[] im2){
		//here we start with already fft'd data and don't output the crosscorr
		//now normalize the fft pixels and conjugate multiply
		float[] real1=real.clone();
		float[] im1=im.clone();
		for(int i=0;i<width*height;i++){
			float length1=(float)Math.sqrt(real1[i]*real1[i]+im1[i]*im1[i]);
			float a=0.0f,b=0.0f;
			if(length1>=(float)0.00001){
				a=real1[i]/length1; b=im1[i]/length1;
			}
			float length2=(float)Math.sqrt(real2[i]*real2[i]+im2[i]*im2[i]);
			float c=0.0f,d=0.0f;
			if(length2>=(float)0.00001){
				c=real2[i]/length2; d=im2[i]/length2;
			}
			real1[i]=a*c+b*d;
			im1[i]=a*d-b*c;
		}
		//fft.dorealfft2D(real1,im1,true); //this is the phase correlation--the calling program can do this
		return new float[][]{real1,im1};
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
			real1[j]/=(float)(avg1*avg2)*(width*height);
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
