/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class autocorr3D{
	public po4realfft3D fft; // this is the 3D fft class used for the
	// autocorrelation
	public int width,height,slices;

	/*
	 * this class does a 3D spatial autocorrelation of an ImageJ ImageStack
	 * object the image must be a power of 2 size in width, height, and # of
	 * slices the image must be real once the class is constructed, it can be
	 * used to do any 3D autocorrelations on images of the same size Copyright
	 * Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public autocorr3D(int width1,int height1,int slices1){
		width=width1;
		height=height1;
		slices=slices1;
		fft=new po4realfft3D(width,height,slices);
	}

	public Object[] doautocorr3D(Object[] stack,boolean doshiftxycenter){
		// do shiftxycenter shifts the data so that the zero point is in the xy
		// center of the image
		// the returned Object array is an array of slices
		Object[] real=new Object[slices];
		Object[] im=new Object[slices];
		double avg=0.0;
		if(stack[0] instanceof float[]){
			for(int i=0;i<slices;i++){
				real[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real[i])[j]=((float[])stack[i])[j];
					avg+=((float[])real[i])[j];
				}
				im[i]=new float[width*height];
			}
		}else{
			for(int i=0;i<slices;i++){
				real[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real[i])[j]=((short[])stack[i])[j]&0xffff;
					avg+=((float[])real[i])[j];
				}
				im[i]=new float[width*height];
			}
		}
		avg/=((double)width*(double)height*slices);
		fft.dorealfft3D(real,im,false);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				float a=((float[])real[i])[j];
				float b=((float[])im[i])[j];
				((float[])real[i])[j]=a*a+b*b;
			}
			im[i]=new float[width*height];
		}
		fft.dorealfft3D(real,im,true);
		float tempfloat=(float)width*(float)height*slices*(float)(avg*avg);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				((float[])real[i])[j]/=tempfloat;
				((float[])real[i])[j]-=1.0f;
			}
			if(doshiftxycenter){
				real[i]=shiftxycenter((float[])real[i]);
			}
		}
		return real;
	}

	private float[] shiftxycenter(float[] data){
		float[] temp=new float[width*height];
		for(int i=0;i<height;i++){
			int dumy=i+height/2;
			if(dumy>=height){
				dumy-=height;
			}
			for(int j=0;j<width;j++){
				int dumx=j+width/2;
				if(dumx>=width){
					dumx-=width;
				}
				temp[j+i*width]=data[dumx+width*dumy];
			}
		}
		return temp;
	}

}
