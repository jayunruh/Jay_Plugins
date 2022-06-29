/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class crosscorr3D{
	public po4realfft3D fft;
	public int width,height,slices;

	/*
	 * this class does a 3D spatial crosscorrelation of two ImageJ ImageStack
	 * objects the images must be a power of 2 size in width, height, and # of
	 * slices the images must be real once the class is constructed, it can be
	 * used to do any 3D crosscorrelations on images of the same size Copyright
	 * Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public crosscorr3D(int width1,int height1,int slices1){
		width=width1;
		height=height1;
		slices=slices1;
		fft=new po4realfft3D(width,height,slices);
	}

	public Object[] docrosscorr3D(Object[] stack1,Object[] stack2,boolean doshiftxycenter,boolean doshiftzcenter){
		Object[] temp=docrosscorr3D(stack1,stack2,doshiftxycenter);
		if(doshiftzcenter){
			Object[] temp2=new Object[temp.length];
			for(int i=0;i<temp.length;i++){
				int x=i+temp.length/2;
				if(x>=temp.length)
					x-=temp.length;
				temp2[i]=temp[x];
			}
			return temp2;
		}else{
			return temp;
		}
	}

	public Object[] docrosscorr3D(Object[] stack1,Object[] stack2,boolean doshiftxycenter){
		Object[] real1=new Object[slices];
		Object[] im1=new Object[slices];
		Object[] real2=new Object[slices];
		Object[] im2=new Object[slices];
		double avg1=0.0;
		double avg2=0.0;
		if(stack1[0] instanceof float[]){
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				real2[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real1[i])[j]=((float[])stack1[i])[j];
					avg1+=((float[])real1[i])[j];
					((float[])real2[i])[j]=((float[])stack2[i])[j];
					avg2+=((float[])real2[i])[j];
				}
				im1[i]=new float[width*height];
				im2[i]=new float[width*height];
			}
		}else{
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				real2[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real1[i])[j]=((short[])stack1[i])[j]&0xffff;
					avg1+=((float[])real1[i])[j];
					((float[])real2[i])[j]=((short[])stack2[i])[j]&0xffff;
					avg2+=((float[])real2[i])[j];
				}
				im1[i]=new float[width*height];
				im2[i]=new float[width*height];
			}
		}
		avg1/=((double)width*(double)height*slices);
		avg2/=((double)width*(double)height*slices);
		fft.dorealfft3D(real1,im1,false);
		fft.dorealfft3D(real2,im2,false);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				float a=((float[])real1[i])[j];
				float b=((float[])im1[i])[j];
				float c=((float[])real2[i])[j];
				float d=((float[])im2[i])[j];
				((float[])real1[i])[j]=a*c+b*d;
				((float[])im1[i])[j]=b*c-a*d;
			}
		}
		fft.dorealfft3D(real1,im1,true);
		float tempfloat=(float)width*(float)height*slices*(float)(avg1*avg2);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				((float[])real1[i])[j]/=tempfloat;
				((float[])real1[i])[j]-=1.0f;
			}
			if(doshiftxycenter){
				real1[i]=shiftxycenter((float[])real1[i]);
			}
		}
		return real1;
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
