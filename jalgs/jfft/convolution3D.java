/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class convolution3D{
	public po4realfft3D fft;
	public int width,height,slices;

	/*
	 * this class does a 3D spatial crosscorrelation of two ImageJ ImageStack
	 * objects the images must be a power of 2 size in width, height, and # of
	 * slices the images must be real once the class is constructed, it can be
	 * used to do any 3D crosscorrelations on images of the same size Copyright
	 * Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public convolution3D(int width1,int height1,int slices1){
		this(width1,height1,slices1,0,0,0);
	}

	public convolution3D(int width1,int height1,int slices1,int fftindex1,int fftindex2,int fftindex3){
		width=width1;
		height=height1;
		slices=slices1;
		fft=new po4realfft3D(width,height,slices,fftindex1,fftindex2,fftindex3);
	}

	public Object[] convolve3D(Object[] stack1,Object[] stack2){
		Object[] real1=new Object[slices];
		Object[] im1=new Object[slices];
		Object[] real2=new Object[slices];
		Object[] im2=new Object[slices];
		if(stack1[0] instanceof float[]){
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				real2[i]=new float[width*height];
				System.arraycopy(stack1[i],0,real1[i],0,width*height);
				System.arraycopy(stack2[i],0,real2[i],0,width*height);
				im1[i]=new float[width*height];
				im2[i]=new float[width*height];
			}
		}else{
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				real2[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real1[i])[j]=((short[])stack1[i])[j]&0xffff;
					((float[])real2[i])[j]=((short[])stack2[i])[j]&0xffff;
				}
				im1[i]=new float[width*height];
				im2[i]=new float[width*height];
			}
		}
		fft.dorealfft3D(real1,im1,false);
		fft.dorealfft3D(real2,im2,false);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				float a=((float[])real1[i])[j];
				float b=((float[])im1[i])[j];
				float c=((float[])real2[i])[j];
				float d=((float[])im2[i])[j];
				((float[])real1[i])[j]=a*c-b*d;
				((float[])im1[i])[j]=b*c+a*d;
			}
		}
		fft.dorealfft3D(real1,im1,true);
		return real1;
	}

	public Object[] convolve3D(Object[] stack1,Object[] rstack2,Object[] istack2){
		Object[] real1=new Object[slices];
		Object[] im1=new Object[slices];
		if(stack1[0] instanceof float[]){
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				System.arraycopy(stack1[i],0,real1[i],0,width*height);
				im1[i]=new float[width*height];
			}
		}else{
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real1[i])[j]=((short[])stack1[i])[j]&0xffff;
				}
				im1[i]=new float[width*height];
			}
		}
		fft.dorealfft3D(real1,im1,false);
		for(int i=0;i<slices;i++){
			for(int j=0;j<width*height;j++){
				float a=((float[])real1[i])[j];
				float b=((float[])im1[i])[j];
				float c=((float[])rstack2[i])[j];
				float d=((float[])istack2[i])[j];
				((float[])real1[i])[j]=a*c-b*d;
				((float[])im1[i])[j]=b*c+a*d;
			}
		}
		fft.dorealfft3D(real1,im1,true);
		return real1;
	}

	public Object[] convolve3D_refl(Object[] stack1,Object[] rstack2,Object[] istack2){
		// here we convolve with a reflected rstack2 and istack2
		// for EM deconvolution
		Object[] real1=new Object[slices];
		Object[] im1=new Object[slices];
		if(stack1[0] instanceof float[]){
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				System.arraycopy(stack1[i],0,real1[i],0,width*height);
				im1[i]=new float[width*height];
			}
		}else{
			for(int i=0;i<slices;i++){
				real1[i]=new float[width*height];
				for(int j=0;j<width*height;j++){
					((float[])real1[i])[j]=((short[])stack1[i])[j]&0xffff;
				}
				im1[i]=new float[width*height];
			}
		}
		fft.dorealfft3D(real1,im1,false);
		for(int i=0;i<slices;i++){
			int refi=slices-i;
			if(i==0)
				refi=0;
			for(int j=0;j<height;j++){
				int refj=height-j;
				if(j==0)
					refj=0;
				refj*=width;
				int offset=j*width;
				for(int k=0;k<width;k++){
					int refk=width-k;
					if(k==0)
						refk=0;
					float a=((float[])real1[i])[k+offset];
					float b=((float[])im1[i])[k+offset];
					float c=((float[])rstack2[refi])[refk+refj];
					float d=((float[])istack2[refi])[refk+refj];
					((float[])real1[i])[k+offset]=a*c-b*d;
					((float[])im1[i])[k+offset]=b*c+a*d;
				}
			}
		}
		fft.dorealfft3D(real1,im1,true);
		return real1;
	}

	public Object[] convolve3D(Object[] rstack1,Object[] istack1,Object[] rstack2,Object[] istack2){
		Object[] tempr=new Object[slices];
		Object[] tempi=new Object[slices];
		for(int i=0;i<slices;i++){
			tempr[i]=new float[width*height];
			tempi[i]=new float[width*height];
			for(int j=0;j<width*height;j++){
				float a=((float[])rstack1[i])[j];
				float b=((float[])istack1[i])[j];
				float c=((float[])rstack2[i])[j];
				float d=((float[])istack2[i])[j];
				((float[])tempr[i])[j]=a*c-b*d;
				((float[])tempi[i])[j]=b*c+a*d;
			}
		}
		fft.dorealfft3D(tempr,tempi,true);
		return tempr;
	}

	public Object[][] convolve3Dinverse(Object[] rstack1,Object[] istack1,Object[] rstack2,Object[] istack2){
		// here we basically just multiply to two complex distributions
		Object[][] temp=new Object[2][slices];
		for(int i=0;i<slices;i++){
			temp[0][i]=new float[width*height];
			temp[1][i]=new float[width*height];
			for(int j=0;j<width*height;j++){
				float a=((float[])rstack1[i])[j];
				float b=((float[])istack1[i])[j];
				float c=((float[])rstack2[i])[j];
				float d=((float[])istack2[i])[j];
				((float[])temp[0][i])[j]=a*c-b*d;
				((float[])temp[1][i])[j]=b*c+a*d;
			}
		}
		return temp;
	}

}
