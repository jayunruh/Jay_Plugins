/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class po4realfft3D{
	public po4tworealfft r2fft; // here is the tworealfft used for the z
	// dimension
	public po4cfft2D cfft2D; // here is the complex 2D fft used for the xy
	// dimesions
	public float[] cosvals;
	public int length1,length2,length3;
	public int[] bitorder;
	public int nthreads=1;

	/*
	 * here is a basic real power of 4 3Dfft once the po4realfft3D class has
	 * been created it can be used again for other ffts of the same size data
	 * must be power of 2 size in width, height, and number of slices Copyright
	 * Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public po4realfft3D(int length11,int length12,int length13){
		this(length11,length12,length13,0,0,0);
	}

	public po4realfft3D(int length11,int length12,int length13,int fftindex1,int fftindex2,int fftindex3){
		length1=length11;
		length2=length12;
		length3=length13;
		cfft2D=new po4cfft2D(length1,length2,fftindex1,fftindex2);
		// only create new cosvals and bitorder arrays if needed
		if(length3==length1){
			cosvals=cfft2D.cosvals1;
			bitorder=cfft2D.bitorder1;
		}else{
			if(length3==length2){
				cosvals=cfft2D.cosvals2;
				bitorder=cfft2D.bitorder2;
			}else{
				initcosvals();
				initrevbits(fftindex3);
			}
		}
		r2fft=new po4tworealfft(length3,bitorder,cosvals,fftindex3);
	}

	public void dorealfft3D(Object[] real,Object[] im,boolean inverse){
		ExecutorService executor=null;
		if(nthreads>1)
			executor=Executors.newFixedThreadPool(nthreads);
		if(!inverse){
			// start by doing the real z ffts
			// copy the data into arrays in bit reversed order
			for(int i=0;i<length1*length2;i+=2){
				if(nthreads>1){
					Runnable worker=new fft1D(real,im,i,false,r2fft);
					executor.execute(worker);
				}else{
					float[] tempfloat=new float[length3];
					float[] tempfloat2=new float[length3];
					float[] tempfloat3=new float[length3];
					float[] tempfloat4=new float[length3];
					for(int j=0;j<length3;j++){
						tempfloat[j]=((float[])real[bitorder[j]])[i];
						tempfloat2[j]=((float[])real[bitorder[j]])[i+1];
					}
					r2fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,false,false);
					for(int j=0;j<length3;j++){
						((float[])real[j])[i]=tempfloat[j];
						((float[])real[j])[i+1]=tempfloat2[j];
						((float[])im[j])[i]=tempfloat3[j];
						((float[])im[j])[i+1]=tempfloat4[j];
					}
				}
			}
			if(nthreads>1){
				executor.shutdown();
				// Wait until all threads are finished
				while(!executor.isTerminated()){
				}
			}
			if(nthreads>1)
				executor=Executors.newFixedThreadPool(nthreads);
			// now do the xy ffts (these are complex)
			for(int i=0;i<length3;i++){
				if(nthreads>1){
					Runnable worker=new fft2D((float[])real[i],(float[])im[i],false,cfft2D);
					executor.execute(worker);
				}else{
					cfft2D.docfft2D((float[])real[i],(float[])im[i],false);
				}
			}
			if(nthreads>1){
				executor.shutdown();
				// Wait until all threads are finished
				while(!executor.isTerminated()){
				}
			}
		}else{
			// here start with the xy complex ffts
			for(int i=0;i<length3;i++){
				if(nthreads>1){
					Runnable worker=new fft2D((float[])real[i],(float[])im[i],true,cfft2D);
					executor.execute(worker);
				}else{
					cfft2D.docfft2D((float[])real[i],(float[])im[i],true);
				}
			}
			if(nthreads>1){
				executor.shutdown();
				// Wait until all threads are finished
				while(!executor.isTerminated()){
				}
			}
			if(nthreads>1)
				executor=Executors.newFixedThreadPool(nthreads);
			// now the z tworeal ffts
			for(int i=0;i<length1*length2;i+=2){
				if(nthreads>1){
					Runnable worker=new fft1D(real,im,i,true,r2fft);
					executor.execute(worker);
				}else{
					float[] tempfloat=new float[length3];
					float[] tempfloat2=new float[length3];
					float[] tempfloat3=new float[length3];
					float[] tempfloat4=new float[length3];
					for(int j=0;j<length3;j++){
						tempfloat[j]=((float[])real[bitorder[j]])[i];
						tempfloat2[j]=((float[])real[bitorder[j]])[i+1];
						tempfloat3[j]=((float[])im[bitorder[j]])[i];
						tempfloat4[j]=((float[])im[bitorder[j]])[i+1];
					}
					r2fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,true,false);
					for(int j=0;j<length3;j++){
						((float[])real[j])[i]=tempfloat[j];
						((float[])real[j])[i+1]=tempfloat2[j];
					}
				}
			}
			if(nthreads>1){
				executor.shutdown();
				// Wait until all threads are finished
				while(!executor.isTerminated()){
				}
			}
		}
	}

	private void initcosvals(){
		cosvals=new float[length3];
		for(int i=0;i<length3;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length3);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

	private void initrevbits(int fftindex){
		bitorder=fftutils.get_bitrev_order(length3,fftindex+4);
	}

}

class fft2D implements Runnable{
	float[] real,im;
	po4cfft2D cfft;
	boolean inverse;

	public fft2D(float[] real,float[] im,boolean inverse,po4cfft2D cfft){
		this.cfft=(po4cfft2D)cfft.clone();
		this.real=real;
		this.im=im;
		this.inverse=inverse;
	}

	public void run(){
		cfft.docfft2D(real,im,inverse);
	}

}

class fft1D implements Runnable{
	Object[] real,im;
	int index;
	po4tworealfft fft;
	boolean inverse;

	public fft1D(Object[] real,Object[] im,int index,boolean inverse,po4tworealfft fft){
		this.fft=(po4tworealfft)fft.clone();
		this.real=real;
		this.im=im;
		this.inverse=inverse;
		this.index=index;
	}

	public void run(){
		// copy the data out in bitrev order
		int[] bitorder=fft.bitorder;
		int length=real.length;
		if(!inverse){
			float[] tempfloat=new float[length];
			float[] tempfloat2=new float[length];
			float[] tempfloat3=new float[length];
			float[] tempfloat4=new float[length];
			for(int j=0;j<length;j++){
				tempfloat[j]=((float[])real[bitorder[j]])[index];
				tempfloat2[j]=((float[])real[bitorder[j]])[index+1];
			}
			fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,inverse,false);
			for(int j=0;j<length;j++){
				((float[])real[j])[index]=tempfloat[j];
				((float[])real[j])[index+1]=tempfloat2[j];
				((float[])im[j])[index]=tempfloat3[j];
				((float[])im[j])[index+1]=tempfloat4[j];
			}
		}else{
			float[] tempfloat=new float[length];
			float[] tempfloat2=new float[length];
			float[] tempfloat3=new float[length];
			float[] tempfloat4=new float[length];
			for(int j=0;j<length;j++){
				tempfloat[j]=((float[])real[bitorder[j]])[index];
				tempfloat2[j]=((float[])real[bitorder[j]])[index+1];
				tempfloat3[j]=((float[])im[bitorder[j]])[index];
				tempfloat4[j]=((float[])im[bitorder[j]])[index+1];
			}
			fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,inverse,false);
			for(int j=0;j<length;j++){
				((float[])real[j])[index]=tempfloat[j];
				((float[])real[j])[index+1]=tempfloat2[j];
			}
		}
	}

}
