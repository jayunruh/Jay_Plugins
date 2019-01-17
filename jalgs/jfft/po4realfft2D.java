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

public class po4realfft2D implements Cloneable{
	public po4tworealfft r2fft; // here is the tworealfft used for this
	// calculation
	public po4fft cfft;
	public float[] cosvals1,cosvals2;
	public int length1,length2;
	public int[] bitorder1,bitorder2;

	/*
	 * here is a basic real power of 4 2Dfft once the po4realfft2D class has
	 * been created it can be used again for other ffts of the same size data
	 * must be power of 2 size in width and height Copyright Jay Unruh Stowers
	 * Institute for Medical Research 4/25/08
	 */

	public po4realfft2D(int length11,int length12){
		length1=length11;
		length2=length12;
		initcosvals();
		initrevbits(0,0);
		r2fft=new po4tworealfft(length2,bitorder2,cosvals2);
		cfft=new po4fft(length1,bitorder1,cosvals1);
	}

	public po4realfft2D(int length11,int length12,int fftindex1,int fftindex2){
		length1=length11;
		length2=length12;
		initcosvals();
		initrevbits(fftindex1,fftindex2);
		r2fft=new po4tworealfft(length2,bitorder2,cosvals2,fftindex2);
		// cfft=new po4fft(length1,bitorder1,cosvals1,fftindex1);
		cfft=fftutils.construct_fft(length1,fftindex1,bitorder1,cosvals1);
	}

	public po4realfft2D(int length11,int length12,float[] cosvals11,float[] cosvals21){
		length1=length11;
		length2=length12;
		cosvals1=cosvals11;
		cosvals2=cosvals21;
		initrevbits(4,4);
		r2fft=new po4tworealfft(length2,bitorder2,cosvals2);
		cfft=new po4fft(length1,bitorder1,cosvals1);
	}


	
	public po4realfft2D(int length11,int length12,float[] cosvals11,float[] cosvals21,int fftindex1,int fftindex2){
		length1=length11;
		length2=length12;
		cosvals1=cosvals11;
		cosvals2=cosvals21;
		initrevbits(fftindex1,fftindex2);
		r2fft=new po4tworealfft(length2,bitorder2,cosvals2,fftindex2);
		// cfft=new po4fft(length1,bitorder1,cosvals1);
		cfft=fftutils.construct_fft(length1,fftindex1,bitorder1,cosvals1);
	}
	
	public po4realfft2D(int length11,int length12,float[] cosvals11,float[] cosvals21,po4tworealfft r2fft,po4fft cfft){
		length1=length11;
		length2=length12;
		cosvals1=cosvals11;
		cosvals2=cosvals21;
		//initrevbits(fftindex1,fftindex2);
		bitorder1=cfft.bitorder;
		bitorder2=r2fft.bitorder;
		this.r2fft=r2fft;
		// cfft=new po4fft(length1,bitorder1,cosvals1);
		this.cfft=cfft;
	}
	
	public Object clone(){
		return new po4realfft2D(length1,length2,cosvals1.clone(),cosvals2.clone(),(po4tworealfft)r2fft.clone(),(po4fft)cfft.clone());
	}

	public void dorealfft2D(float[] real,float[] im,boolean inverse){
		if(!inverse){
			// start by doing the vertical fft
			// copy the data into arrays in bit reversed order
			for(int i=0;i<length1;i+=2){
				float[] tempfloat=new float[length2];
				float[] tempfloat2=new float[length2];
				float[] tempfloat3=new float[length2];
				float[] tempfloat4=new float[length2];
				for(int j=0;j<length2;j++){
					tempfloat[j]=real[bitorder2[j]*length1+i];
					tempfloat2[j]=real[bitorder2[j]*length1+i+1];
				}
				r2fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,false,false);
				for(int j=0;j<length2;j++){
					real[j*length1+i]=tempfloat[j];
					real[j*length1+i+1]=tempfloat2[j];
					im[j*length1+i]=tempfloat3[j];
					im[j*length1+i+1]=tempfloat4[j];
				}
			}
			// now do the horizontal ffts (these are complex)
			for(int i=0;i<length2;i++){
				float[] tempfloat=new float[length1];
				float[] tempfloat2=new float[length1];
				for(int j=0;j<length1;j++){
					tempfloat[j]=real[i*length1+bitorder1[j]];
					tempfloat2[j]=im[i*length1+bitorder1[j]];
				}
				cfft.dopo4fft(tempfloat,tempfloat2,false,false);
				for(int j=0;j<length1;j++){
					real[i*length1+j]=tempfloat[j];
					im[i*length1+j]=tempfloat2[j];
				}
			}
		}else{
			// here start with the horizontal complex ffts
			for(int i=0;i<length2;i++){
				float[] tempfloat=new float[length1];
				float[] tempfloat2=new float[length1];
				for(int j=0;j<length1;j++){
					tempfloat[j]=real[i*length1+bitorder1[j]];
					tempfloat2[j]=im[i*length1+bitorder1[j]];
				}
				cfft.dopo4fft(tempfloat,tempfloat2,true,false);
				for(int j=0;j<length1;j++){
					real[i*length1+j]=tempfloat[j];
					im[i*length1+j]=tempfloat2[j];
				}
			}
			// now the vertical tworeal ffts
			for(int i=0;i<length1;i+=2){
				float[] tempfloat=new float[length2];
				float[] tempfloat2=new float[length2];
				float[] tempfloat3=new float[length2];
				float[] tempfloat4=new float[length2];
				for(int j=0;j<length2;j++){
					tempfloat[j]=real[bitorder2[j]*length1+i];
					tempfloat2[j]=real[bitorder2[j]*length1+i+1];
					tempfloat3[j]=im[bitorder2[j]*length1+i];
					tempfloat4[j]=im[bitorder2[j]*length1+i+1];
				}
				r2fft.dotworealfft(tempfloat,tempfloat2,tempfloat3,tempfloat4,true,false);
				for(int j=0;j<length2;j++){
					real[j*length1+i]=tempfloat[j];
					real[j*length1+i+1]=tempfloat2[j];
				}
			}
		}
	}
	


	private void initcosvals(){
		cosvals1=new float[length1];
		for(int i=0;i<length1;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length1);
			cosvals1[i]=(float)Math.cos(dumdbl);
		}
		if(length2!=length1){
			cosvals2=new float[length2];
			for(int i=0;i<length2;i++){
				double dumdbl=Math.PI*2.0*((double)i/(double)length2);
				cosvals2[i]=(float)Math.cos(dumdbl);
			}
		}else{
			cosvals2=cosvals1;
		}
	}

	private void initrevbits(int fftindex1,int fftindex2){
		bitorder1=fftutils.get_bitrev_order(length1,fftindex1+4);
		if(length2==length1)
			bitorder2=bitorder1;
		else
			bitorder2=fftutils.get_bitrev_order(length2,fftindex2+4);
	}

}


