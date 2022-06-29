/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4cfft2D implements Cloneable{
	public po4fft cfft1,cfft2; // these are the 1D fft classes used for the 2D
	// fft
	public float[] cosvals1,cosvals2;
	public int length1,length2;
	public int[] bitorder1,bitorder2;

	/*
	 * This class does complex 2D ffts once the po4cfft2D class has been created
	 * it can be used again for other 2D ffts of the same size width and height
	 * must be power of 2 size Copyright Jay Unruh Stowers Institute for Medical
	 * Research 4/25/08
	 */

	public po4cfft2D(int length11,int length12){
		this(length11,length12,0,0);
	}

	public po4cfft2D(int length11,int length12,int fftindex1,int fftindex2){
		length1=length11;
		length2=length12;
		initcosvals();
		initrevbits(fftindex1,fftindex2);
		cfft1=fftutils.construct_fft(length1,fftindex1,bitorder1,cosvals1);
		cfft2=fftutils.construct_fft(length2,fftindex2,bitorder2,cosvals2);
	}

	public po4cfft2D(int length11,int length12,float[] cosvals11,float[] cosvals21){
		this(length11,length12,0,0,cosvals11,cosvals21);
	}

	public po4cfft2D(int length11,int length12,int fftindex1,int fftindex2,float[] cosvals11,float[] cosvals21){
		length1=length11;
		length2=length12;
		cosvals1=cosvals11;
		cosvals2=cosvals21;
		initrevbits(fftindex1,fftindex2);
		cfft1=fftutils.construct_fft(length1,fftindex1,bitorder1,cosvals1);
		cfft2=fftutils.construct_fft(length2,fftindex2,bitorder2,cosvals2);
	}

	public po4cfft2D(int length11,int length12,int fftindex1,int fftindex2,float[] cosvals11,float[] cosvals21,int[] bitorder11,int[] bitorder21){
		length1=length11;
		length2=length12;
		cosvals1=cosvals11;
		cosvals2=cosvals21;
		bitorder1=bitorder11;
		bitorder2=bitorder21;
		cfft1=fftutils.construct_fft(length1,fftindex1,bitorder1,cosvals1);
		cfft2=fftutils.construct_fft(length2,fftindex2,bitorder2,cosvals2);
	}

	public Object clone(){
		return new po4cfft2D(length1,length2,cfft1.fftindex,cfft2.fftindex,cosvals1.clone(),cosvals2.clone());
	}

	public void docfft2D(float[] real,float[] im,boolean inverse){
		// start by doing the horizontal fft
		// copy the data into arrays in bit reversed order
		for(int i=0;i<length2;i++){
			float[] tempfloat=new float[length1];
			float[] tempfloat2=new float[length1];
			for(int j=0;j<length1;j++){
				tempfloat[j]=real[i*length1+bitorder1[j]];
				tempfloat2[j]=im[i*length1+bitorder1[j]];
			}
			cfft1.dopo4fft(tempfloat,tempfloat2,inverse,false);
			for(int j=0;j<length1;j++){
				real[i*length1+j]=tempfloat[j];
				im[i*length1+j]=tempfloat2[j];
			}
		}
		// now the vertical fft
		for(int i=0;i<length1;i++){
			float[] tempfloat=new float[length2];
			float[] tempfloat2=new float[length2];
			for(int j=0;j<length2;j++){
				tempfloat[j]=real[bitorder2[j]*length1+i];
				tempfloat2[j]=im[bitorder2[j]*length1+i];
			}
			cfft2.dopo4fft(tempfloat,tempfloat2,inverse,false);
			for(int j=0;j<length2;j++){
				real[j*length1+i]=tempfloat[j];
				im[j*length1+i]=tempfloat2[j];
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
