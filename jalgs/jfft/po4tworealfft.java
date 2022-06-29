/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4tworealfft implements Cloneable{
	public po4fft cfft;
	public float[] cosvals;
	public int length;
	public int[] bitorder;
	public int fftindex;

	/*
	 * here is a basic tworeal power of 4 it does the fft of two real data sets
	 * or the inverse fft of two imaginary data sets that result in real outputs
	 * once the po4tworealfft class has been created it can be used again for
	 * other ffts of the same size data must be power of 2 length Copyright Jay
	 * Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public po4tworealfft(int size){
		this(size,0);
	}

	public po4tworealfft(int size,int[] bitorder1){
		this(size,bitorder1,0);
	}

	public po4tworealfft(int size,int[] bitorder1,float[] cosvals1){
		this(size,bitorder1,cosvals1,0);
	}

	public po4tworealfft(int size,int index){
		length=size;
		initcosvals();
		initrevbits();
		fftindex=index;
		cfft=fftutils.construct_fft(length,fftindex,bitorder,cosvals);
	}

	public po4tworealfft(int size,int[] bitorder1,int index){
		length=size;
		initcosvals();
		bitorder=bitorder1;
		fftindex=index;
		cfft=fftutils.construct_fft(length,fftindex,bitorder,cosvals);
	}

	public po4tworealfft(int size,int[] bitorder1,float[] cosvals1,int index){
		length=size;
		cosvals=cosvals1;
		bitorder=bitorder1;
		fftindex=index;
		cfft=fftutils.construct_fft(length,fftindex,bitorder,cosvals);
	}

	public Object clone(){
		return new po4tworealfft(length,bitorder,cosvals,fftindex);
	}

	public void dotworealfft(float[] real1,float[] real2,float[] im1,float[] im2,boolean inverse,boolean bitrev){
		// for the forward transform, interlace the data and then unscramble the
		// results
		// this algorithm is based on two things:
		// First, the fourier transform of f = a + ib is F = A + iB
		// Second, for purely real f, F(n) = F(N-n)*, and for purely imaginary
		// f, F(n) = -(F(N-n)*)
		if(!inverse){
			cfft.dopo4fft(real1,real2,false,bitrev);
			for(int i=1;i<length/2;i++){
				float dum1=real1[i];
				float dum2=real1[length-i];
				float dum3=real2[i];
				float dum4=real2[length-i];
				real1[i]=0.5f*(dum1+dum2);
				real1[length-i]=real1[i];
				real2[i]=0.5f*(dum3+dum4);
				real2[length-i]=real2[i];
				im1[i]=0.5f*(dum3-dum4);
				im1[length-i]=-im1[i];
				im2[i]=0.5f*(-dum1+dum2);
				im2[length-i]=-im2[i];
			}
			im1[0]=0.0f;
			im2[0]=0.0f;
			im1[length/2]=0.0f;
			im2[length/2]=0.0f;
		}else{
			// here we "scramble" the data and then do the inverse transform
			for(int i=0;i<length;i++){
				real1[i]-=im2[i];
				real2[i]+=im1[i];
			}
			cfft.dopo4fft(real1,real2,true,bitrev);
		}
	}

	private void initcosvals(){
		cosvals=new float[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

	private void initrevbits(){
		bitorder=fftutils.get_bitrev_order(length,fftindex+4);
	}

}
