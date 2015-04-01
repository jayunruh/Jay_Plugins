/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.algutils;

public class fftutils{
	// here we have static utility methods for ffts
	public static int[] get_bitrev_order(int length,int divider){
		// this routine returns an integer array in bit reversed order for ffts
		// divider length stretches are preserved in the data
		int lengthd=length/divider;
		int[] tempbitorder=new int[lengthd];
		for(int i=0;i<lengthd;i++){
			tempbitorder[i]=i;
		}
		int j=0;
		int k;
		for(int i=0;i<lengthd;i++){
			if(j>i){
				tempbitorder[j]=i;
				tempbitorder[i]=j;
			}
			k=lengthd/2;
			while(k>=1&&j>=k){
				j-=k;
				k/=2;
			}
			j+=k;
		}
		int[] bitorder=new int[length];
		for(int i=0;i<lengthd;i++){
			for(k=0;k<divider;k++){
				bitorder[k+i*divider]=tempbitorder[i]+k*lengthd;
			}
		}
		return bitorder;
	}

	public static int[] initrevbits(int length){
		// this has been modified for power of 4 ffts
		int[] bitorder=new int[length];
		for(int i=0;i<length;i++){
			bitorder[i]=i;
		}
		int j=0;
		int k;
		for(int i=0;i<length;i++){
			if(j>i){
				bitorder[j]=i;
				bitorder[i]=j;
			}
			k=length/2;
			while(k>=1&&j>=k){
				j-=k;
				k/=2;
			}
			j+=k;
		}
		for(int i=0;i<length;i+=4){
			int temp=bitorder[i+1];
			bitorder[i+1]=bitorder[i+2];
			bitorder[i+2]=temp;
		}
		return bitorder;
	}

	public static int[] initrevbits2(int length){
		// here is the original for power of 2 length data
		int[] bitorder=new int[length];
		for(int i=0;i<length;i++){
			bitorder[i]=i;
		}
		int j=0;
		int k;
		for(int i=0;i<length;i++){
			if(j>i){
				bitorder[j]=i;
				bitorder[i]=j;
			}
			k=length/2;
			while(k>=1&&j>=k){
				j-=k;
				k/=2;
			}
			j+=k;
		}
		return bitorder;
	}

	public static po4fft construct_fft(int length,int fftindex){
		switch(fftindex){
		case 0:
			return new po4fft(length);
		case 1:
			return new po5fft(length);
		case 2:
			return new po6fft(length);
		case 3:
			return new po7fft(length);
		}
		return new poxfft(length,fftindex+4);
	}

	public static po4fft construct_fft(int length,int fftindex,float[] cosvals){
		switch(fftindex){
		case 0:
			return new po4fft(length,cosvals);
		case 1:
			return new po5fft(length,cosvals);
		case 2:
			return new po6fft(length,cosvals);
		case 3:
			return new po7fft(length,cosvals);
		}
		return new poxfft(length,cosvals,fftindex+4);
	}

	public static po4fft construct_fft(int length,int fftindex,int[] bitorder,float[] cosvals){
		switch(fftindex){
		case 0:
			return new po4fft(length,bitorder,cosvals);
		case 1:
			return new po5fft(length,bitorder,cosvals);
		case 2:
			return new po6fft(length,bitorder,cosvals);
		case 3:
			return new po7fft(length,bitorder,cosvals);
		}
		return new poxfft(length,bitorder,cosvals,fftindex+4);
	}

	public static int trim_length(int length,boolean over,int divider){
		// this static method gives the closest po2/divider length without going
		// over unless over is true
		// if length is already po2/divider, it is returned
		if(length<divider){
			if(!over)
				return 0;
			else
				return divider;
		}
		double p2length=Math.log((double)length/(double)divider)/Math.log(2.0);
		int under=(int)p2length;
		if(p2length==under)
			return length;
		if(!over)
			return divider*(int)Math.pow(2.0,under);
		else
			return divider*(int)Math.pow(2.0,under+1);
	}

	public static int[] get_best_index(int length,boolean over){
		int mintl=trim_length(length,over,4);
		int min=Math.abs(mintl-length);
		int minindex=0;
		for(int i=1;i<4;i++){
			int divider=i+4;
			int tl=trim_length(length,over,divider);
			int dist=Math.abs(tl-length);
			if(dist<min){
				min=dist;
				mintl=tl;
				minindex=i;
			}
		}
		return new int[]{minindex,mintl};
	}

	public static int[] get_best_index(int length,boolean over,int maxdivider){
		// I don't recommend max divider greater than 19
		int mintl=trim_length(length,over,4);
		int min=Math.abs(mintl-length);
		int minindex=4;
		for(int i=5;i<=maxdivider;i++){
			int divider=i;
			int tl=trim_length(length,over,divider);
			int dist=Math.abs(tl-length);
			if(dist<min){
				min=dist;
				mintl=tl;
				minindex=i;
			}
		}
		return new int[]{minindex-4,mintl};
	}
	
	public static float[] complex_amp(Object real,Object im){
		//here real and im are real and imaginary arrays of arbitrary type
		//return type is always float
		float[] rconv=algutils.convert_arr_float2(real);
		float[] iconv=algutils.convert_arr_float2(im);
		float[] amp=new float[rconv.length];
		for(int i=0;i<rconv.length;i++){
			amp[i]=(float)Math.sqrt(rconv[i]*rconv[i]+iconv[i]*iconv[i]);
		}
		return amp;
	}
	
	public static Object[] complex_amp(Object[] real,Object[] im){
		//return type is always float
		Object[] amp=new Object[real.length];
		for(int i=0;i<real.length;i++){
			amp[i]=complex_amp(real[i],im[i]);
		}
		return amp;
	}
	
	public static float[] complex_phase(Object real,Object im){
		//here real and im are real and imaginary arrays of arbitrary type
		//return type is always float
		float[] rconv=algutils.convert_arr_float2(real);
		float[] iconv=algutils.convert_arr_float2(im);
		float[] phase=new float[rconv.length];
		for(int i=0;i<rconv.length;i++){
			phase[i]=(float)Math.atan2(iconv[i],rconv[i]);
		}
		return phase;
	}
	
	public static Object[] complex_phase(Object[] real,Object[] im){
		//return type is always float
		Object[] phase=new Object[real.length];
		for(int i=0;i<real.length;i++){
			phase[i]=complex_phase(real[i],im[i]);
		}
		return phase;
	}

}
