/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class poxfft extends po4fft implements Cloneable{
	// public int[] bitorder;
	// public float[] cosvals;
	// public int length;
	// this needs more work
	public int fftindex=4;
	public int kernelsize;
	public float[] c; // these are the coefficients of the x sized fft kernel
	public float[] s; //

	/*
	 * here is a basic complex power of 2/x complex fft (x,2x,4x,8x,...) once
	 * the poxfft class has been created it can be used again for other ffts of
	 * the same size data must be correct length note that bit reversal is made
	 * optional to allow for the calling class to do the bit reversal ahead of
	 * time. This can save extra copying If the bitreversal is done ahead of
	 * time, it is wise to initialize this class with the current bitorder array
	 * to save memory this same logic applies to all of the fft classes
	 */

	public poxfft(int size,int kernelsize){
		super();
		length=size;
		this.kernelsize=kernelsize;
		initkernel();
		initrevbits();
		initcosvals();
	}

	public poxfft(int size,int[] bitorder1,float[] cosvals1,int kernelsize){
		length=size;
		bitorder=bitorder1;
		cosvals=cosvals1;
		this.kernelsize=kernelsize;
		initkernel();
	}

	public poxfft(int size,float[] cosvals1,int kernelsize){
		length=size;
		cosvals=cosvals1;
		this.kernelsize=kernelsize;
		initkernel();
		initrevbits();
	}
	
	public Object clone(){
		return new poxfft(length,bitorder.clone(),cosvals.clone(),kernelsize);
	}

	public void dopo4fft(float[] data,boolean inverse){
		float[] r=new float[length];
		float[] im=new float[length];
		for(int i=0;i<length;i++){
			r[i]=data[bitorder[i]+bitorder[i]];
			im[i]=data[bitorder[i]+bitorder[i]+1];
		}
		dopo4fft(r,im,inverse,false);
		for(int i=0;i<length;i++){
			data[i+i]=r[i];
			data[i+i+1]=im[i];
		}
	}

	public void dopo4fft(double[] data,boolean inverse){
		float[] r=new float[length];
		float[] im=new float[length];
		for(int i=0;i<length;i++){
			r[i]=(float)data[bitorder[i]+bitorder[i]];
			im[i]=(float)data[bitorder[i]+bitorder[i]+1];
		}
		dopo4fft(r,im,inverse,false);
		for(int i=0;i<length;i++){
			data[i+i]=r[i];
			data[i+i+1]=im[i];
		}
	}

	public void dopo4fft(float[] r,float[] im,boolean inverse,boolean bitrev){
		// split the data into almost bit reversed real and imaginary arrays
		// (except for kernels)
		int ld2=length/2;
		int ld4=length/4;
		float[] rt,imt;
		if(bitrev){
			rt=new float[length];
			imt=new float[length];
			for(int i=0;i<length;i++){
				rt[i]=r[i];
				imt[i]=im[i];
			}
			for(int i=0;i<length;i++){
				r[i]=rt[bitorder[i]];
				im[i]=imt[bitorder[i]];
			}
		}
		rt=new float[kernelsize];
		imt=new float[kernelsize];
		// do the initial five point transforms
		int k2size=(int)(0.5*kernelsize);
		if(!inverse){
			for(int i=0;i<length;i+=kernelsize){
				for(int j=0;j<kernelsize;j++){
					rt[j]=r[i+j];
					r[i+j]=0.0f;
					imt[j]=im[i+j];
					im[i+j]=0.0f;
				}
				for(int j=0;j<kernelsize;j++){ // here we iterate over
					// the frequencies
					int k=j;
					if(k>k2size)
						k-=kernelsize;
					int sp=Math.abs(k);
					if(k>=0){
						for(int t=0;t<kernelsize;t++){ // here we
							// integrate
							// over the time
							// points
							r[i+j]+=c[(t*sp)%kernelsize]*rt[t]+s[(t*sp)%kernelsize]*imt[t];
							im[i+j]+=c[(t*sp)%kernelsize]*imt[t]-s[(t*sp)%kernelsize]*rt[t];
						}
					}else{
						for(int t=0;t<kernelsize;t++){ // here we
							// integrate
							// over the time
							// points
							r[i+j]+=c[(t*sp)%kernelsize]*rt[t]-s[(t*sp)%kernelsize]*imt[t];
							im[i+j]+=s[(t*sp)%kernelsize]*rt[t]+c[(t*sp)%kernelsize]*imt[t];
						}
					}
				}
			}
		}else{
			for(int i=0;i<length;i+=kernelsize){
				for(int j=0;j<kernelsize;j++){
					rt[j]=r[i+j];
					r[i+j]=0.0f;
					imt[j]=im[i+j];
					im[i+j]=0.0f;
				}
				for(int j=0;j<kernelsize;j++){ // here we iterate over
					// the frequencies
					int k=j;
					if(k>k2size)
						k-=kernelsize;
					int sp=Math.abs(k);
					if(k>=0){
						for(int t=0;t<kernelsize;t++){ // here we
							// integrate
							// over the time
							// points
							r[i+j]+=c[(t*sp)%kernelsize]*rt[t]-s[(t*sp)%kernelsize]*imt[t];
							im[i+j]+=s[(t*sp)%kernelsize]*rt[t]+c[(t*sp)%kernelsize]*imt[t];
						}
					}else{
						for(int t=0;t<kernelsize;t++){ // here we
							// integrate
							// over the time
							// points
							r[i+j]+=c[(t*sp)%kernelsize]*rt[t]+s[(t*sp)%kernelsize]*imt[t];
							im[i+j]+=c[(t*sp)%kernelsize]*imt[t]-s[(t*sp)%kernelsize]*rt[t];
						}
					}
				}
			}
		}
		rt=null;
		imt=null;
		// now combine these into ever increasing size transforms to get the
		// final fft
		int increment=kernelsize;
		if(!inverse){
			while(increment<length){
				for(int i=0;i<length;i+=(increment+increment)){
					for(int k=0;k<increment;k++){
						int trigindex=k*(ld2/increment);
						float tempcos=cosvals[trigindex];
						float tempsin;
						if(trigindex<ld4){
							tempsin=cosvals[trigindex+ld4+ld2];
						}else{
							tempsin=cosvals[trigindex-ld4];
						}
						int temp1=i+k;
						int temp2=temp1+increment;
						float flt1=tempcos*r[temp2]+tempsin*im[temp2];
						float flt2=tempcos*im[temp2]-tempsin*r[temp2];
						r[temp2]=r[temp1]-flt1;
						im[temp2]=im[temp1]-flt2;
						r[temp1]+=flt1;
						im[temp1]+=flt2;
					}
				}
				increment+=increment;
			}
		}else{
			while(increment<length){
				for(int i=0;i<length;i+=(increment+increment)){
					for(int k=0;k<increment;k++){
						int trigindex=k*(ld2/increment);
						float tempcos=cosvals[trigindex];
						float tempsin;
						if(trigindex<ld4){
							tempsin=cosvals[trigindex+ld4+ld2];
						}else{
							tempsin=cosvals[trigindex-ld4];
						}
						int temp1=i+k;
						int temp2=temp1+increment;
						float flt1=tempcos*r[temp2]-tempsin*im[temp2];
						float flt2=tempcos*im[temp2]+tempsin*r[temp2];
						r[temp2]=r[temp1]-flt1;
						im[temp2]=im[temp1]-flt2;
						r[temp1]+=flt1;
						im[temp1]+=flt2;
					}
				}
				increment+=increment;
			}
			for(int i=0;i<length;i++){
				r[i]/=length;
				im[i]/=length;
			}
		}
	}

	public void dopo4fft(double[] r,double[] im,boolean inverse,boolean bitrev){
	}

	private void initrevbits(){
		bitorder=fftutils.get_bitrev_order(length,kernelsize);
	}

	private void initcosvals(){
		cosvals=new float[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

	private void initkernel(){
		c=new float[kernelsize];
		s=new float[kernelsize];
		for(int i=0;i<kernelsize;i++){
			c[i]=(float)Math.cos(2.0*Math.PI*i/kernelsize);
			s[i]=(float)Math.sin(2.0*Math.PI*i/kernelsize);
		}
	}

}
