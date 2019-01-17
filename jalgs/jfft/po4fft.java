/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4fft implements Cloneable{
	public int[] bitorder;
	public float[] cosvals;
	public int length;
	public int fftindex=0;

	/*
	 * here is a basic complex power of 4 complex fft (actually powers of 2)
	 * once the po4fft class has been created it can be used again for other
	 * ffts of the same size data must be power of 2 length note that bit
	 * reversal is made optional to allow for the calling class to do the bit
	 * reversal ahead of time. This can save extra copying If the bitreversal is
	 * done ahead of time, it is wise to initialize this class with the current
	 * bitorder array to save memory this same logic applies to all of the fft
	 * classes
	 */

	public po4fft(){
	} // empty constructor

	public po4fft(int size){
		length=size;
		initrevbits();
		initcosvals();
	}

	public po4fft(int size,int[] bitorder1,float[] cosvals1){
		length=size;
		bitorder=bitorder1;
		cosvals=cosvals1;
	}

	public po4fft(int size,float[] cosvals1){
		length=size;
		cosvals=cosvals1;
		initrevbits();
	}
	
	public Object clone(){
		return new po4fft(length,bitorder.clone(),cosvals.clone());
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
		// here is a power of 4 fft
		// split the data into almost bit reversed real and imaginary arrays
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
		rt=new float[4];
		imt=new float[4];
		// do the initial four point transforms
		if(!inverse){
			for(int i=0;i<length;i+=4){
				for(int j=0;j<4;j++){
					rt[j]=r[i+j];
					imt[j]=im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3];
				r[i+1]=rt[0]-rt[2]+imt[1]-imt[3];
				im[i+1]=rt[3]-rt[1]+imt[0]-imt[2];
				r[i+2]=rt[0]-rt[1]+rt[2]-rt[3];
				im[i+2]=imt[0]-imt[1]+imt[2]-imt[3];
				r[i+3]=rt[0]-rt[2]-imt[1]+imt[3];
				im[i+3]=-rt[3]+rt[1]+imt[0]-imt[2];
			}
		}else{
			for(int i=0;i<length;i+=4){
				for(int j=0;j<4;j++){
					rt[j]=r[i+j];
					imt[j]=im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3];
				r[i+1]=rt[0]-rt[2]-imt[1]+imt[3];
				im[i+1]=-rt[3]+rt[1]+imt[0]-imt[2];
				r[i+2]=rt[0]-rt[1]+rt[2]-rt[3];
				im[i+2]=imt[0]-imt[1]+imt[2]-imt[3];
				r[i+3]=rt[0]-rt[2]+imt[1]-imt[3];
				im[i+3]=rt[3]-rt[1]+imt[0]-imt[2];
			}
		}
		rt=null;
		imt=null;
		// now combine these into ever increasing size transforms to get the
		// final fft
		int increment=4;
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
		// here is a power of 4 fft
		// split the data into almost bit reversed real and imaginary arrays
		int ld2=length/2;
		int ld4=length/4;
		float[] rt,imt;
		if(bitrev){
			rt=new float[length];
			imt=new float[length];
			for(int i=0;i<length;i++){
				rt[i]=(float)r[i];
				imt[i]=(float)im[i];
			}
			for(int i=0;i<length;i++){
				r[i]=rt[bitorder[i]];
				im[i]=imt[bitorder[i]];
			}
		}
		rt=new float[4];
		imt=new float[4];
		// do the initial four point transforms
		if(!inverse){
			for(int i=0;i<length;i+=4){
				for(int j=0;j<4;j++){
					rt[j]=(float)r[i+j];
					imt[j]=(float)im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3];
				r[i+1]=rt[0]-rt[2]+imt[1]-imt[3];
				im[i+1]=rt[3]-rt[1]+imt[0]-imt[2];
				r[i+2]=rt[0]-rt[1]+rt[2]-rt[3];
				im[i+2]=imt[0]-imt[1]+imt[2]-imt[3];
				r[i+3]=rt[0]-rt[2]-imt[1]+imt[3];
				im[i+3]=-rt[3]+rt[1]+imt[0]-imt[2];
			}
		}else{
			for(int i=0;i<length;i+=4){
				for(int j=0;j<4;j++){
					rt[j]=(float)r[i+j];
					imt[j]=(float)im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3];
				r[i+1]=rt[0]-rt[2]-imt[1]+imt[3];
				im[i+1]=-rt[3]+rt[1]+imt[0]-imt[2];
				r[i+2]=rt[0]-rt[1]+rt[2]-rt[3];
				im[i+2]=imt[0]-imt[1]+imt[2]-imt[3];
				r[i+3]=rt[0]-rt[2]+imt[1]-imt[3];
				im[i+3]=rt[3]-rt[1]+imt[0]-imt[2];
			}
		}
		rt=null;
		imt=null;
		// now combine these into ever increasing size transforms to get the
		// final fft
		int increment=4;
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
						float flt1=tempcos*(float)r[temp2]+tempsin*(float)im[temp2];
						float flt2=tempcos*(float)im[temp2]-tempsin*(float)r[temp2];
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
						float flt1=tempcos*(float)r[temp2]-tempsin*(float)im[temp2];
						float flt2=tempcos*(float)im[temp2]+tempsin*(float)r[temp2];
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

	private void initrevbits(){
		bitorder=fftutils.get_bitrev_order(length,4);
	}

	private void initcosvals(){
		cosvals=new float[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

}
