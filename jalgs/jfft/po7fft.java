/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po7fft extends po4fft implements Cloneable{
	// public int[] bitorder;
	// public float[] cosvals;
	// public int length;
	// this needs more work
	public int fftindex=3;
	public static float[] c={1.0f,0.623490f,-0.222521f,-0.900969f,-0.900969f,-0.222521f,0.623490f};
	public static float[] s={0.0f,0.781831f,0.974928f,0.433884f,-0.433884f,-0.974928f,-0.781831f};

	/*
	 * here is a basic complex power of 2/7 complex fft (7,14,28,56,...) once
	 * the po5fft class has been created it can be used again for other ffts of
	 * the same size data must be correct length note that bit reversal is made
	 * optional to allow for the calling class to do the bit reversal ahead of
	 * time. This can save extra copying If the bitreversal is done ahead of
	 * time, it is wise to initialize this class with the current bitorder array
	 * to save memory this same logic applies to all of the fft classes
	 */

	public po7fft(int size){
		super();
		length=size;
		initrevbits();
		initcosvals();
	}

	public po7fft(int size,int[] bitorder1,float[] cosvals1){
		length=size;
		bitorder=bitorder1;
		cosvals=cosvals1;
	}

	public po7fft(int size,float[] cosvals1){
		length=size;
		cosvals=cosvals1;
		initrevbits();
	}
	
	public Object clone(){
		return new po7fft(length,bitorder.clone(),cosvals.clone());
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
		rt=new float[7];
		imt=new float[7];
		// do the initial five point transforms

		if(!inverse){
			for(int i=0;i<length;i+=7){
				for(int j=0;j<7;j++){
					rt[j]=r[i+j];
					imt[j]=im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3]+rt[4]+rt[5]+rt[6];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3]+imt[4]+imt[5]+imt[6];
				r[i+1]=rt[0]+c[1]*rt[1]+c[2]*rt[2]+c[3]*rt[3]+c[4]*rt[4]+c[5]*rt[5]+c[6]*rt[6]+0.0f+s[1]*imt[1]+s[2]*imt[2]+s[3]*imt[3]+s[4]*imt[4]+s[5]*imt[5]+s[6]*imt[6];
				im[i+1]=0.0f-s[1]*rt[1]-s[2]*rt[2]-s[3]*rt[3]-s[4]*rt[4]-s[5]*rt[5]-s[6]*rt[6]+imt[0]+c[1]*imt[1]+c[2]*imt[2]+c[3]*imt[3]+c[4]*imt[4]+c[5]*imt[5]+c[6]*imt[6];
				r[i+2]=rt[0]+c[2]*rt[1]+c[4]*rt[2]+c[6]*rt[3]+c[1]*rt[4]+c[3]*rt[5]+c[5]*rt[6]+0.0f+s[2]*imt[1]+s[4]*imt[2]+s[6]*imt[3]+s[1]*imt[4]+s[3]*imt[5]+s[5]*imt[6];
				im[i+2]=0.0f-s[2]*rt[1]-s[4]*rt[2]-s[6]*rt[3]-s[1]*rt[4]-s[3]*rt[5]-s[5]*rt[6]+imt[0]+c[2]*imt[1]+c[4]*imt[2]+c[6]*imt[3]+c[1]*imt[4]+c[3]*imt[5]+c[5]*imt[6];
				r[i+3]=rt[0]+c[3]*rt[1]+c[6]*rt[2]+c[2]*rt[3]+c[5]*rt[4]+c[1]*rt[5]+c[4]*rt[6]+0.0f+s[3]*imt[1]+s[6]*imt[2]+s[2]*imt[3]+s[5]*imt[4]+s[1]*imt[5]+s[4]*imt[6];
				im[i+3]=0.0f-s[3]*rt[1]-s[6]*rt[2]-s[2]*rt[3]-s[5]*rt[4]-s[1]*rt[5]-s[4]*rt[6]+imt[0]+c[3]*imt[1]+c[6]*imt[2]+c[2]*imt[3]+c[5]*imt[4]+c[1]*imt[5]+c[4]*imt[6];
				r[i+4]=rt[0]+c[3]*rt[1]+c[6]*rt[2]+c[2]*rt[3]+c[5]*rt[4]+c[1]*rt[5]+c[4]*rt[6]-0.0f-s[3]*imt[1]-s[6]*imt[2]-s[2]*imt[3]-s[5]*imt[4]-s[1]*imt[5]-s[4]*imt[6];
				im[i+4]=0.0f+s[3]*rt[1]+s[6]*rt[2]+s[2]*rt[3]+s[5]*rt[4]+s[1]*rt[5]+s[4]*rt[6]+imt[0]+c[3]*imt[1]+c[6]*imt[2]+c[2]*imt[3]+c[5]*imt[4]+c[1]*imt[5]+c[4]*imt[6];
				r[i+5]=rt[0]+c[2]*rt[1]+c[4]*rt[2]+c[6]*rt[3]+c[1]*rt[4]+c[3]*rt[5]+c[5]*rt[6]-0.0f-s[2]*imt[1]-s[4]*imt[2]-s[6]*imt[3]-s[1]*imt[4]-s[3]*imt[5]-s[5]*imt[6];
				im[i+5]=0.0f+s[2]*rt[1]+s[4]*rt[2]+s[6]*rt[3]+s[1]*rt[4]+s[3]*rt[5]+s[5]*rt[6]+imt[0]+c[2]*imt[1]+c[4]*imt[2]+c[6]*imt[3]+c[1]*imt[4]+c[3]*imt[5]+c[5]*imt[6];
				r[i+6]=rt[0]+c[1]*rt[1]+c[2]*rt[2]+c[3]*rt[3]+c[4]*rt[4]+c[5]*rt[5]+c[6]*rt[6]-0.0f-s[1]*imt[1]-s[2]*imt[2]-s[3]*imt[3]-s[4]*imt[4]-s[5]*imt[5]-s[6]*imt[6];
				im[i+6]=0.0f+s[1]*rt[1]+s[2]*rt[2]+s[3]*rt[3]+s[4]*rt[4]+s[5]*rt[5]+s[6]*rt[6]+imt[0]+c[1]*imt[1]+c[2]*imt[2]+c[3]*imt[3]+c[4]*imt[4]+c[5]*imt[5]+c[6]*imt[6];
			}
		}else{
			for(int i=0;i<length;i+=7){
				for(int j=0;j<7;j++){
					rt[j]=r[i+j];
					imt[j]=im[i+j];
				}
				r[i]=rt[0]+rt[1]+rt[2]+rt[3]+rt[4]+rt[5]+rt[6];
				im[i]=imt[0]+imt[1]+imt[2]+imt[3]+imt[4]+imt[5]+imt[6];
				r[i+1]=rt[0]+c[1]*rt[1]+c[2]*rt[2]+c[3]*rt[3]+c[4]*rt[4]+c[5]*rt[5]+c[6]*rt[6]-0.0f-s[1]*imt[1]-s[2]*imt[2]-s[3]*imt[3]-s[4]*imt[4]-s[5]*imt[5]-s[6]*imt[6];
				im[i+1]=0.0f+s[1]*rt[1]+s[2]*rt[2]+s[3]*rt[3]+s[4]*rt[4]+s[5]*rt[5]+s[6]*rt[6]+imt[0]+c[1]*imt[1]+c[2]*imt[2]+c[3]*imt[3]+c[4]*imt[4]+c[5]*imt[5]+c[6]*imt[6];
				r[i+2]=rt[0]+c[2]*rt[1]+c[4]*rt[2]+c[6]*rt[3]+c[1]*rt[4]+c[3]*rt[5]+c[5]*rt[6]-0.0f-s[2]*imt[1]-s[4]*imt[2]-s[6]*imt[3]-s[1]*imt[4]-s[3]*imt[5]-s[5]*imt[6];
				im[i+2]=0.0f+s[2]*rt[1]+s[4]*rt[2]+s[6]*rt[3]+s[1]*rt[4]+s[3]*rt[5]+s[5]*rt[6]+imt[0]+c[2]*imt[1]+c[4]*imt[2]+c[6]*imt[3]+c[1]*imt[4]+c[3]*imt[5]+c[5]*imt[6];
				r[i+3]=rt[0]+c[3]*rt[1]+c[6]*rt[2]+c[2]*rt[3]+c[5]*rt[4]+c[1]*rt[5]+c[4]*rt[6]-0.0f-s[3]*imt[1]-s[6]*imt[2]-s[2]*imt[3]-s[5]*imt[4]-s[1]*imt[5]-s[4]*imt[6];
				im[i+3]=0.0f+s[3]*rt[1]+s[6]*rt[2]+s[2]*rt[3]+s[5]*rt[4]+s[1]*rt[5]+s[4]*rt[6]+imt[0]+c[3]*imt[1]+c[6]*imt[2]+c[2]*imt[3]+c[5]*imt[4]+c[1]*imt[5]+c[4]*imt[6];
				r[i+4]=rt[0]+c[3]*rt[1]+c[6]*rt[2]+c[2]*rt[3]+c[5]*rt[4]+c[1]*rt[5]+c[4]*rt[6]+0.0f+s[3]*imt[1]+s[6]*imt[2]+s[2]*imt[3]+s[5]*imt[4]+s[1]*imt[5]+s[4]*imt[6];
				im[i+4]=0.0f-s[3]*rt[1]-s[6]*rt[2]-s[2]*rt[3]-s[5]*rt[4]-s[1]*rt[5]-s[4]*rt[6]+imt[0]+c[3]*imt[1]+c[6]*imt[2]+c[2]*imt[3]+c[5]*imt[4]+c[1]*imt[5]+c[4]*imt[6];
				r[i+5]=rt[0]+c[2]*rt[1]+c[4]*rt[2]+c[6]*rt[3]+c[1]*rt[4]+c[3]*rt[5]+c[5]*rt[6]+0.0f+s[2]*imt[1]+s[4]*imt[2]+s[6]*imt[3]+s[1]*imt[4]+s[3]*imt[5]+s[5]*imt[6];
				im[i+5]=0.0f-s[2]*rt[1]-s[4]*rt[2]-s[6]*rt[3]-s[1]*rt[4]-s[3]*rt[5]-s[5]*rt[6]+imt[0]+c[2]*imt[1]+c[4]*imt[2]+c[6]*imt[3]+c[1]*imt[4]+c[3]*imt[5]+c[5]*imt[6];
				r[i+6]=rt[0]+c[1]*rt[1]+c[2]*rt[2]+c[3]*rt[3]+c[4]*rt[4]+c[5]*rt[5]+c[6]*rt[6]+0.0f+s[1]*imt[1]+s[2]*imt[2]+s[3]*imt[3]+s[4]*imt[4]+s[5]*imt[5]+s[6]*imt[6];
				im[i+6]=0.0f-s[1]*rt[1]-s[2]*rt[2]-s[3]*rt[3]-s[4]*rt[4]-s[5]*rt[5]-s[6]*rt[6]+imt[0]+c[1]*imt[1]+c[2]*imt[2]+c[3]*imt[3]+c[4]*imt[4]+c[5]*imt[5]+c[6]*imt[6];
			}
		}
		rt=null;
		imt=null;
		// now combine these into ever increasing size transforms to get the
		// final fft
		int increment=7;
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
		bitorder=fftutils.get_bitrev_order(length,7);
	}

	private void initcosvals(){
		cosvals=new float[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

}
