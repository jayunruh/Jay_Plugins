/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4cfft3D{
	public po4fft cfft;
	public po4cfft2D cfft2D;
	public float[] cosvals;
	public int length1,length2,length3;
	public int[] bitorder;

	/*
	 * This class does complex 3D ffts once the po4cfft3D class has been created
	 * it can be used again for other 3D ffts of the same size width, height,
	 * and number of slices must be power of two size Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */

	public po4cfft3D(int length11,int length12,int length13){
		length1=length11;
		length2=length12;
		length3=length13;
		cfft2D=new po4cfft2D(length1,length2);
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
				initrevbits();
			}
		}
		cfft=new po4fft(length3,bitorder,cosvals);
	}

	public void docfft3D(Object[] real,Object[] im,boolean inverse){
		// the input objects arrays must contain 1D float arrays representing
		// each z "slice"
		if(!inverse){
			// start by doing the z ffts
			// copy the data into arrays in bit reversed order
			for(int i=0;i<length1*length2;i++){
				float[] tempfloat=new float[length3];
				float[] tempfloat2=new float[length3];
				for(int j=0;j<length3;j++){
					tempfloat[j]=((float[])real[bitorder[j]])[i];
					tempfloat2[j]=((float[])im[bitorder[j]])[i];
				}
				cfft.dopo4fft(tempfloat,tempfloat2,false,false);
				for(int j=0;j<length3;j++){
					((float[])real[j])[i]=tempfloat[j];
					((float[])im[j])[i]=tempfloat2[j];
				}
			}
			// now do the xy ffts (these are complex)
			for(int i=0;i<length3;i++){
				cfft2D.docfft2D((float[])real[i],(float[])im[i],false);
			}
		}else{
			// here start with the xy complex ffts
			for(int i=0;i<length3;i++){
				cfft2D.docfft2D((float[])real[i],(float[])im[i],true);
			}
			// now the z tworeal ffts
			for(int i=0;i<length1*length2;i++){
				float[] tempfloat=new float[length3];
				float[] tempfloat2=new float[length3];
				for(int j=0;j<length3;j++){
					tempfloat[j]=((float[])real[bitorder[j]])[i];
					tempfloat2[j]=((float[])im[bitorder[j]])[i];
				}
				cfft.dopo4fft(tempfloat,tempfloat2,true,false);
				for(int j=0;j<length3;j++){
					((float[])real[j])[i]=tempfloat[j];
					((float[])im[j])[i]=tempfloat2[j];
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

	private void initrevbits(){
		// this has been modified for power of 4 ffts
		bitorder=new int[length3];
		for(int i=0;i<length3;i++){
			bitorder[i]=i;
		}
		int j=0;
		int k;
		for(int i=0;i<length3;i++){
			if(j>i){
				bitorder[j]=i;
				bitorder[i]=j;
			}
			k=length3/2;
			while(k>=1&&j>=k){
				j-=k;
				k/=2;
			}
			j+=k;
		}
		for(int i=0;i<length3;i+=4){
			int temp=bitorder[i+1];
			bitorder[i+1]=bitorder[i+2];
			bitorder[i+2]=temp;
		}
	}

}
