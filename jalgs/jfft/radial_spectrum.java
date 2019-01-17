/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.bessel;

public class radial_spectrum{
	public po4realfft2D fft;
	public bessel b;
	public int size;

	public radial_spectrum(int size){
		this.size=size;
		fft=new po4realfft2D(size,size);
		b=new bessel();
	}

	public float[] get_spectrum(float[] image,float minr,float maxr,float dr){
		int nr=(int)((maxr-minr)/dr);
		float[] spectrum=new float[nr];
		float[] rimage=new float[size*size];
		float[] iimage=new float[size*size];
		System.arraycopy(image,0,rimage,0,size*size);
		fft.dorealfft2D(rimage,iimage,false);
		int counter=0;
		for(double r=minr;r<maxr;r+=dr){
			for(int i=0;i<size;i++){
				for(int j=0;j<size;j++){
					double k=2.0*Math.PI*Math.sqrt(j*j+i*i)/size;
					int t=j+i*size;
					// spectrum[counter]+=(rimage[t]*rimage[t]+iimage[t]*iimage[t])*(float)(r*b.besselval(k*r,0));
					if(k>0.0){
						spectrum[counter]+=(rimage[t]*rimage[t]+iimage[t]*iimage[t])*(float)(r*b.besselval(k*r,1)/k);
					}
				}
			}
			if(r>0.0){
				spectrum[counter]/=(float)r;
			}
			counter++;
		}
		return spectrum;
	}

	public float[] get_spectrum(float[] image){
		return get_spectrum(image,0.0f,size/2,1.0f);
	}

}
