/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4realfftd{
	public po4fftd cfft; // complex fft class used for this calculation
	public double[] cosvals;
	public int length;

	/*
	 * here is a basic real power of 4 fft once the po4realfft class has been
	 * created it can be used again for other ffts of the same size data must be
	 * power of 2 length for the forward transform, the imaginary vector should
	 * be empty for the inverse transform, the imaginary vector will be returned
	 * empty Copyright Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public po4realfftd(int size){
		length=size;
		initcosvals();
		double[] cosvals2=new double[length/2];
		for(int i=0;i<length/2;i++){
			cosvals2[i]=cosvals[i+i];
		}
		cfft=new po4fftd(length/2,cosvals2);
	}

	private void initcosvals(){
		cosvals=new double[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=Math.cos(dumdbl);
		}
	}

	public void realfft(double[] real,double[] im,boolean inverse){
		// this algorithm finds the FFT of the even and odd components using the
		// twofft approach
		// these are then converted to the final FFT using
		// F(n)=F(e,n)+exp(2*pi*i*n/N)*F(o,n)
		int ld2=length/2;
		int ld4=length/4;
		if(!inverse){
			// fft the already "interlaced" data just like in twofft
			cfft.dopo4fft(real,false);
			double[] reven=new double[ld2];
			double[] rodd=new double[ld2];
			double[] ieven=new double[ld2];
			double[] iodd=new double[ld2];
			reven[0]=real[0];
			rodd[0]=real[1];
			reven[ld4]=real[ld2];
			rodd[ld4]=real[ld2+1];
			for(int i=1;i<ld4;i++){
				int i2=i+i;
				reven[i]=0.5*(real[i2]+real[length-i2]);
				reven[ld2-i]=reven[i];
				rodd[i]=0.5*(real[i2+1]+real[length-i2+1]);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5*(real[i2+1]-real[length-i2+1]);
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5*(-real[i2]+real[length-i2]);
				iodd[ld2-i]=-iodd[i];
			}
			ieven[0]=0.0f;
			iodd[0]=0.0f;
			ieven[ld4]=0.0f;
			iodd[ld4]=0.0f;
			// now that we have the even and odd transforms, get the final
			// transform
			real[0]=reven[0]+rodd[0];
			im[0]=0.0f;
			real[ld2]=reven[0]-rodd[0];
			for(int i=1;i<ld2;i++){
				double tempsin=-cosvals[i+ld4];
				real[i]=reven[i]+cosvals[i]*rodd[i]+tempsin*iodd[i];
				real[length-i]=real[i];
				im[i]=ieven[i]+cosvals[i]*iodd[i]-tempsin*rodd[i];
				im[length-i]=-im[i];
			}
		}else{
			// start by obtaining the even and odd transforms
			double[] reven=new double[ld2];
			double[] rodd=new double[ld2];
			double[] ieven=new double[ld2];
			double[] iodd=new double[ld2];
			reven[0]=0.5*(real[0]+real[ld2]);
			reven[ld4]=real[ld4];
			rodd[0]=0.5*(real[0]-real[ld2]);
			rodd[ld4]=-im[ld4];
			for(int i=1;i<ld4;i++){
				double tempsin=-cosvals[i+ld4];
				double rpl=real[i]+real[i+ld2];
				double rmn=real[i]-real[i+ld2];
				double ipl=im[i]+im[i+ld2];
				double imn=im[i]-im[i+ld2];
				reven[i]=0.5*rpl;
				reven[ld2-i]=reven[i];
				rodd[i]=0.5*(cosvals[i]*rmn-tempsin*imn);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5*ipl;
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5*(cosvals[i]*imn+tempsin*rmn);
				iodd[ld2-i]=-iodd[i];
			}
			// combine them into the appropriate complex data set
			for(int i=0;i<ld2;i++){
				real[2*i]=reven[i]-iodd[i];
				real[2*i+1]=rodd[i]+ieven[i];
			}
			cfft.dopo4fft(real,true);
		}
	}

}
