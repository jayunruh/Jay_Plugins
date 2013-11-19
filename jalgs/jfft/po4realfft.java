/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class po4realfft{
	public po4fft cfft; // complex fft class used for this calculation
	public float[] cosvals;
	public int length;
	public int fftindex;

	/*
	 * here is a basic real power of 4 fft once the po4realfft class has been
	 * created it can be used again for other ffts of the same size data must be
	 * power of 2 length for the forward transform, the imaginary vector should
	 * be empty for the inverse transform, the imaginary vector will be returned
	 * empty Copyright Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */

	public po4realfft(int size){
		length=size;
		initcosvals();
		float[] cosvals2=new float[length/2];
		for(int i=0;i<length/2;i++){
			cosvals2[i]=cosvals[i+i];
		}
		cfft=new po4fft(length/2,cosvals2);
	}

	public po4realfft(int size,int index){
		// here we can use an alternative cfft
		// indices are 0:po4 1:po5 2:po6 3:po7
		// only works with even length data sets
		length=size;
		initcosvals();
		float[] cosvals2=new float[length/2];
		for(int i=0;i<length/2;i++){
			cosvals2[i]=cosvals[i+i];
		}
		fftindex=index;
		cfft=fftutils.construct_fft(length/2,fftindex,cosvals2);
	}

	private void initcosvals(){
		cosvals=new float[length];
		for(int i=0;i<length;i++){
			double dumdbl=Math.PI*2.0*((double)i/(double)length);
			cosvals[i]=(float)Math.cos(dumdbl);
		}
	}

	public void realfft(float[] real,float[] im,boolean inverse){
		// this algorithm finds the FFT of the even and odd components using the
		// twofft approach
		// these are then converted to the final FFT using
		// F(n)=F(e,n)+exp(2*pi*i*n/N)*F(o,n)
		int ld2=length/2;
		int ld4=length/4;
		if(!inverse){
			// fft the already "interlaced" data just like in twofft
			cfft.dopo4fft(real,false);
			float[] reven=new float[ld2];
			float[] rodd=new float[ld2];
			float[] ieven=new float[ld2];
			float[] iodd=new float[ld2];
			reven[0]=real[0];
			rodd[0]=real[1];
			reven[ld4]=real[ld2];
			rodd[ld4]=real[ld2+1];
			for(int i=1;i<ld4;i++){
				int i2=i+i;
				reven[i]=0.5f*(real[i2]+real[length-i2]);
				reven[ld2-i]=reven[i];
				rodd[i]=0.5f*(real[i2+1]+real[length-i2+1]);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5f*(real[i2+1]-real[length-i2+1]);
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5f*(-real[i2]+real[length-i2]);
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
				float tempsin=-cosvals[i+ld4];
				real[i]=reven[i]+cosvals[i]*rodd[i]+tempsin*iodd[i];
				real[length-i]=real[i];
				im[i]=ieven[i]+cosvals[i]*iodd[i]-tempsin*rodd[i];
				im[length-i]=-im[i];
			}
		}else{
			// start by obtaining the even and odd transforms
			float[] reven=new float[ld2];
			float[] rodd=new float[ld2];
			float[] ieven=new float[ld2];
			float[] iodd=new float[ld2];
			reven[0]=0.5f*(real[0]+real[ld2]);
			reven[ld4]=real[ld4];
			rodd[0]=0.5f*(real[0]-real[ld2]);
			rodd[ld4]=-im[ld4];
			for(int i=1;i<ld4;i++){
				float tempsin=-cosvals[i+ld4];
				float rpl=real[i]+real[i+ld2];
				float rmn=real[i]-real[i+ld2];
				float ipl=im[i]+im[i+ld2];
				float imn=im[i]-im[i+ld2];
				reven[i]=0.5f*rpl;
				reven[ld2-i]=reven[i];
				rodd[i]=0.5f*(cosvals[i]*rmn-tempsin*imn);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5f*ipl;
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5f*(cosvals[i]*imn+tempsin*rmn);
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
			float[] reven=new float[ld2];
			float[] rodd=new float[ld2];
			float[] ieven=new float[ld2];
			float[] iodd=new float[ld2];
			reven[0]=(float)real[0];
			rodd[0]=(float)real[1];
			reven[ld4]=(float)real[ld2];
			rodd[ld4]=(float)real[ld2+1];
			for(int i=1;i<ld4;i++){
				int i2=i+i;
				reven[i]=0.5f*((float)real[i2]+(float)real[length-i2]);
				reven[ld2-i]=reven[i];
				rodd[i]=0.5f*((float)real[i2+1]+(float)real[length-i2+1]);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5f*((float)real[i2+1]-(float)real[length-i2+1]);
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5f*(-(float)real[i2]+(float)real[length-i2]);
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
				float tempsin=-cosvals[i+ld4];
				real[i]=reven[i]+cosvals[i]*rodd[i]+tempsin*iodd[i];
				real[length-i]=real[i];
				im[i]=ieven[i]+cosvals[i]*iodd[i]-tempsin*rodd[i];
				im[length-i]=-im[i];
			}
		}else{
			// start by obtaining the even and odd transforms
			float[] reven=new float[ld2];
			float[] rodd=new float[ld2];
			float[] ieven=new float[ld2];
			float[] iodd=new float[ld2];
			reven[0]=0.5f*((float)real[0]+(float)real[ld2]);
			reven[ld4]=(float)real[ld4];
			rodd[0]=0.5f*((float)real[0]-(float)real[ld2]);
			rodd[ld4]=-(float)im[ld4];
			for(int i=1;i<ld4;i++){
				float tempsin=-cosvals[i+ld4];
				float rpl=(float)real[i]+(float)real[i+ld2];
				float rmn=(float)real[i]-(float)real[i+ld2];
				float ipl=(float)im[i]+(float)im[i+ld2];
				float imn=(float)im[i]-(float)im[i+ld2];
				reven[i]=0.5f*rpl;
				reven[ld2-i]=reven[i];
				rodd[i]=0.5f*(cosvals[i]*rmn-tempsin*imn);
				rodd[ld2-i]=rodd[i];
				ieven[i]=0.5f*ipl;
				ieven[ld2-i]=-ieven[i];
				iodd[i]=0.5f*(cosvals[i]*imn+tempsin*rmn);
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
