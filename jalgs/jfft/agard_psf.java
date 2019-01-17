/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

import jalgs.bessel;

public class agard_psf{
	//here we implement the algorithm of Agard 1984 Ann. Rev. Biophys. Bioeng. vol 13, pp 191-219.

	public static String[] manufacturers={"Zeiss","Olympus","Leica","Nikon"};
	public static float[] tllengths={165000.0f,180000.0f,200000.0f,200000.0f};
	public int maxsize;
	public int halfsize=128;
	public float halfsize2;
	public float res=0.05f;
	public float[] rzpsf,ctf2D;

	public agard_psf(double mag,double NA,double lem,double n,int man){
		double rresolution=res;
		maxsize=halfsize*2;
		halfsize2=res*(float)halfsize;
		double qresolution=1.0/(rresolution*maxsize);
		double zresolution=res;
		double alpha=Math.asin(NA/n);
		double fc=(2.0*NA)/lem;
		int fclim=(int)(fc/qresolution);
		if(fclim<(fc/qresolution)){
			fclim++;
		}
		double tllength=tllengths[man];
		double df=tllength/mag;
		bessel b=new bessel();
		po4realfft fft=new po4realfft(halfsize*2);
		rzpsf=new float[halfsize*halfsize];
		ctf2D=new float[halfsize*halfsize];
		for(int i=0;i<halfsize;i++){
			float[] ctf=new float[halfsize];
			double z=zresolution*i;
			double temp=Math.cos(alpha);
			double w=-df-z*temp+Math.sqrt(df*df+2.0*df*z+z*z*temp*temp);
			for(int j=0;j<fclim;j++){
				double q=qresolution*j;
				double temp2=(8.0*Math.PI*w*(1.0-q/fc)*(q/fc))/lem;
				double temp3;
				if(temp2!=0.0){
					temp=temp2;
					temp3=(2.0*b.besselval(temp2,1))/temp;
				}else{
					temp3=1.0;
				}
				double beta=Math.acos(q/fc);
				ctf[j]=(float)(((2.0*beta-Math.sin(2.0*beta))*temp3)/Math.PI);
				ctf2D[j+i*halfsize]=ctf[j];
			}
			// the psf is just the zero order hankel transfrom of the 2D ctf
			double twopi=2.0*Math.PI;
			for(int j=0;j<halfsize;j++){
				double r=rresolution*j;
				for(int k=0;k<fclim;k++){
					double q=qresolution*k;
					rzpsf[j+i*halfsize]+=(float)(ctf[k]*q*b.besselval(twopi*q*r,0));
				}
			}
		}
		for(int i=0;i<fclim;i++){
			float[] temp=new float[halfsize*2];
			temp[0]=ctf2D[i];
			for(int j=1;j<halfsize;j++){
				temp[j]=ctf2D[i+j*halfsize];
				temp[2*halfsize-j]=temp[j];
			}
			fft.realfft(temp,new float[2*halfsize],false);
			for(int j=0;j<halfsize;j++){
				ctf2D[i+j*halfsize]=temp[j];
			}
		}
	}

	/*****************
	 * gets the psf value at a specified x, y, and z distance from the psf center
	 * @param x: x distance from the psf center in microns
	 * @param y: y distance from the psf center in microns
	 * @param z: axial distance from the psf center in microns
	 * @return
	 */
	public float getvalue(float x,float y,float z){
		return getvalue((float)Math.sqrt(x*x+y*y),z);
	}
	
	/**********************
	 * gets the relative radial psf value
	 * @param r: radial distance from the psf center in microns
	 * @param z: axial distance from the psf center in microns
	 * @return
	 */
	public float getvalue(float r,float z){
		float r2=r/res; //radial distance now in pixels
		float z2=z/res; //axial distance now in pixels
		if(r2<0.0f) r2=-r2; //make them positive
		if(z2<0.0f) z2=-z2; //make them positive
		if(r2>=(halfsize-1)) return 0.0f;
		if(z2>=(halfsize-1)) return 0.0f;
		int prevr=(int)r2;
		float rrem=r2-prevr;
		int prevz=(int)z2;
		float prevval=rzpsf[prevr+prevz*halfsize]+rrem*(rzpsf[prevr+1+prevz*halfsize]-rzpsf[prevr+prevz*halfsize]);
		float nextval=rzpsf[prevr+(prevz+1)*halfsize]+rrem*(rzpsf[prevr+1+(prevz+1)*halfsize]-rzpsf[prevr+(prevz+1)*halfsize]);
		float zrem=z2-prevz;
		return prevval+zrem*(nextval-prevval);
	}

	/**********************
	 * this draws the xy cross section of the psf at zpos distance from the center
	 * @param image: the image to draw on
	 * @param width
	 * @param height
	 * @param x: the psf center in x in microns
	 * @param y: the psf center in y in microns
	 * @param dx: pixel size
	 * @param zpos: distance from psf center in microns
	 * @param maxsize: our cutoff distance for the edge of the psf in microns
	 * @param intensity: the intensity multiplier for the psf
	 */
	public void draw_psf(float[] image,int width,int height,float x,float y,float dx,float zpos,float maxsize,float intensity){
		if(intensity==0.0f) return;
		if(Float.isNaN(intensity)) return;
		float maxsize2=maxsize;
		if(maxsize2>halfsize2) maxsize2=halfsize2;
		int xstart=(int)((x-0.5f*maxsize2)/dx);
		if(xstart<0) xstart=0;
		if(xstart>(width-1)) return;
		int xend=1+(int)((x+0.5f*maxsize2)/dx);
		if(xend>(width-1)) xend=width-1;
		if(xend<0) return;
		int ystart=(int)((y-0.5f*maxsize2)/dx);
		if(ystart<0) ystart=0;
		if(ystart>(height-1)) return;
		int yend=1+(int)((y+0.5f*maxsize2)/dx);
		if(yend>(height-1)) yend=height-1;
		if(xend<0) return;
		for(int i=ystart;i<yend;i++){
			float ypos=dx*i-y;
			for(int j=xstart;j<xend;j++){
				float xpos=dx*j-x;
				float temp=getvalue(xpos,ypos,zpos);
				if(!Float.isNaN(temp)) image[j+i*width]+=intensity*temp; //why would I get nan values here?
				//image[j+i*width]+=intensity*temp;
			}
		}
	}

	/**********************
	 * 
	 * @param stack: the image stack, must be float arrays
	 * @param width
	 * @param height
	 * @param slices
	 * @param x: the psf center in x (in microns)
	 * @param y: the psf center in y (in microns)
	 * @param z: the psf center in z (in microns);
	 * @param dx: xy pixel size of the stack
	 * @param dz: z spacing of the stack
	 * @param maxsizex: our cutoff for the edge of the psf (in microns)
	 * @param maxsizez: our cutoff for the edge of the psf (in microns)
	 * @param intensity: the intensity multiplier for the psf
	 */
	public void draw_psf(Object[] stack,int width,int height,int slices,float x,float y,float z,float dx,float dz,float maxsizex,float maxsizez,float intensity){
		if(intensity==0.0f) return;
		if(Float.isNaN(intensity)) return;
		float maxsizez2=maxsizez;
		if(maxsizez2>halfsize2) maxsizez2=halfsize2;
		int zstart=(int)((z-0.5f*maxsizez2)/dz);
		if(zstart<0) zstart=0;
		if(zstart>(slices-1)) return;
		int zend=(int)((z+0.5f*maxsizez2)/dz);
		if(zend>(slices-1)) zend=slices-1;
		if(zend<0) return;
		for(int i=zstart;i<zend;i++){
			float zpos=dz*i-z;
			draw_psf((float[])stack[i],width,height,x,y,dx,zpos,maxsizex,intensity);
		}
	}

}
