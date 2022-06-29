/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.MemoryImageSource;

public class sim_calcintensity{
	rngs random;
	double[] model_gaus;
	double stdev,ratio;
	int maxdist,maxdistz;
	sim_setup set;

	public sim_calcintensity(sim_setup set1){
		set=set1;
		random=set.random;
		stdev=0.5*set.w0pixels;
		ratio=set.z0/set.w0;
		maxdist=(int)(2.0f*set.w0pixels);
		maxdistz=(int)(2.0f*set.z0pixels);
		init_gaus();
	}

	public void change_setup(sim_setup set1){
		set=set1;
		stdev=0.5*set.w0pixels;
		ratio=set.z0/set.w0;
		maxdist=(int)(2.0f*set.w0pixels);
		maxdistz=(int)(2.0f*set.z0pixels);
		init_gaus();
	}

	public float[][] calc_image(sim_multispecies sm,float[] back){
		// here we calculate an image at once from the molecule data
		// multiple images are calculated for multiple channels
		int nchannels=sm.sslist[0].bright2.length;
		sim_particle[] sps=sm.sps;
		float[][] fintensity=new float[nchannels][set.boxpixels*set.boxpixels];

		for(int i=0;i<sps.length;i++){
			boolean bleached=true;
			for(int j=0;j<nchannels;j++)
				if(sps[i].bleach_state[j]!=0){
					bleached=false;
					break;
				}
			if(!bleached){
				double zr=0.0;
				float multiplierz=1.0f;
				if(set.confineindex!=1&&set.confineindex!=3){
					zr=Math.abs(sps[i].coords[2]-set.fboxpixels/2.0f);
				}
				if(zr<=maxdistz){
					if(set.confineindex!=1&&set.confineindex!=3){
						multiplierz=(float)getinterpgaus(zr/ratio);
					}
					int starty=(int)sps[i].coords[1]-maxdist;
					int endy=(int)sps[i].coords[1]+maxdist;
					for(int pointy=starty;pointy<endy;pointy++){
						double yr=Math.abs(sps[i].coords[1]-pointy);
						float multiplier=(float)getinterpgaus(yr);
						int startx=(int)sps[i].coords[0]-maxdist;
						int endx=(int)sps[i].coords[0]+maxdist;
						int pointy2=pointy;
						if(pointy2<0){
							pointy2+=set.boxpixels;
						}
						if(pointy2>=set.boxpixels){
							pointy2-=set.boxpixels;
						}
						for(int pointx=startx;pointx<endx;pointx++){
							double xr=Math.abs(sps[i].coords[0]-pointx);
							float multiplierx=(float)getinterpgaus(xr);
							int pointx2=pointx;
							if(pointx2<0){
								pointx2+=set.boxpixels;
							}
							if(pointx2>=set.boxpixels){
								pointx2-=set.boxpixels;
							}
							for(int j=0;j<nchannels;j++){
								float intensity=sps[i].get_brightness(sm,j)*multiplierx*multiplier*multiplierz;
								fintensity[j][pointx2+pointy2*set.boxpixels]+=intensity;
								sps[i].photonsused[j]+=intensity;
							}
						}
					}
				}
			}
		}
		for(int i=0;i<nchannels;i++){
			float tempback=back[i]*set.frametime2;
			for(int j=0;j<set.boxpixels*set.boxpixels;j++){
				fintensity[i][j]+=tempback;
			}
		}
		if(set.noiseindex==2){
			return fintensity;
		}else{
			float[][] photons=new float[nchannels][set.boxpixels*set.boxpixels];
			for(int i=0;i<nchannels;i++){
				for(int j=0;j<set.boxpixels*set.boxpixels;j++){
					photons[i][j]=set.random.poidev(fintensity[i][j]);
				}
			}
			if(set.noiseindex==0){
				return photons;
			}else{
				float[][] aint=new float[nchannels][set.boxpixels*set.boxpixels];
				for(int i=0;i<nchannels;i++){
					for(int j=0;j<set.boxpixels*set.boxpixels;j++){
						aint[i][j]=(float)set.random.analogdev((int)photons[i][j],set.gain,set.offset,set.readstdev);
					}
				}
				return aint;
			}
		}
	}

	public float[][] calc_line(sim_multispecies sm,float[] back){
		return calc_line(sm,back,set.boxpixels/2);
	}

	public float[][] calc_line(sim_multispecies sm,float[] back,int pointy){
		// here we calculate an line at once from the molecule data
		// multiple images are calculated for multiple channels
		int nchannels=sm.sslist[0].bright2.length;
		sim_particle[] sps=sm.sps;
		float[][] fintensity=new float[nchannels][set.boxpixels];

		for(int i=0;i<sps.length;i++){
			boolean bleached=true;
			for(int j=0;j<nchannels;j++)
				if(sps[i].bleach_state[j]!=0){
					bleached=false;
					break;
				}
			if(!bleached){
				double zr=0.0;
				float multiplierz=1.0f;
				if(set.confineindex!=1&&set.confineindex!=3){
					zr=Math.abs(sps[i].coords[2]-set.fboxpixels/2.0f);
				}
				if(zr<=maxdistz){
					if(set.confineindex!=1&&set.confineindex!=3){
						multiplierz=(float)getinterpgaus(zr/ratio);
					}
					int starty=(int)sps[i].coords[1]-maxdist;
					int endy=(int)sps[i].coords[1]+maxdist;
					if(pointy>=starty&&pointy<endy){
						double yr=Math.abs(sps[i].coords[1]-pointy);
						float multiplier=(float)getinterpgaus(yr);
						int startx=(int)sps[i].coords[0]-maxdist;
						int endx=(int)sps[i].coords[0]+maxdist;
						for(int pointx=startx;pointx<endx;pointx++){
							double xr=Math.abs(sps[i].coords[0]-pointx);
							float multiplierx=(float)getinterpgaus(xr);
							int pointx2=pointx;
							if(pointx2<0){
								pointx2+=set.boxpixels;
							}
							if(pointx2>=set.boxpixels){
								pointx2-=set.boxpixels;
							}
							for(int j=0;j<nchannels;j++){
								float intensity=sps[i].get_brightness(sm,j)*multiplierx*multiplier*multiplierz;
								fintensity[j][pointx2]+=intensity;
								sps[i].photonsused[j]+=intensity;
							}
						}
					}
				}
			}
		}
		for(int i=0;i<nchannels;i++){
			float tempback=back[i]*set.frametime2;
			for(int j=0;j<set.boxpixels;j++){
				fintensity[i][j]+=tempback;
			}
		}
		if(set.noiseindex==2){
			return fintensity;
		}else{
			float[][] photons=new float[nchannels][set.boxpixels];
			for(int i=0;i<nchannels;i++){
				for(int j=0;j<set.boxpixels;j++){
					photons[i][j]=set.random.poidev(fintensity[i][j]);
				}
			}
			if(set.noiseindex==0){
				return photons;
			}else{
				float[][] aint=new float[nchannels][set.boxpixels];
				for(int i=0;i<nchannels;i++){
					for(int j=0;j<set.boxpixels;j++){
						aint[i][j]=(float)set.random.analogdev((int)photons[i][j],set.gain,set.offset,set.readstdev);
					}
				}
				return aint;
			}
		}
	}

	public float[] calc_point(sim_multispecies sm,float[] back){
		return calc_point(sm,back,set.boxpixels/2,set.boxpixels/2);
	}

	public float[] calc_point_tpe(sim_multispecies sm,float[] back){
		return calc_point_tpe(sm,back,set.boxpixels/2,set.boxpixels/2);
	}

	public float[] calc_point(sim_multispecies sm,float[] back,int pointx){
		return calc_point(sm,back,pointx,set.boxpixels/2);
	}

	public float[] calc_point_tpe(sim_multispecies sm,float[] back,int pointx){
		return calc_point_tpe(sm,back,pointx,set.boxpixels/2);
	}

	public float[] calc_point(sim_multispecies sm,float[] back,int pointx,int pointy){
		// here we calculate a point from the molecule data
		// multiple points are calculated for multiple channels
		int nchannels=sm.sslist[0].bright2.length;
		sim_particle[] sps=sm.sps;
		float[] fintensity=new float[nchannels];
		for(int i=0;i<sps.length;i++){
			boolean bleached=true;
			for(int j=0;j<nchannels;j++)
				if(sps[i].bleach_state[j]!=0){
					bleached=false;
					break;
				}
			if(!bleached){
				double zr=0.0;
				float multiplierz=1.0f;
				if(set.confineindex!=1&&set.confineindex!=3){
					zr=Math.abs(sps[i].coords[2]-set.fboxpixels/2.0f);
				}
				if(zr<=maxdistz){
					if(set.confineindex!=1&&set.confineindex!=3){
						multiplierz=(float)getinterpgaus(zr/ratio);
					}
					int starty=(int)sps[i].coords[1]-maxdist;
					int endy=(int)sps[i].coords[1]+maxdist;
					if(pointy>=starty&&pointy<endy){
						double yr=Math.abs(sps[i].coords[1]-pointy);
						float multiplier=(float)getinterpgaus(yr);
						int startx=(int)sps[i].coords[0]-maxdist;
						int endx=(int)sps[i].coords[0]+maxdist;
						if(pointx>=startx&&pointx<endx){
							double xr=Math.abs(sps[i].coords[0]-pointx);
							float multiplierx=(float)getinterpgaus(xr);
							for(int j=0;j<nchannels;j++){
								float intensity=sps[i].get_brightness(sm,j)*multiplierx*multiplier*multiplierz;
								fintensity[j]+=intensity;
								sps[i].photonsused[j]+=intensity;
							}
						}else{
							// wrap around detection in the x direction
							if(startx<0&&pointx>=(startx+set.boxpixels)){
								double xr=Math.abs(sps[i].coords[0]-(pointx-set.boxpixels));
								float multiplierx=(float)getinterpgaus(xr);
								for(int j=0;j<nchannels;j++){
									float intensity=sps[i].get_brightness(sm,j)*multiplierx*multiplier*multiplierz;
									fintensity[j]+=intensity;
									sps[i].photonsused[j]+=intensity;
								}
							}else{
								if(endx>=set.boxpixels&&pointx<(endx-set.boxpixels)){
									double xr=Math.abs(sps[i].coords[0]-(pointx+set.boxpixels));
									float multiplierx=(float)getinterpgaus(xr);
									for(int j=0;j<nchannels;j++){
										float intensity=sps[i].get_brightness(sm,j)*multiplierx*multiplier*multiplierz;
										fintensity[j]+=intensity;
										sps[i].photonsused[j]+=intensity;

									}
								}
							}
						}
					}
				}
			}
		}
		for(int i=0;i<nchannels;i++){
			float tempback=back[i]*set.frametime2;
			fintensity[i]+=tempback;
		}
		if(set.noiseindex==2){
			return fintensity;
		}else{
			float[] photons=new float[nchannels];
			for(int i=0;i<nchannels;i++){
				photons[i]=set.random.poidev(fintensity[i]);
			}
			if(set.noiseindex==0){
				return photons;
			}else{
				float[] aint=new float[nchannels];
				for(int i=0;i<nchannels;i++){
					aint[i]=(float)set.random.analogdev((int)photons[i],set.gain,set.offset,set.readstdev);
				}
				return aint;
			}
		}
	}

	public float[] calc_point_tpe(sim_multispecies sm,float[] back,int pointx,int pointy){
		// here we calculate a point from the molecule data
		// multiple points are calculated for multiple channels
		int nchannels=sm.sslist[0].bright2.length;
		sim_particle[] sps=sm.sps;
		float[] fintensity=new float[nchannels];
		for(int i=0;i<sps.length;i++){
			boolean bleached=true;
			for(int j=0;j<nchannels;j++)
				if(sps[i].bleach_state[j]!=0){
					bleached=false;
					break;
				}
			if(!bleached){
				double zr=0.0;
				if(set.confineindex!=1&&set.confineindex!=3){
					zr=Math.abs(sps[i].coords[2]-set.fboxpixels/2.0f);
				}
				if(zr<=maxdistz){
					int starty=(int)sps[i].coords[1]-maxdist;
					int endy=(int)sps[i].coords[1]+maxdist;
					if(pointy>=starty&&pointy<endy){
						double yr=Math.abs(sps[i].coords[1]-pointy);
						int startx=(int)sps[i].coords[0]-maxdist;
						int endx=(int)sps[i].coords[0]+maxdist;
						if(pointx>=startx&&pointx<endx){
							double xr=Math.abs(sps[i].coords[0]-pointx);
							for(int j=0;j<nchannels;j++){
								float intensity=sps[i].get_brightness(sm,j)*(float)get_tpe_psf(xr,yr,zr);
								fintensity[j]+=intensity;
								sps[i].photonsused[j]+=intensity;
							}
						}else{
							// wrap around detection in the x direction
							if(startx<0&&pointx>=(startx+set.boxpixels)){
								double xr=Math.abs(sps[i].coords[0]-(pointx-set.boxpixels));
								for(int j=0;j<nchannels;j++){
									float intensity=sps[i].get_brightness(sm,j)*(float)get_tpe_psf(xr,yr,zr);
									fintensity[j]+=intensity;
									sps[i].photonsused[j]+=intensity;
								}
							}else{
								if(endx>=set.boxpixels&&pointx<(endx-set.boxpixels)){
									double xr=Math.abs(sps[i].coords[0]-(pointx+set.boxpixels));
									for(int j=0;j<nchannels;j++){
										float intensity=sps[i].get_brightness(sm,j)*(float)get_tpe_psf(xr,yr,zr);
										fintensity[j]+=intensity;
										sps[i].photonsused[j]+=intensity;
									}
								}
							}
						}
					}
				}
			}
		}
		for(int i=0;i<nchannels;i++){
			float tempback=back[i]*set.frametime2;
			fintensity[i]+=tempback;
		}
		if(set.noiseindex==2){
			return fintensity;
		}else{
			float[] photons=new float[nchannels];
			for(int i=0;i<nchannels;i++){
				photons[i]=set.random.poidev(fintensity[i]);
			}
			if(set.noiseindex==0){
				return photons;
			}else{
				float[] aint=new float[nchannels];
				for(int i=0;i<nchannels;i++){
					aint[i]=(float)set.random.analogdev((int)photons[i],set.gain,set.offset,set.readstdev);
				}
				return aint;
			}
		}
	}

	public Image calc_demo_image(sim_multispecies sm,float maxintensity){
		// this method produces a demo image with molecules shown in 2D as dots
		// the molecules in channel 0 are depicted as green and in channel 1 are
		// depicted as red
		// a gaussian is drawn in the center scaled to maximum intensity
		// its color is related to the intensity in channel 0 and channel 1
		int nchannels=sm.sslist[0].bright2.length;
		sim_particle[] sps=sm.sps;
		int[][] intensity=new int[2][set.boxpixels*set.boxpixels];
		float[] signal=calc_point(sm,new float[nchannels]);
		int starty=(int)(set.boxpixels/2.0f)-maxdist;
		int endy=(int)(set.boxpixels/2.0f)+maxdist;
		for(int pointy=starty;pointy<endy;pointy++){
			double yr=Math.abs(set.boxpixels/2.0f-pointy);
			float multiplier=(float)getinterpgaus(yr);
			int startx=starty;
			int endx=endy;
			for(int pointx=startx;pointx<endx;pointx++){
				double xr=Math.abs(set.boxpixels/2.0f-pointx);
				float multiplierx=(float)getinterpgaus(xr);
				int gintensity=(int)(((signal[0]*multiplierx*multiplier)/maxintensity)*255.0f);
				// int
				// gintensity=(int)(((1.0f*multiplierx*multiplier)/maxintensity)*255.0f);
				if(gintensity>255){
					gintensity=255;
				}
				int rintensity=0;
				if(nchannels>1){
					rintensity=(int)(((signal[1]*multiplierx*multiplier)/maxintensity)*255.0f);
					if(rintensity>255){
						rintensity=255;
					}
				}
				intensity[0][pointx+set.boxpixels*pointy]=gintensity;
				intensity[1][pointx+set.boxpixels*pointy]=rintensity;
			}
		}
		for(int i=0;i<sps.length;i++){
			int xvalue=(int)sps[i].coords[0];
			int yvalue=(int)sps[i].coords[1];
			float b1=sps[i].get_brightness(sm,0);
			float b2=sps[i].get_brightness(sm,1);
			if(b1>0.0f){
				intensity[0][xvalue+yvalue*set.boxpixels]=255;
			}
			if(b2>0.0f){
				intensity[1][xvalue+yvalue*set.boxpixels]=255;
			}
		}
		// BufferedImage retimage=new
		// BufferedImage(set.boxpixels,set.boxpixels,BufferedImage.TYPE_INT_ARGB);
		int[] rgb=new int[set.boxpixels*set.boxpixels];
		for(int i=0;i<set.boxpixels*set.boxpixels;i++){
			rgb[i]=0xff000000+(intensity[1][i]<<16)+(intensity[0][i]<<8);
		}
		MemoryImageSource source=new MemoryImageSource(set.boxpixels,set.boxpixels,rgb,0,set.boxpixels);
		return Toolkit.getDefaultToolkit().createImage(source);
		// retimage.setRGB(0,0,set.boxpixels,set.boxpixels,rgb,0,1);
		// return retimage;
		// return rgb;
	}

	private double get_tpe_psf(double xr,double yr,double zr){
		double[] vals=get_tpe_w0_amp(2.0*stdev,zr);
		return vals[1]*getinterpgaus(xr,vals[0])*getinterpgaus(yr,vals[0]);
	}

	private double[] get_tpe_w0_amp(double w0,double lambda,double z){
		// note that lambda here must be in pixel units
		double[] retvals=new double[2];
		// 0.7071
		retvals[0]=0.7071*w0*Math.sqrt(1+lambda*lambda*z*z/(Math.PI*Math.PI*w0*w0*w0*w0));
		retvals[1]=0.4053*Math.pow(w0*0.7071/retvals[0],4.0);
		return retvals;
	}

	private double[] get_tpe_w0_amp(double w0,double z){
		return get_tpe_w0_amp(w0,0.8*set.fboxpixels/set.boxsize,z);
	}

	private double getinterpgaus(double r){
		int rp=(int)(r*10.0);
		double rem=r*10.0-rp;
		if(rp<999){
			return rem*(model_gaus[rp+1]-model_gaus[rp])+model_gaus[rp];
		}else{
			return 0.0;
		}
	}

	private double getinterpgaus(double r,double tempw0){
		return getinterpgaus(r*2.0*stdev/tempw0);
	}

	private void init_gaus(){
		model_gaus=new double[1000];
		for(int i=0;i<1000;i++){
			double r=0.1*i;
			model_gaus[i]=Math.exp(-(r*r)/(2.0*stdev*stdev));
		}
	}

}
