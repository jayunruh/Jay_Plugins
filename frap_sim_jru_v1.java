/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import jalgs.jsim.*;
import ij.measure.*;
import jalgs.jfit.*;
import jalgs.jfft.*;

public class frap_sim_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we simulate FRAP recovery in a 2D plane
		//binding can occur from the cytoplasm
		//simulation starts with experimental distribution
		float carpetlengthsec=260.0f;
		int carpetwidth=96;
		float pixelsize=0.088f;
		float steptime=0.1f; //units are seconds
		int carpetlength=(int)(carpetlengthsec/steptime);
		float prefrapsec=260.0f;
		//float prefrapsec=0.0f;
		int prefrap=(int)(prefrapsec/steptime);
		float Dlow=0.053f; //values for wt
		float Dhigh=0.0061f;
		//float Dlow=0.053f; //values for S185D
		//float Dhigh=0.022f;
		//float Dlow=0.088f; //values for S185D ratio matching
		//float Dhigh=0.01f;
		float w0=0.2f;
		float z0=0.8f;
		//float Dlow=0.0075f;
		//float Dhigh=0.0075f;
		//float koff=0.004f;
		//float kon=25.0f*koff;
		float kon=0.0f;
		float koff=0.0f;
		float koff2=koff*steptime;
		float kon2=kon*steptime;
		boolean convolve=true;
		float[] initial=new float[carpetwidth*carpetwidth];
		float[] bound=new float[carpetwidth*carpetwidth];
		int halfcw=carpetwidth/2;
		float patchsize=0.5f; //patch halfsize
		float patchdist=7.0f;
		float patchdist2=13.0f;
		float patchangle1=(float)Math.toRadians(60.0f);
		float patchangle2=(float)Math.toRadians(30.0f);
		float capsize=15.0f; //actually stdev or half width of ifrap region
		int npatches=19;
		//int fudge=(int)((float)patchdist/(float)Math.sqrt(2.0));
		//int fudge=(int)patchdist;
		//float[][] patchcoords={{halfcw+fudge,halfcw},{halfcw,halfcw+fudge},{halfcw-fudge,halfcw},{halfcw,halfcw-fudge},{halfcw,halfcw}};
		//float[][] patchcoords={{halfcw-fudge,halfcw-fudge},{halfcw+fudge,halfcw+fudge},{halfcw-fudge,halfcw+fudge},{halfcw+fudge,halfcw-fudge},{halfcw,halfcw}};
		float[][] patchcoords=new float[npatches][2];
		for(int i=0;i<6;i++){
			patchcoords[i]=rotate_vector(patchdist,(float)i*patchangle1,halfcw,halfcw);
		}
		for(int i=0;i<12;i++){
			patchcoords[i+6]=rotate_vector(patchdist2,(float)i*patchangle2,halfcw,halfcw);
		}
		patchcoords[18]=new float[]{halfcw,halfcw};
		
		float[] Dprofile=getDprofile(carpetwidth,carpetwidth,patchcoords,patchsize,Dlow,Dhigh);
		ImageStack stack=new ImageStack(carpetwidth,carpetwidth);
		float[] carpet=new float[carpetwidth*(prefrap+carpetlength)];
		//create the initial prefrap distribution
		for(int i=0;i<carpetwidth;i++){
			for(int j=0;j<carpetwidth;j++){
				initial[j+i*carpetwidth]=1.0f; //uncomment this for equilibrium measurements
				for(int k=0;k<npatches;k++){
					float xdist2=(float)Math.abs((float)i-patchcoords[k][0]);
					float ydist2=(float)Math.abs((float)j-patchcoords[k][1]);
					double dist2=(double)(xdist2*xdist2+ydist2*ydist2);
					//float tempval=(float)Math.exp(-dist2/(2.0*(double)(patchsize*patchsize)));
					float tempval=2.0f*(float)patchsize;
					tempval=(float)Math.pow(tempval*tempval/((float)dist2+tempval*tempval),0.8);
					initial[j+i*carpetwidth]+=10.0f*tempval;
				}
			}
		}
		//Dprofile=(new manipulate_quads()).shiftxycenter(Dprofile,carpetwidth,carpetwidth);
		new ImagePlus("D_profile",new FloatProcessor(carpetwidth,carpetwidth,Dprofile,null)).show();
		//finite_differences fd=new finite_differences(0,10,pixelsize,steptime,Dprofile,carpetwidth,carpetwidth);
		gauss_convolution fd=new gauss_convolution(0,pixelsize,steptime,Dprofile,carpetwidth,carpetwidth);
		//new ImagePlus("Diffusion Profile1",new FloatProcessor((float[][])fd.profiles[0])).show();
		//new ImagePlus("Diffusion Profile2",new FloatProcessor((float[][])fd.profiles[1])).show();
		for(int i=0;i<prefrap;i++){
			Object image=fd.steptime(initial);
			initial=(float[])image;
			handle_binding(initial,bound,carpetwidth,kon2,koff2,patchcoords,patchsize);
			float[] output=add_images(initial,bound);
			stack.addSlice("",output);
			if(convolve){
				float[] convolved=convolve_psf(output,carpetwidth,w0/pixelsize,z0/pixelsize);
				System.arraycopy(convolved,0,carpet,i*carpetwidth,carpetwidth);
			}
			IJ.showProgress(i,prefrap);
			if(IJ.escapePressed()){break;}
		}
		//frap the sample
		for(int i=0;i<carpetwidth;i++){
			float ydist=(float)Math.abs(i-halfcw);
			for(int j=0;j<carpetwidth;j++){
				float xdist=(float)Math.abs(j-halfcw);
				//float xdist=(float)Math.abs(j-carpetwidth/5);
				//float xdist=(float)j;
				//if(xdist>halfcw) xdist=(float)(carpetwidth-j);
				float dist2=xdist*xdist+ydist*ydist;
				if(dist2>(capsize*capsize)){
					initial[j+i*carpetwidth]=0.0f;
					bound[j+i*carpetwidth]=0.0f;
				}
			}
		}
		stack.addSlice("",add_images(initial,bound));
		if(convolve){
			float[] convolved=convolve_psf(add_images(initial,bound),carpetwidth,w0/pixelsize,z0/pixelsize);
			System.arraycopy(convolved,0,carpet,prefrap*carpetwidth,carpetwidth);
		}
		//fd=new gauss_convolution(0,pixelsize,steptime,0.013f,carpetwidth,carpetwidth); //use only for uniform dispersion
		for(int i=1;i<carpetlength;i++){
			Object image=fd.steptime(initial);
			initial=(float[])image;
			handle_binding(initial,bound,carpetwidth,kon2,koff2,patchcoords,patchsize);
			float[] output=add_images(initial,bound);
			stack.addSlice("",output);
			if(convolve){
				float[] convolved=convolve_psf(output,carpetwidth,w0/pixelsize,z0/pixelsize);
				System.arraycopy(convolved,0,carpet,(i+prefrap)*carpetwidth,carpetwidth);
			}
			IJ.showProgress(i,carpetlength);
			if(IJ.escapePressed()){break;}
		}
		if(convolve){
			ImagePlus imp=new ImagePlus("Sim FRAP Carpet",new FloatProcessor(carpetwidth,carpetlength+prefrap,carpet,null));
			jutils.set_psize(imp,pixelsize);
			imp.show();
			//also create the central profile for fitting
			int analhalfwidth=(int)capsize;
			float[] profile=new float[carpetlength+prefrap];
			for(int i=0;i<(carpetlength+prefrap);i++){
				for(int j=(halfcw-analhalfwidth);j<(halfcw+analhalfwidth);j++){
					profile[i]+=carpet[i*carpetwidth+j];
				}
				/*for(int j=(carpetwidth/5-analhalfwidth);j<(carpetwidth/5+analhalfwidth);j++){
					profile[i]+=carpet[i*carpetwidth+j];
				}*/
				//for(int j=0;j<analhalfwidth;j++) profile[i]+=carpet[i*carpetwidth+j];
				//for(int j=(carpetwidth-analhalfwidth);j<carpetwidth;j++) profile[i]+=carpet[i*carpetwidth+j];
				profile[i]/=(float)(2*analhalfwidth);
			}
			float ptratio=carpet[(prefrap-1)*carpetwidth+halfcw]/carpet[(prefrap-1)*carpetwidth+halfcw+3];
			IJ.log("P T Ratio = "+ptratio);
			int timebin=20;
			int newlength=(int)((float)profile.length/(float)timebin);
			float[] binnedprofile=new float[newlength];
			float[] xvals=new float[newlength];
			for(int i=0;i<newlength;i++){
				for(int j=0;j<timebin;j++){
					binnedprofile[i]+=profile[j+i*timebin];
				}
				binnedprofile[i]/=(float)timebin;
				xvals[i]=steptime*(float)timebin*(float)i;
			}
			new PlotWindow4("Binned Temporal Profile","time (s)","Intensity",xvals,binnedprofile).draw();
			IJ.run("batch FRAP fit jru v1", "frames_before=130 frames_to=260 minimum=2.00000 maximum=768.00000 min_tau_ratio=2.00000");
		} else {
			ImagePlus imp=new ImagePlus("Sim FRAP",stack);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(1,1,carpetlength+prefrap);
			jutils.set_psize(imp,pixelsize);
			imp.getCalibration().frameInterval=steptime;
			imp.show();
		}
	}

	public float[] getDprofile(int width,int height,float[][] patchcoords,float patchsize,float Dlow,float Dhigh){
		float[] Dprofile=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				Dprofile[j+i*width]=Dlow;
				for(int k=0;k<patchcoords.length;k++){
					float xdist2=(float)Math.abs((float)i-patchcoords[k][0]);
					float ydist2=(float)Math.abs((float)j-patchcoords[k][1]);
					double dist2=(double)(xdist2*xdist2+ydist2*ydist2);
					float tempval=2.0f*(float)patchsize;
					tempval=(float)Math.pow(tempval*tempval/((float)dist2+tempval*tempval),0.8);
					Dprofile[j+i*width]+=(Dhigh-Dlow)*tempval;
					if(Dprofile[j+i*width]<Dhigh) Dprofile[j+i*width]=Dhigh;
					if(Dprofile[j+i*width]>Dlow) Dprofile[j+i*width]=Dlow;
				}
			}
		}
		//make the D matrix so it contains 32 unique values or less
		float nlevels=32.0f;
		for(int i=0;i<width*height;i++){
			int scaled=(int)(nlevels*(Dprofile[i]-Dhigh)/(Dlow-Dhigh));
			float scaled2=((float)scaled)/nlevels;
			Dprofile[i]=Dhigh+scaled2*(Dlow-Dhigh);
		}
		return Dprofile;
	}

	public float[] rotate_vector(float r,float theta,float xcenter,float ycenter){
		return new float[]{xcenter+r*(float)Math.cos(theta),ycenter+r*(float)Math.sin(theta)};
	}

	public float[] add_images(float[] img1,float[] img2){
		float[] temp=new float[img1.length];
		for(int i=0;i<img1.length;i++){
			temp[i]=img1[i]+img2[i];
		}
		return temp;
	}

	public void handle_binding(float[] image,float[] bound,int width,float kon2,float koff2,float[][] patchcoords,float patchsize){
		if(kon2==0.0f && koff2==0.0f) return;
		float[] temp=image.clone();
		float[] tempb=bound.clone();
		for(int i=0;i<patchcoords.length;i++){
			int yc=(int)patchcoords[i][1];
			int xc=(int)patchcoords[i][0];
			for(int j=(yc-(int)patchsize);j<=(yc+(int)patchsize);j++){
				for(int k=(xc-(int)patchsize);k<=(xc+(int)patchsize);k++){
					float change=kon2*temp[k+j*width]-koff2*tempb[k+j*width];
					bound[k+j*width]+=change;
					image[k+j*width]-=change;
				}
			}
		}
	}

	public float[] convolve_psf(float[] image,int maxsize,float w02,float z02){
		int xhalfsize=(int)(w02*1.5f);
		int zhalfsize=(int)(z02*1.5f);
		gausfunc gf=new gausfunc();
		float[] zprofile=gf.get_norm_func(-zhalfsize,2*zhalfsize,1.0,z02*0.5f);
		float[] xprofile=gf.get_norm_func(-xhalfsize,2*xhalfsize,1.0,w02*0.5f);
		float[] conv=new float[maxsize];
		int center=maxsize/2;
		for(int i=0;i<maxsize;i++){ //carpet coordinates
			for(int j=i-xhalfsize;j<(i+xhalfsize);j++){ //profile coordinates
				int xpos=j;
				if(xpos<0) xpos+=maxsize;
				if(xpos>=maxsize) xpos-=maxsize;
				for(int k=center-zhalfsize;k<(center+zhalfsize);k++){
					conv[i]+=image[xpos+k*maxsize]*xprofile[j-i+xhalfsize]*zprofile[k-center+zhalfsize];
				}
			}
		}
		return conv;
	}

}
