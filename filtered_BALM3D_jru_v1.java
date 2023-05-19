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
import ij.plugin.frame.*;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import jalgs.jseg.*;

public class filtered_BALM3D_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we implement a 3D filtered BALM method (no fitting)
		//based on Munck et al Nature Methods 2012
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int pframes=imp.getNFrames();
		int nframes=pframes;
		int nslices=imp.getNSlices();
		int nchannels=imp.getNChannels();
		GenericDialog gd=new GenericDialog("Options");
		boolean blurtime=true;
		gd.addCheckbox("Blur_time",true);
		float timesigma=0.8f;
		gd.addNumericField("Time_blur_sigma(frames)",timesigma,5,15,null);
		float filtersigma1=0.08f;
		gd.addNumericField("Spatial_sigma(um)",filtersigma1,5,15,null);
		float filtersigmaz1=0.24f;
		gd.addNumericField("Z_sigma(um)",filtersigmaz1,5,15,null);
		float filterratio=2.0f;
		gd.addNumericField("Mexican_hat_ratio",filterratio,5,15,null);
		float zratio=(float)jutils.get_pdepth(imp)/(float)jutils.get_psize(imp);
		gd.addNumericField("Z_ratio",zratio,5,15,null);
		gd.addNumericField("#_of_frames",pframes,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		blurtime=gd.getNextBoolean();
		timesigma=(float)gd.getNextNumber();
		filtersigma1=(float)gd.getNextNumber();
		filtersigmaz1=(float)gd.getNextNumber();
		filterratio=(float)gd.getNextNumber();
		zratio=(float)gd.getNextNumber();
		pframes=(int)gd.getNextNumber();
		
		float filtersigma=filtersigma1/(float)jutils.get_psize(imp);
		float zsize=zratio*(float)jutils.get_psize(imp);
		float filtersigmaz=filtersigmaz1/zsize;
		//start by filtering the data in the time domain with a gaussian stdev of 0.8
		ImageStack stack=imp.getStack();
		Object[][][] newslices=new Object[nchannels][pframes][nslices];
		for(int i=0;i<pframes;i++){
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchannels;k++){
					newslices[k][i][j]=new float[width*height];
				}
			}
		}
		gausfunc gf=new gausfunc();
		float[] symblurfunc=jsmooth.generate_symblurfunc(timesigma,gf);
		for(int i=0;i<width*height;i++){
			for(int j=0;j<nchannels;j++){
				for(int k=0;k<nslices;k++){
					Object[] tseries=jutils.get3DTSeries(stack,k,j,pframes,nslices,nchannels);
					float[] temp=algutils.convert_arr_float(algutils.get_stack_col(tseries,width,height,i,pframes));
					if(blurtime) jsmooth.convsym1D(temp,symblurfunc);
					set_t_column(newslices,temp,j,k,pframes,i,width,height);
				}
			}
		}
		//now calculate the average of all images
		float[][] avgprofile=new float[nchannels][pframes];
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<pframes;j++){
				for(int k=0;k<nslices;k++){
					avgprofile[i][j]+=jstatistics.getstatistic("Avg",newslices[i][j][k],null);
				}
				avgprofile[i][j]/=nslices;
			}
		}
		//and calculate the sum projection for each channel and slice
		float[][][] sumproj=new float[nchannels][nslices][];
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<nslices;j++){
				Object[] tseries=get_t_series(newslices,i,j,pframes);
				sumproj[i][j]=algutils.get_stack_proj_stat("Sum",tseries,width,height,pframes,null);
			}
		}
		//now calculate the differential image in place
		Object[][] recon=new Object[nchannels][nslices];
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<nslices;j++) recon[i][j]=new float[width*height];
		}
		Object[][][] diff=new Object[nchannels][pframes-1][nslices];
		for(int c=0;c<nchannels;c++){
			for(int i=0;i<(pframes-1);i++){
				float alpha=avgprofile[c][i+1]/avgprofile[c][i];
				for(int j=0;j<nslices;j++){
					//float[] tdiff=new float[width*height];
					float[] curr=(float[])newslices[c][i][j];
					float[] next=(float[])newslices[c][i+1][j];
					for(int k=0;k<width*height;k++) curr[k]=(float)Math.abs(alpha*curr[k]-next[k]);
					diff[c][i][j]=curr;
				}
				//now the mexican hat
				jsmooth.mexhat3D(diff[c][i],filtersigma,filtersigmaz,filterratio,width,height,gf);
				//now truncate at zero and add to the reconstruction
				for(int j=0;j<nslices;j++){
					float[] tdiff=(float[])diff[c][i][j];
					float[] trecon=(float[])recon[c][j];
					for(int k=0;k<width*height;k++){
						if(tdiff[k]<0.0f) tdiff[k]=0.0f;
						trecon[k]+=tdiff[k];
					}
				}
				IJ.showProgress(i,pframes-1);
			}
			for(int i=0;i<nslices;i++){
				float[] trecon=(float[])recon[c][i];
				for(int j=0;j<width*height;j++) trecon[j]*=(float)Math.sqrt(sumproj[c][i][j]);
			}
		}
		ImageStack outstack=new ImageStack(width,height);
		for(int i=0;i<nslices;i++){
			for(int j=0;j<nchannels;j++){
				outstack.addSlice("",recon[j][i]);
				outstack.addSlice("",sumproj[j][i]);
			}
		}
		ImageStack Dstack=new ImageStack(width,height);
		for(int i=0;i<(pframes-1);i++){
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchannels;k++){
					Dstack.addSlice("",diff[k][i][j]);
				}
			}
		}
		ImagePlus pimpimp=jutils.create_hyperstack("filtered BALM reconstruction",outstack,1,nslices,nchannels*2,true,null);
		pimpimp.copyScale(imp);
		pimpimp.show();
		ImagePlus dimp=jutils.create_hyperstack("Filtered Differential Stack",Dstack,pframes-1,nslices,nchannels,true,null);
		dimp.copyScale(imp);
		dimp.show();
	}

	public void set_t_column(Object[][][] stack,float[] source,int chan,int slice,int frames,int index,int width,int height){
		Object[] tseries=get_t_series(stack,chan,slice,frames);
		algutils.set_stack_col(tseries,width,height,index,frames,source);
	}

	public Object[] get_t_series(Object[][][] stack,int chan,int slice,int frames){
		Object[] substack=new Object[frames];
		for(int i=0;i<frames;i++) substack[i]=stack[chan][i][slice];
		return substack;
	}

}
