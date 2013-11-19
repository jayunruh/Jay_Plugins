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
import jalgs.jfit.*;
import jguis.*;
import jalgs.*;

public class pca_spectral_unmixing_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_Components",2,0);
		gd.addCheckbox("Output_Variances",false);
		gd.addCheckbox("Output_All_Components",false);
		gd.addCheckbox("Get_Spectra_from_Data",true);
		gd.addNumericField("#Data_Points_for_Spectra",5,0);
		gd.addCheckbox("Unmix",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int nspectra=(int)gd.getNextNumber();
		boolean outvar=gd.getNextBoolean();
		boolean outcomp=gd.getNextBoolean();
		boolean get_from_data=gd.getNextBoolean();
		int navg=(int)gd.getNextNumber();
		boolean unmix=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int nchannels=stack.getSize();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int npixels=width*height;
		float[][] y=new float[nchannels][npixels]; //this is the data stack
		for(int i=0;i<nchannels;i++){
			y[i]=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
		}
		float[] avgint=new float[nchannels];
		float[] avgy=new float[npixels];
		for(int i=0;i<npixels;i++){
			for(int j=0;j<nchannels;j++){
				avgint[j]+=y[j][i];
				avgy[i]+=y[j][i];
			}
		}
		for(int i=0;i<npixels;i++) avgy[i]/=(float)nchannels;
		float[][] corry=new float[npixels][nchannels];
		for(int i=0;i<nchannels;i++){
			avgint[i]/=(float)npixels;
			for(int j=0;j<npixels;j++){
				corry[j][i]=y[i][j]-avgint[i];
				//corry[j][i]=y[i][j];
			}
		}
		float[] w=new float[nchannels];
		float[][] v=new float[nchannels][nchannels];
		float[][] u=jsvd.svdcmp(corry,w,v);
		v=transpose_neg(v);
		if(outvar){new PlotWindow4("Variances","x","y",w).draw();}
		if(outcomp){new PlotWindow4("Components","x","y",v,null).draw();}
		int[] order=jsort.get_javasort_order(w);
		float[][] spectra=new float[nspectra][nchannels];
		if(!get_from_data){
			for(int i=0;i<nspectra;i++){
				float[] temp=v[order[nchannels-i-1]].clone();
				//for(int j=0;j<nchannels;j++) temp[j]+=avgint[j];
				spectra[i]=normalize2(temp);
			}
		} else {
			RoiManager rman=RoiManager.getInstance();
			if(rman==null)
				rman=new RoiManager();
			rman.runCommand("show all");
			for(int i=0;i<nspectra;i++){
				float[] dist=distance_matrix(v[order[nchannels-i-1]],corry,avgy);
				int[] minindices=getnmin(dist,navg);
				int besty=(int)((float)minindices[0]/(float)width);
				int bestx=minindices[0]-besty*width;
				IJ.log("best point "+i+" , "+minindices[0]);
				rman.addRoi(new PointRoi(bestx,besty));
				for(int k=0;k<navg;k++){
					float[] temp=new float[nchannels];
					for(int j=0;j<nchannels;j++) temp[j]=y[j][minindices[k]];
					temp=normalize2(temp);
					for(int j=0;j<nchannels;j++) spectra[i][j]+=temp[j]/(float)navg;
				}
			}
		}
		new PlotWindow4("Spectra","Channel","Intensity",spectra,null).draw();
		if(unmix){
			WindowManager.setCurrentWindow(imp.getWindow());
			StringBuffer args=new StringBuffer("what_would_you_like_to_do?=Unmix use_plot_spectra spectra=Spectra how_many_species_are_present="+nspectra+" truncate_negative_values species_1=series0");
			for(int i=1;i<nspectra;i++){
				args.append(" species_"+(i+1)+"=series"+i);
			}
			IJ.run("linear unmixing jru v1",args.toString());
		}
		/*if(u!=null){
			IJ.log("w");
			IJ.log(table_tools.print_float_array(w));
			IJ.log("v");
			IJ.log(table_tools.print_float_array(v));
			IJ.log("u");
			IJ.log(table_tools.print_float_array(u));
		}*/
	}

	public float[][] transpose_neg(float[][] input){
		float[][] output=new float[input[0].length][input.length];
		for(int i=0;i<input.length;i++){
			for(int j=0;j<input[0].length;j++){
				output[j][i]=-input[i][j];
			}
		}
		return output;
	}

	public int[] getnmin(float[] vec,int nmin){
		int[] temp=new int[nmin];
		boolean[] taken=new boolean[vec.length];
		for(int i=0;i<nmin;i++){
			boolean first=true;
			float min=0.0f; int minindex=0;
			for(int j=0;j<vec.length;j++){
				if(vec[j]>=0.0f && (!taken[j])){
					if(first){
						min=vec[j]; minindex=j; first=false;
					} else {
						if(vec[j]<min){min=vec[j]; minindex=j;}
					}
				}
			}
			taken[minindex]=true;
			temp[i]=minindex;
		}
		return temp;
	}

	public float[] distance_matrix(float[] target,float[][] ycorr,float[] yavg){
		float[] dist=new float[ycorr.length];
		for(int i=0;i<ycorr.length;i++){
			if(yavg[i]>0.0f){
				dist[i]=norm_dist(target,ycorr[i]);
			} else {
				dist[i]=-1.0f;
			}
		}
		return dist;
	}

	public float norm_dist(float[] vec1,float[] vec2){
		float[] normvec1=normalize(vec1);
		float[] normvec2=normalize(vec2);
		float dist=0.0f;
		for(int i=0;i<normvec1.length;i++){
			dist+=(normvec1[i]-normvec2[i])*(normvec1[i]-normvec2[i]);
		}
		return (float)Math.sqrt(dist);
	}

	public float[] normalize(float[] vec1){
		float[] vec=vec1.clone();
		float sum=0.0f;
		for(int i=0;i<vec.length;i++) sum+=(vec[i]*vec[i]);
		sum=(float)Math.sqrt(sum);
		for(int i=0;i<vec.length;i++) vec[i]/=sum;
		return vec;
	}

	public float[] normalize2(float[] vec1){
		float[] vec=vec1.clone();
		float sum=0.0f;
		for(int i=0;i<vec.length;i++) sum+=vec[i];
		for(int i=0;i<vec.length;i++) vec[i]/=sum;
		return vec;
	}

}
