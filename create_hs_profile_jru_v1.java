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
import ij.text.*;
import ij.util.*;
import java.io.*;
import jalgs.*;
import jguis.*;

public class create_hs_profile_jru_v1 implements PlugIn {
	//this plugin calculates some statistic for each image in a stack (or region of interest) and plots it
	int index,projindex,histbins;
	float histstart,histend;

	public void run(String arg) {
		String[] options=jstatistics.stats;
		GenericDialog gd = new GenericDialog("Options");
		index=1;
		gd.addChoice("Profile_Statistic?",options,options[index]);
		projindex=1;
		gd.addChoice("Projection_Statistic?",options,options[projindex]);
		gd.addCheckbox("Z_Profile_(Time_projection)?",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		index=gd.getNextChoiceIndex();
		projindex=gd.getNextChoiceIndex();
		boolean zprof=gd.getNextBoolean();
		String statistic=options[index];
		String projstat=options[projindex];
		//get the image and its info
		ImagePlus imp = WindowManager.getCurrentImage();
		int height = imp.getHeight();
		int width = imp.getWidth();
		Rectangle r=null;
		boolean[] mask=null;
		Roi roi = imp.getRoi();
		if(roi!=null){
			if(roi.getType()==0){
				r=roi.getBounds();
			} else {
				mask=jutils.roi2mask(roi,width,height);
			}
		}
		ImageStack stack = imp.getStack();
		int slices = imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int pts=frames;
		if(zprof){pts=slices;}
		float[][] spectral_data = new float[channels][pts];
		jstatistics js=new jstatistics();
		float[] histoptions=null;
		/*if(statistic=="Mode" || projstat=="Mode"){
			histoptions=gethistoptions();
			if(histoptions==null){return;}
		}*/
		if(zprof){
			for(int i=0;i<channels;i++){
				for(int j=0;j<slices;j++){
					Object proj=null;
					if(frames==1) proj=jutils.get3DSlice(stack,0,j,i,frames,slices,channels);
					else proj=jutils.get3DProjTStat(stack,j,i,frames,slices,channels,projstat);
					if(mask==null) spectral_data[i][j]=js.getstatistic(statistic,proj,width,height,r,null);
					else spectral_data[i][j]=js.getstatistic(statistic,proj,width,height,mask,null);
				}
			}
		} else {
			for(int i=0;i<channels;i++){
				for(int j=0;j<frames;j++){
					Object proj=null;
					if(slices==1) proj=jutils.get3DSlice(stack,j,0,i,frames,slices,channels);
					else proj=jutils.get3DProjZStat(stack,j,i,frames,slices,channels,projstat);
					if(mask==null) spectral_data[i][j]=js.getstatistic(statistic,proj,width,height,r,null);
					else spectral_data[i][j]=js.getstatistic(statistic,proj,width,height,mask,null);
				}
			}
		}
		if(slices==1 && zprof){
			StringBuffer sb=new StringBuffer();
			sb.append(statistic+" ");
			sb.append("ch1="+spectral_data[0][0]);
			for(int i=1;i<channels;i++){
				sb.append(", ch"+(i+1)+"="+spectral_data[i][0]);
			}
			IJ.log(sb.toString());
		}
		else {
			PlotWindow4 plot = new PlotWindow4(statistic+" Profile","slice","intensity",spectral_data,null);
			plot.draw();
		}
	}

	/*private float[] gethistoptions(){
		GenericDialog gd=new GenericDialog("Histogram Options");
		gd.addNumericField("Histogram Bins",histbins,0);
		gd.addNumericField("Histogram Start",histstart,5,10,null);
		gd.addNumericField("Histogram End",histend,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		float[] temp=new float[3];
		histbins=(int)gd.getNextNumber();
		histstart=(float)gd.getNextNumber();
		histend=(float)gd.getNextNumber();
		temp[0]=(float)histbins; temp[1]=histstart; temp[2]=histend;
		return temp;
	}*/

}
