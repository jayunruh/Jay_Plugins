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
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import ij.measure.*;
import java.io.*;
import jguis.*;
import jalgs.jfit.*;

public class FRET_heat_map_jru_v2 implements PlugIn {
	//this version does a better job of calculating correct donor and acceptor intensities
	//the acceptor intensity is corrected for bleedthrough using the after images assuming all acceptor is bleached
	//the donor intensity is corrected for FRET by using the after intensity again assuming all acceptor is bleached
	//the acceptor intensity still contains the contribution due to FRET which requires the direct excitation efficiency to calculate

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Start_Before Slice",1,0);
		gd.addNumericField("End_Before Slice",4,0);
		gd.addNumericField("Start_After Slice",5,0);
		gd.addNumericField("End_After Slice",8,0);
		gd.addNumericField("Donor_Channel",1,0);
		gd.addNumericField("Acceptor_Channel",2,0);
		gd.addCheckbox("Use_Trend",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int bfstart=(int)gd.getNextNumber();
		int bfend=(int)gd.getNextNumber();
		int afstart=(int)gd.getNextNumber();
		int afend=(int)gd.getNextNumber();
		int dchannel=(int)gd.getNextNumber();
		int achannel=(int)gd.getNextNumber();
		boolean trend=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi backroi=imp.getRoi();
		ImageStack stack=imp.getStack();
		int width=imp.getWidth();
		int height=imp.getHeight();
		float background=0.0f;
		float abackground=0.0f;
		if(backroi!=null){
			background=get_background(stack.getPixels(dchannel),width,height,backroi);
			abackground=get_background(stack.getPixels(achannel),width,height,backroi);
		}
		int nchannels=imp.getNChannels();
		int nslices=imp.getNFrames();
		int bfslices=bfend-bfstart+1;
		int afslices=afend-afstart+1;
		float[] bfxvals=new float[bfslices];
		for(int i=0;i<bfslices;i++) bfxvals[i]=(float)i;
		float[] afxvals=new float[afslices];
		for(int i=0;i<afslices;i++) afxvals[i]=(float)i;
		float[] before=new float[width*height];
		float[] after=new float[width*height];
		float[] acc=new float[width*height];
		float[] accafter=new float[width*height];
		float[] sum=new float[width*height];
		for(int j=0;j<width*height;j++){
			float[] temp=new float[bfslices];
			for(int i=(bfstart-1);i<bfend;i++){
				Object image=stack.getPixels(i*nchannels+dchannel);
				if(image instanceof float[]){
					temp[i-bfstart+1]=((float[])image)[j];
				} else {
					temp[i-bfstart+1]=(float)(((short[])image)[j]&0xffff);
				}
			}
			before[j]=trend_data(temp,bfxvals,true,!trend);
		}
		for(int j=0;j<width*height;j++){
			float[] temp=new float[bfslices];
			for(int i=(bfstart-1);i<bfend;i++){
				Object image=stack.getPixels(i*nchannels+achannel);
				if(image instanceof float[]){
					temp[i-bfstart+1]=((float[])image)[j];
				} else {
					temp[i-bfstart+1]=(float)(((short[])image)[j]&0xffff);
				}
			}
			acc[j]=trend_data(temp,bfxvals,true,!trend);
			temp=new float[afslices];
			for(int i=(afstart-1);i<afend;i++){
				Object image=stack.getPixels(i*nchannels+achannel);
				if(image instanceof float[]){
					temp[i-afstart+1]=((float[])image)[j];
				} else {
					temp[i-afstart+1]=(float)(((short[])image)[j]&0xffff);
				}
			}
			accafter[j]=trend_data(temp,afxvals,false,!trend);
		}
		for(int j=0;j<width*height;j++){
			float[] temp=new float[afslices];
			for(int i=(afstart-1);i<afend;i++){
				Object image=stack.getPixels(i*nchannels+dchannel);
				if(image instanceof float[]){
					temp[i-afstart+1]=((float[])image)[j];
				} else {
					temp[i-afstart+1]=(float)(((short[])image)[j]&0xffff);
				}
				after[j]=trend_data(temp,afxvals,false,!trend);
			}
		}
		float taccafter=0.0f;
		float tdafter=0.0f;
		for(int i=0;i<width*height;i++){
			before[i]-=background;
			after[i]-=background;
			acc[i]-=abackground;
			accafter[i]-=abackground;
			taccafter+=accafter[i];
			tdafter+=after[i];
		}
		float fbleed=taccafter/tdafter;
		IJ.log("fraction bleedthrough = "+fbleed);
		for(int i=0;i<width*height;i++){
			acc[i]-=before[i]*fbleed;
		}
		ImagePlus imp1=new ImagePlus("Before",new FloatProcessor(width,height,before,null));
		imp1.copyScale(imp);
		//imp1.show();
		ImagePlus imp2=new ImagePlus("After",new FloatProcessor(width,height,after,null));
		imp2.copyScale(imp);
		//imp2.show();
		ImagePlus imp3=new ImagePlus("Acceptor",new FloatProcessor(width,height,acc,null));
		
		final Hist2DWindow cw=new Hist2DWindow();
		cw.init(imp1,imp2,imp3,imp2,2);
		Hist2DWindow.launch_frame("Interactive 2D Histogram",cw);
	}

	float get_background(Object image,int width,int height,Roi backroi){
		float sum=0.0f;
		int area=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(backroi.contains(j,i)){
					if(image instanceof float[]){
						sum+=((float[])image)[j+i*width];
					} else {
						sum+=(float)(((short[])image)[j+i*width]&0xffff);
					}
					area++;
				}
			}
		}
		return sum/(float)area;
	}

	float trend_data(float[] data,float[] xvals,boolean bf,boolean avg){
		if(avg){
			float sum=0.0f;
			for(int i=0;i<data.length;i++){
				sum+=data[i];
			}
			return sum/(float)data.length;
		} else {
			float[] coef=(new linleastsquares()).get_amp_offset(xvals,data,true);
			if(bf){
				return coef[1]+coef[0]*xvals[xvals.length-1];
			} else {
				return coef[1]+coef[0]*xvals[0];
			}
		}
	}

}
