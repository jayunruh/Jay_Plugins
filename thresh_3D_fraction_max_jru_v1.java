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
import jguis.*;

public class thresh_3D_fraction_max_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we threshhold multiple z stack, each time calculating the threshhold as
		//the percentage of maximum for the stack
		//the maximum is calculated independently for each channel
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Normalization_Statistic",jstatistics.stats,jstatistics.stats[2]);
		gd2.addCheckbox("New_Image",true);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		String stat=jstatistics.stats[gd2.getNextChoiceIndex()];
		boolean newimg=gd2.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageProcessor ip=imp.getProcessor();
		double minthresh=ip.getMinThreshold();
		if(minthresh==ImageProcessor.NO_THRESHOLD){
			if(!IJ.isMacro()) IJ.showMessage("Threshold is not set");
			//return;
		}
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int currchan=imp.getChannel()-1;
		int currframe=imp.getFrame()-1;
		float[][] maxval=new float[channels][frames];
		for(int i=0;i<frames;i++){
			for(int j=0;j<channels;j++){
				float[] zprofile=new float[slices];
				for(int k=0;k<slices;k++){
					zprofile[k]=jutils.get3DSliceStat(stack,i,k,j,frames,slices,channels,stat);
				}
				maxval[j][i]=jstatistics.getstatistic(stat,zprofile,null);
			}
		}
		float fraction=0.1f;
		if(minthresh!=ImageProcessor.NO_THRESHOLD){
			fraction=(float)minthresh/maxval[currchan][currframe];
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Fraction of 3D "+stat+"?",fraction,5,15,null);
		gd.addCheckbox("Use Upper Threshold?",false);
		gd.addNumericField("Upper: Fraction of 3D "+stat+"?",fraction,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		fraction=(float)gd.getNextNumber();
		boolean upper=gd.getNextBoolean();
		float ufraction=(float)gd.getNextNumber();
		ImageStack threshstack=new ImageStack(width,height);
		int counter=0;
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					float thresh=maxval[k][i]*fraction;
					float uthresh=maxval[k][i]*ufraction;
					Object threshed=threshhold_image(jutils.get3DSlice(stack,i,j,k,frames,slices,channels),thresh,upper,uthresh);
					if(newimg) threshstack.addSlice("",threshed);
					else jutils.set3DSlice(stack,threshed,i,j,k,frames,slices,channels);
					counter++;
					IJ.showProgress(counter,frames*slices*channels);
				}
			}
		}
		if(newimg){
			ImagePlus imp2=new ImagePlus("Mask",threshstack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(channels,slices,frames);
			imp2.show();
		} else {
			imp.setStack(stack);
			imp.updateAndDraw();
		}
	}

	Object threshhold_image(Object pixels,float threshhold,boolean upper,float uthreshold){
		if(pixels instanceof float[]){
			float[] temppix=(float[])pixels;
			int length=temppix.length;
			float[] retvals=new float[length];
			for(int i=0;i<length;i++){
				if(temppix[i]>threshhold){
					retvals[i]=1.0f;
					if(upper && temppix[i]>=uthreshold) retvals[i]=0.0f;
				}
			}
			return retvals;
		} else if(pixels instanceof short[]) {
			short[] temppix=(short[])pixels;
			int length=temppix.length;
			short[] retvals=new short[length];
			for(int i=0;i<length;i++){
				float temp=(float)(temppix[i]&0xffff);
				if(temp>threshhold){
					retvals[i]=(short)1;
					if(upper && temppix[i]>=uthreshold) retvals[i]=(short)0;
				}
			}
			return retvals;
		} else {
			byte[] temppix=(byte[])pixels;
			int length=temppix.length;
			byte[] retvals=new byte[length];
			for(int i=0;i<length;i++){
				float temp=(float)(temppix[i]&0xff);
				if(temp>threshhold){
					retvals[i]=(byte)255;
					if(upper && temppix[i]>=uthreshold) retvals[i]=(byte)0;
				}
			}
			return retvals;
		}
	}
}
