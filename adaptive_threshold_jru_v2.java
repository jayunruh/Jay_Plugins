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
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class adaptive_threshold_jru_v2 implements PlugIn, threshinterface {
	float multiplier;
	String method;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] methods=AutoThresholder.getMethods();
		gd.addChoice("Threshold_Method",methods,methods[0]);
		gd.addNumericField("Multiplier",1.0,5,15,null);
		gd.addNumericField("Thresh_Region_Size",25,0);
		gd.addNumericField("Thresh_Resolution",12,0);
		gd.addCheckbox("Output_Threshold",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		method=methods[gd.getNextChoiceIndex()];
		multiplier=(float)gd.getNextNumber();
		int regsize=(int)gd.getNextNumber();
		int skip=(int)gd.getNextNumber();
		boolean outthresh=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack mask=new ImageStack(width,height);
		ImageStack threshstack=new ImageStack(width,height);
		for(int i=0;i<size;i++){
			float[] temp=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
			float[] thresh=jthresh.getadaptivethresh(temp,width,height,this,regsize,skip);
			byte[] mask2=new byte[width*height];
			for(int j=0;j<width*height;j++) if(temp[j]>=thresh[j]) mask2[j]=(byte)255;
			mask.addSlice("",mask2);
			if(outthresh) threshstack.addSlice("",thresh);
		}
		new ImagePlus("Threshholded",mask).show();
		if(outthresh) new ImagePlus("Threshhold",threshstack).show();
	}

	public float calcthresh(float[] subset){
		float[] minmax=jstatistics.getminmax(subset);
		int[] hist=jstatistics.fhistogram(subset,new float[]{256.0f,minmax[0],minmax[1]});
		int index=(new AutoThresholder()).getThreshold(method,hist);
		float thresh=((float)index/255.0f)*(minmax[1]-minmax[0])+minmax[0];
		return multiplier*thresh;
	}

}
