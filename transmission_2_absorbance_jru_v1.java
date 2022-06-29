/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;
import jguis.*;

public class transmission_2_absorbance_jru_v1 implements PlugIn {
	//this plugin converts an image from transmission into absorbance
	//0 is assumed to be no transmission (camera offset must be 0)
	//the roi is assumed to contain the white background (100% transmission)

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		Roi roi=imp.getRoi();
		boolean[] mask=jutils.roi2mask(roi,width,height);
		Object[] stack=jutils.stack2array(imp.getStack());
		float[] white=getSpectrum(stack,width,height,mask,"Avg");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Camera Offset",0,5,15,null);
		gd.addNumericField("Max OD",3.0f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float off=(float)gd.getNextNumber();
		float maxod=(float)gd.getNextNumber();
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<stack.length;i++){
			retstack.addSlice("",convert2OD(stack[i],off,white[i],maxod));
		}
		jutils.create_hyperstack("Converted OD",retstack,imp).show();
	}

	public float[] convert2OD(Object image,float offset,float white,float maxod){
		float[] temp=algutils.convert_arr_float2(image);
		float[] OD=new float[temp.length];
		for(int i=0;i<temp.length;i++){
			float T=(temp[i]-offset)/(white-offset);
			if(T<=0.0f) OD[i]=maxod;
			else OD[i]=(float)(-Math.log10(T));
		}
		return OD;
	}

	public float[] getSpectrum(Object[] stack,int width,int height,boolean[] mask,String stat){
		float[] spectrum=new float[stack.length];
		for(int i=0;i<stack.length;i++){
			spectrum[i]=jstatistics.getstatistic(stat,stack[i],width,height,mask,null);
		}
		return spectrum;
	}

}
