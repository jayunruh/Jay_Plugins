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
import jalgs.jseg.*;
import jguis.*;

public class divide_by_blurred_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		float blurstdev=10.0f;
		float mindiv=1.0f;
		float multiplier=1.0f;
		boolean convfloat=true;
		boolean makenew=true;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Blur_Stdev",blurstdev,5,15,null);
		gd.addNumericField("Minimum_Division_Value",mindiv,5,15,null);
		gd.addNumericField("Multiplier",multiplier,5,15,null);
		gd.addCheckbox("Create_New",makenew);
		gd.addCheckbox("Convert_to_float",convfloat);
		gd.showDialog(); if(gd.wasCanceled()) return;
		blurstdev=(float)gd.getNextNumber();
		mindiv=(float)gd.getNextNumber();
		multiplier=(float)gd.getNextNumber();
		makenew=gd.getNextBoolean();
		convfloat=gd.getNextBoolean();
		Object[] stack2=jutils.stack2array(stack);
		boolean isfloat=(stack2[0] instanceof float[]);
		//if we are converting something to float we must make a new image
		if(!isfloat){
			if(convfloat){
				makenew=true;
			}
		}
		ImageStack stack3=new ImageStack(width,height);
		for(int i=0;i<stack2.length;i++){
			if(!makenew){
				if(isfloat){
					div_by_smooth((float[])stack2[i],width,height,blurstdev,mindiv,multiplier);
				} else {
					float[] temp=algutils.convert_arr_float(stack2[i]); //this copies the array
					div_by_smooth(temp,width,height,blurstdev,mindiv,multiplier);
					Object temp1=algutils.convert_array(temp,algutils.get_array_type(stack2[i]));
					temp=null;
					System.arraycopy(temp1,0,stack2[i],0,width*height);
					temp1=null;
				}
			} else {
				float[] temp=algutils.convert_arr_float(stack2[i]); //this copies the array
				div_by_smooth(temp,width,height,blurstdev,mindiv,multiplier);
				if(convfloat || isfloat){
					stack3.addSlice("",temp);
				} else {
					Object temp1=algutils.convert_array(temp,algutils.get_array_type(stack2[i]));
					temp=null;
					stack3.addSlice("",temp1);
				}
			}
			IJ.showProgress(i,stack2.length);
		}
		if(makenew){
			ImagePlus imp2=jutils.create_hyperstack("Divided_by_Blurred",stack3,imp);
			imp2.show();
		} else {
			imp.updateAndDraw();
		}
	}

	public static void div_by_smooth(float[] image,int width,int height,float blurstdev,float mindiv,float multiplier){
		//float[] smoothed=image.clone();
		//jsmooth.blur2D(smoothed,blurstdev,width,height);
		float[] smoothed=jutils.gaussian_blur(image,blurstdev,width,height);
		for(int i=0;i<width*height;i++){
			float temp=smoothed[i];
			if(temp<mindiv) temp=mindiv;
			image[i]/=(temp/multiplier);
		}
		smoothed=null;
	}

}
