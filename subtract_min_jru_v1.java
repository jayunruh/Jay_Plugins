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

public class subtract_min_jru_v1 implements PlugIn {
	//this plugin subtracts the desired minimum value from each pixel, setting negative values to zero

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageProcessor ip=imp.getProcessor();
		GenericDialog gd = new GenericDialog("Set min");
		float min=(float)ip.getMin();
		gd.addNumericField("Min",min,5,20,null);
		gd.addCheckbox("Subtract_All_Slices",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		min = (float)gd.getNextNumber();
		boolean allslices=gd.getNextBoolean();

		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		if(allslices){
			if(ip instanceof FloatProcessor){
				for(int j=0;j<slices;j++){
					float[] pixels=(float[])stack.getPixels(j+1);
					for(int i=0;i<pixels.length;i++){
						if(Float.isNaN(pixels[i])){pixels[i]=0.0f;}
						if(pixels[i]<min){pixels[i]=0.0f;}
						else{pixels[i]-=min;}
					}
				}
			} else {
				for(int j=0;j<slices;j++){
					short[] pixels=(short[])stack.getPixels(j+1);
					for(int i=0;i<pixels.length;i++){
						float temp=(float)(pixels[i]&0xffff);
						if(temp<min){pixels[i]=(short)0;}
						else{pixels[i]=(short)(temp-min);}
					}
				}
			}
		} else {
			if(ip instanceof FloatProcessor){
				float[] pixels=(float[])ip.getPixels();
				for(int i=0;i<pixels.length;i++){
					if(Float.isNaN(pixels[i])){pixels[i]=0.0f;}
					if(pixels[i]<min){pixels[i]=0.0f;}
					else{pixels[i]-=min;}
				}
			} else {
				short[] pixels=(short[])ip.getPixels();
				for(int i=0;i<pixels.length;i++){
					float temp=(float)(pixels[i]&0xffff);
					if(temp<min){pixels[i]=(short)0;}
					else{pixels[i]=(short)(temp-min);}
				}
			}
		}
		imp.updateAndDraw();
	}

	float getmin(float[] farr){
		float min;
		min=farr[0];
		for(int i=1;i<farr.length;i++){
			if(farr[i]<min){min=farr[i];}
		}
		return min;
	}
	
	float getmax(float[] farr){
		float max;
		max=farr[0];
		for(int i=1;i<farr.length;i++){
			if(farr[i]>max){max=farr[i];}
		}
		return max;
	}

}
