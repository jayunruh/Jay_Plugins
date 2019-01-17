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

public class roi_average_divide_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Subtract_1",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean sub=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageProcessor ip = imp.getProcessor();
		Rectangle r = ip.getRoi();
		ImageStack stack = imp.getStack();
		int slices=stack.getSize();
		ImageStack result_stack=new ImageStack(width,height);
		for(int i=0;i<slices;i++){
			if(ip instanceof FloatProcessor){
				float[] pixels = (float[])stack.getPixels(i+1);
				float[] temp=new float[width*height];
				float avg=0.0f;
				for(int j=0;j<r.height;j++){
					for(int k=0;k<r.width;k++){
						avg+=pixels[r.y*width+r.x+j*width+k]/(float)(r.width*r.height);
					}
				}
				for(int j=0;j<width*height;j++){
					temp[j]=pixels[j]/avg;
					if(sub) temp[j]-=1.0f;
				}
				result_stack.addSlice("",(Object)temp);
			}
			if(ip instanceof ShortProcessor){
				short[] pixels = (short[])stack.getPixels(i+1);
				float[] temp=new float[width*height];
				float avg=0.0f;
				for(int j=0;j<r.height;j++){
					for(int k=0;k<r.width;k++){
						float temp2=pixels[r.y*width+r.x+j*width+k]&0xffff;
						avg+=temp2/((float)(r.width*r.height));
					}
				}
				for(int j=0;j<width*height;j++){
					float temp2=pixels[j]&0xffff;
					temp[j]=temp2/avg;
					if(sub) temp[j]-=1.0f;
				}
				result_stack.addSlice("",(Object)temp);
			}
			IJ.showProgress(i,slices);
		}
		ImagePlus imp2  = new ImagePlus("Detrended",result_stack);
		imp2.show();
	}

}
