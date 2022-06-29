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

public class autoscale_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		float upper=95.0f;
		gd.addNumericField("Upper Threshhold Percentage",upper,5,15,null);
		float lower=5.0f;
		gd.addNumericField("Lower Threshhold Percentage",lower,5,15,null);
		int maxval=235;
		gd.addNumericField("Maximum Value (<=255)",maxval,0);
		int minval=0;
		gd.addNumericField("Minimum Value (>=0)",minval,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		upper=(float)gd.getNextNumber();
		lower=(float)gd.getNextNumber();
		maxval=(int)gd.getNextNumber();
		minval=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		ImageStack stack2=new ImageStack(width,height);
		String[] labels=stack.getSliceLabels();
		for(int i=0;i<slices;i++){
			float[] pixels=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
			byte[] newpixels=autoscale_image(pixels,upper,lower,maxval,minval);
			stack2.addSlice(labels[i],newpixels);
			IJ.showProgress(i,slices);
		}
		ImagePlus imp2=new ImagePlus("Autoscaled Stack",stack2);
		imp2.copyScale(imp);
		imp2.show();
	}

	public byte[] autoscale_image(float[] image,float upper,float lower,int maxval,int minval){
		float[] temp={lower};
		float min=jstatistics.getstatistic("Percentile",image,temp);
		temp[0]=upper;
		float max=jstatistics.getstatistic("Percentile",image,temp);
		byte[] newimage=new byte[image.length];
		float range=(float)(maxval-minval);
		for(int i=0;i<image.length;i++){
			int scaled=minval+(int)(range*(image[i]-min)/(max-min));
			if(scaled<0){scaled=0;}
			if(scaled>255){scaled=255;}
			newimage[i]=(byte)scaled;
		}
		return newimage;
	}

}
