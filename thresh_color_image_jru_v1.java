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
import jguis.*;
import jalgs.*;

public class thresh_color_image_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		Roi roi=imp.getRoi();
		if(roi==null){
			//use the max sum pixel for the first slice and the pixels surrounding it as the background roi
			byte[][] rgb=jutils.intval2rgb((int[])stack.getPixels(1));
			float[] sumproj=algutils.get_stack_proj_stat("Avg",rgb,width,height,3,null);
			float maxpos=jstatistics.getstatistic("maxpos",sumproj,null);
			int maxy=(int)(maxpos/(float)width);
			int maxx=(int)maxpos-maxy*width;
			if(maxx<1) maxx=1;  if(maxx>=(width-1)) maxx=width-2;
			if(maxy<1) maxy=1;  if(maxy>=(height-1)) maxy=height-2;
			roi=new Roi(maxx-1,maxy-1,3,3);
		}
		boolean[] mask=jutils.roi2mask(roi,width,height);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Red_Upper_Limit",1.3,5,15,null);
		gd.addNumericField("Red_Lower_Limit",0,5,15,null);
		gd.addNumericField("Green_Upper_Limit",0.9,5,15,null);
		gd.addNumericField("Green_Lower_Limit",0,5,15,null);
		gd.addNumericField("Blue_Upper_Limit",0.9,5,15,null);
		gd.addNumericField("Blue_Lower_Limit",0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float redup=(float)gd.getNextNumber();
		float reddown=(float)gd.getNextNumber();
		float greenup=(float)gd.getNextNumber();
		float greendown=(float)gd.getNextNumber();
		float blueup=(float)gd.getNextNumber();
		float bluedown=(float)gd.getNextNumber();

		ImageStack threshed=new ImageStack(width,height);
		for(int i=0;i<stack.getSize();i++){
			int[] pix=(int[])stack.getPixels(i+1);
			byte[][] rgb=jutils.intval2rgb(pix);
			float[] backs=jstatistics.getspectrum("Avg",rgb,width,height,mask,null);
			float[] mins=jstatistics.getspectrum("Min",rgb,null);
			float redup2=mins[0]+redup*(backs[0]-mins[0]);
			float reddown2=mins[0]+reddown*(backs[0]-mins[0]);
			float greenup2=mins[1]+greenup*(backs[1]-mins[1]);
			float greendown2=mins[1]+greendown*(backs[1]-mins[1]);
			float blueup2=mins[2]+blueup*(backs[2]-mins[2]);
			float bluedown2=mins[2]+bluedown*(backs[2]-mins[2]);
			IJ.log("slice "+(i+1)+" red thresh = \t"+reddown2+" , "+redup2);
			IJ.log("slice "+(i+1)+" green thresh = \t"+greendown2+" , "+greenup2);
			IJ.log("slice "+(i+1)+" blue thresh = \t"+bluedown2+" , "+blueup2);
			byte[] thresh=new byte[width*height];
			for(int j=0;j<pix.length;j++){
				int[] rgb2={rgb[0][j]&0xff,rgb[1][j]&0xff,rgb[2][j]&0xff};
				if(rgb2[0]>=reddown2 && rgb2[0]<redup2 && rgb2[1]>=greendown2 && rgb2[1]<greenup2 && rgb2[2]>=bluedown2 && rgb2[2]<blueup2){
					thresh[j]=(byte)255;
				}
			}
			threshed.addSlice("",thresh);
		}
		new ImagePlus("Thresholded",threshed).show();
	}

}
