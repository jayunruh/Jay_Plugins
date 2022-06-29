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
import jalgs.*;
import jalgs.jseg.*;

public class track_cells_phase_jru_v1 implements PlugIn {
	//this plugin attempts to track cells in phase contrast using a combination of canny hysteresis edge detection and dilations and erosions
	//lower resolution phase tends to close gaps in membranes

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean all=false;
		gd.addCheckbox("Analyze All Slices",all);
		float thresh1=100.0f;
		gd.addNumericField("Threshhold 1",thresh1,5,10,null);
		float thresh2=50.0f;
		gd.addNumericField("Hysteresis Threshhold",thresh2,5,10,null);
		int minsize=1000;
		gd.addNumericField("Min Cell Size",minsize,0);
		int maxsize=100000;
		gd.addNumericField("Max Cell Size",maxsize,0);
		int ndilations=3;
		gd.addNumericField("Number of Dilations",ndilations,0);
		int nerosions=2;
		gd.addNumericField("Number of Erosions",nerosions,0);
		gd.addCheckbox("Output_Canny",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		all=gd.getNextBoolean();
		thresh1=(float)gd.getNextNumber();
		thresh2=(float)gd.getNextNumber();
		minsize=(int)gd.getNextNumber();
		maxsize=(int)gd.getNextNumber();
		ndilations=(int)gd.getNextNumber();
		nerosions=(int)gd.getNextNumber();
		boolean outcanny1=gd.getNextBoolean();
		
		//start by getting a gradient image in x and y with a sobel edge filter
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		canny cannyclass=new canny(thresh1,thresh2,width,height);
		ImageStack outstack=new ImageStack(width,height);
		ImageStack outcanny=new ImageStack(width,height);
		findblobs3 fb=new findblobs3(width,height);
		binary_processing bp=new binary_processing(width,height);
		if(!all){slices=1;}
		for(int i=0;i<slices;i++){
			float[] temp=(float[])(stack.getProcessor(i+1).convertToFloat()).getPixels();
			byte[] hisedgeimage=cannyclass.find_edges(temp);
			//outcanny.addSlice("",hisedgeimage.clone());
			for(int j=0;j<ndilations;j++){
				bp.dilate(hisedgeimage);
			}
			for(int j=0;j<nerosions;j++){
				bp.erode(hisedgeimage);
			}
			if(outcanny1) outcanny.addSlice("",hisedgeimage.clone());
			bp.fillholes(hisedgeimage);
			float[] blobs=fb.dofindblobs(hisedgeimage);
			int[] limits={minsize,maxsize};
			fb.filter_area(blobs,limits);
			outstack.addSlice("",fb.tobinary(blobs,false));
			//outstack.addSlice("",hisedgeimage);
			IJ.showProgress(i,slices);
		}
		new ImagePlus("Cell Track Image",outstack).show();
		if(outcanny1) new ImagePlus("Canny Image",outcanny).show();
	}

}
