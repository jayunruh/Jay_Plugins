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
import java.awt.Frame;
import java.awt.Color;
import java.awt.AWTEvent;
import ij.text.*;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import jalgs.jseg.*;
import jalgs.*;
import jguis.*;

public class count_max_not_mask_jru_v1 implements PlugIn{
	int slices,currframe,width,height,minsize;
	float threshfraction;
	boolean com;
	String fracstat;
	findblobs fb;
	ImageStack stack;
	ImagePlus imp;

	public void run(String arg) {
		imp=WindowManager.getCurrentImage();
		width=imp.getWidth(); height=imp.getHeight();
		stack=imp.getStack();
		slices=stack.getSize();
		GenericDialog gd5=new GenericDialog("Options");
		gd5.addChoice("Threshhold Statistic",jstatistics.stats,jstatistics.stats[2]);
		gd5.showDialog(); if(gd5.wasCanceled()) return;
		fracstat=jstatistics.stats[gd5.getNextChoiceIndex()];
		threshfraction=2.0f;
		if(fracstat.equals("Max")) threshfraction=0.5f;
		fb=new findblobs(width,height,new float[]{0.0f,10000.0f,0.0f,1000.0f,10,5,10});
		float[] criteria=get_criteria(imp,fb);
		if(criteria==null) return;
		imp.setHideOverlay(true);
		fb.setcriteria(criteria);
		fb.usemaxpt=(!com);
		for(int j=0;j<slices;j++){
			float[] pixels=algutils.convert_arr_float2(stack.getPixels(j+1));
			float max=jstatistics.getstatistic(fracstat,pixels,null);
			fb.thresh=threshfraction*max;
			if(fb.thresh==0.0f) return;
			float[][] blobstats=fb.dofindblobs2(pixels,new float[width*height]);
			IJ.log(""+blobstats.length);
			float[][][] newstats=new float[2][blobstats.length][1];
			for(int i=0;i<blobstats.length;i++){
				newstats[0][i][0]=blobstats[i][0];
				newstats[1][i][0]=blobstats[i][1];
			}
			new PlotWindow4("Positions"+j,"x","y",newstats[0],newstats[1],null).draw();
			IJ.showProgress(j,slices);
		}
	}

	public float[] get_criteria(ImagePlus imp,findblobs fb){
		GenericDialog gd=new GenericDialog("Adjust Threshold");
		gd.addNumericField("Min_separation (pix)",4,0);
		gd.addNumericField("Thresh_fraction",threshfraction,5,15,null);
		gd.addNumericField("Edge_buffer (pix)",10,0);
		gd.addNumericField("Max_Blobs",1000000,0);
		gd.addCheckbox("Center_Of_Mass",true);
		//gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		//criteria are 0minarea, 1maxarea, 2searchd(nu), 3maxblobs, 4thresh, 5minsep, 6edgebuf
		float[] criteria={0.0f,10000.0f,0.0f,1000.0f,0.5f,4.0f,10.0f};
		criteria[5]=(float)gd.getNextNumber();
		threshfraction=(float)gd.getNextNumber();
		criteria[6]=(float)gd.getNextNumber();
		criteria[3]=(float)gd.getNextNumber();
		com=gd.getNextBoolean();
		return criteria;
	}

}
