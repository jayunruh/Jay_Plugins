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
import java.io.*;
import jalgs.*;
import jguis.*;

public class stack_statistics_jru_v2 implements PlugIn {
		//this plugin calculates statistics for each pixel in a stack
		int avgframes,index;
		boolean avgall,batchmode;

	public void run(String arg) {
		GenericDialog gd = new GenericDialog("Options");
		String[] stats=jstatistics.stats;
		gd.addChoice("Statistic",stats,stats[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String stat=stats[gd.getNextChoiceIndex()];
		float[] extras=jutils.getStatsOptions(stat);
		//get the current image and its info
		ImagePlus imp=WindowManager.getCurrentImage();
		String otitle=imp.getTitle();
		int namelen=otitle.lastIndexOf(".");
		String trunctitle;
		if(namelen>0){
			trunctitle=otitle.substring(0,namelen);
		} else{
			trunctitle=otitle;
		}
		ImageStack stack = imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		float[][] result=new float[channels][];
		for(int i=0;i<channels;i++){
			result[i]=jutils.get3DProjZStat(stack,0,i,1,slices*frames,channels,stat);
		}
		ImageStack retstack=new ImageStack(stack.getWidth(),stack.getHeight());
		for(int i=0;i<result.length;i++) retstack.addSlice("",result[i]);
		ImagePlus imp2=new ImagePlus(trunctitle+"-"+stat,retstack);
		imp2.copyScale(imp);
		if(channels>1 && imp.isComposite()){
			CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
			LUT[] lut=((CompositeImage)imp).getLuts();
			ci.setLuts(lut);
			ci.resetDisplayRanges();
			ci.show();
		} else {
			imp2.show();
		}
	}

}
