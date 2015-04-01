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
import jguis.*;

public class mask_stats_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length+1];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		titles[wList.length]="null";
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Image",titles,titles[0]);
		gd2.addChoice("Mask",titles,titles[0]);
		gd2.addChoice("Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd2.addCheckbox("Report_Outside_Stats?",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		int index1=gd2.getNextChoiceIndex();
		int index2=gd2.getNextChoiceIndex();
		int statindex=gd2.getNextChoiceIndex();
		boolean reportoutside=gd2.getNextBoolean();
		String stat=jstatistics.stats[statindex];
		ImagePlus imp1 = WindowManager.getImage(wList[index1]);
		ImagePlus imp2=null;
		if(index2<wList.length) imp2 = WindowManager.getImage(wList[index2]);
		int height = imp1.getHeight();
		int width = imp1.getWidth();
		float[] extras=jutils.getStatsOptions(stat);
		
		Roi roi=imp1.getRoi();
		boolean[] roimask=new boolean[width*height];
		for(int i=0;i<roimask.length;i++) roimask[i]=true;
		if(roi!=null) roimask=jstatistics.poly2mask(roi.getPolygon(),width,height);

		ImageStack maskstack=null;
		int maskslices=1;
		float[] mask=new float[width*height];
		for(int i=0;i<mask.length;i++) mask[i]=1.0f;
		int nis=0; int nnot=0;
		if(imp2!=null){
			//if we have a mask, the roi limits the total analysis size
			maskstack=imp2.getStack();
			maskslices=maskstack.getSize();
			float[] temp=(float[])((maskstack.getProcessor(1).convertToFloat()).getPixels());
			mask=temp.clone();
			for(int i=0;i<width*height;i++){
				if(roimask[i]){
					if(mask[i]>0.0f) nis++;
					else nnot++;
				} else {
					mask[i]=0.0f;
				}
			}
		} else {
			//if we don't have a mask, the roi is the mask
			for(int i=0;i<width*height;i++){
				if(roimask[i]){
					mask[i]=1.0f;
					nis++;
				} else {
					mask[i]=0.0f;
					nnot++;
				}
				roimask[i]=true;
			}
		}

		ImageStack datastack=imp1.getStack();
		int slices=datastack.getSize();
		boolean isfloat=(datastack.getPixels(1) instanceof float[]);
		float[] statis=new float[slices];
		float[] statnot=new float[slices];
		for(int i=0;i<slices;i++){
			if(maskslices>1){
				mask=(float[])((maskstack.getProcessor(i+1).convertToFloat()).getPixels());
				nis=0;
				nnot=0;
				for(int j=0;j<width*height;j++){
					if(roimask[j]){
						if(mask[j]>0.0f) nis++;
						else nnot++;
					}
				}
			}
			float[] tempis=new float[nis];
			float[] tempnot=null;
			if(reportoutside) tempnot=new float[nnot];
			Object temppixels=datastack.getPixels(i+1);
			int counteris=0;
			int counternot=0;
			for(int j=0;j<width*height;j++){
				if(roimask[j]){
					if(mask[j]>0.0f){
						if(isfloat){tempis[counteris]=((float[])temppixels)[j];}
						else{tempis[counteris]=((short[])temppixels)[j]&0xffff;}
						counteris++;
					} else {
						if(reportoutside){
							if(isfloat){tempnot[counternot]=((float[])temppixels)[j];}
							else{tempnot[counternot]=((short[])temppixels)[j]&0xffff;}
							counternot++;
						}
					}
				}
			}
			statis[i]=jstatistics.getstatistic(stat,tempis,extras);
			if(reportoutside) statnot[i]=jstatistics.getstatistic(stat,tempnot,extras);
		}
		if(slices==1){
			IJ.log(stat+" Mask = "+statis[0]);
			if(reportoutside) IJ.log(stat+" Not Mask = "+statnot[0]);
		} else {
			PlotWindow4 pw=new PlotWindow4("Mask "+stat,"slice",stat,statis);
			pw.draw();
			if(reportoutside) pw.addPoints(statnot,true);
		}
	}

}
