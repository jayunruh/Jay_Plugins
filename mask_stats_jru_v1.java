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
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Image",titles,titles[0]);
		gd2.addChoice("Mask",titles,titles[0]);
		gd2.addChoice("Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		int index1=gd2.getNextChoiceIndex();
		int index2=gd2.getNextChoiceIndex();
		int statindex=gd2.getNextChoiceIndex();
		String stat=jstatistics.stats[statindex];
		ImagePlus imp1 = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		int height = imp1.getHeight();
		int width = imp1.getWidth();

		ImageStack maskstack=imp2.getStack();
		int maskslices=maskstack.getSize();
		float[] mask=(float[])((maskstack.getProcessor(1).convertToFloat()).getPixels());
		int nis=0;
		int nnot=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]>0.0f){nis++;}
			else{nnot++;}
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
					if(mask[j]>0.0f){nis++;}
					else{nnot++;}
				}
			}
			float[] tempis=new float[nis];
			float[] tempnot=new float[nnot];
			Object temppixels=datastack.getPixels(i+1);
			int counteris=0;
			int counternot=0;
			for(int j=0;j<width*height;j++){
				if(mask[j]>0.0f){
					if(isfloat){tempis[counteris]=((float[])temppixels)[j];}
					else{tempis[counteris]=((short[])temppixels)[j]&0xffff;}
					counteris++;
				} else {
					if(isfloat){tempnot[counternot]=((float[])temppixels)[j];}
					else{tempnot[counternot]=((short[])temppixels)[j]&0xffff;}
					counternot++;
				}
			}
			statis[i]=jstatistics.getstatistic(stat,tempis,null);
			statnot[i]=jstatistics.getstatistic(stat,tempnot,null);
		}
		if(slices==1){
			IJ.log(stat+" Mask = "+statis[0]);
			IJ.log(stat+" Not Mask = "+statnot[0]);
		} else {
			PlotWindow4 pw=new PlotWindow4("Mask "+stat,"slice",stat,statis);
			pw.draw();
			pw.addPoints(statnot,true);
		}
	}

}
