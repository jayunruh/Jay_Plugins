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
import jguis.*;
import jalgs.*;

public class roi_outside_fill_avg_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi roi=imp.getRoi();
		int width=imp.getWidth();
		int height=imp.getHeight();
		boolean[] mask=jutils.roi2mask(roi,width,height);
		ImageStack stack=imp.getStack();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Fill_Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd2.addCheckbox("Fill_Inside",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		int statindex=gd2.getNextChoiceIndex();
		boolean fillinside=gd2.getNextBoolean();
		if(fillinside) for(int i=0;i<width*height;i++){mask[i]=(!mask[i]);}
		String stat=jstatistics.stats[statindex];
		if(frames==1){
			frames=slices;
			slices=1;
		}
		boolean analyzeall=true;
		if(frames>1 && (slices>1 || channels>1)){
			GenericDialog gd=new GenericDialog("Options");
			gd.addCheckbox("analyze_all_slices_&_channels?",true);
			gd.showDialog();
			if(gd.wasCanceled()){return;}
			analyzeall=gd.getNextBoolean();
		}
		if(analyzeall){
			for(int i=0;i<stack.getSize();i++){
				Object pixels=stack.getPixels(i+1);
				float roiavg=jstatistics.getstatistic(stat,pixels,width,height,mask,null);
				if(pixels instanceof float[]){
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((float[])pixels)[j]=roiavg;
						}
					}
				} else if(pixels instanceof short[]) {
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((short[])pixels)[j]=(short)roiavg;
						}
					}
				} else {
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((byte[])pixels)[j]=(byte)roiavg;
						}
					}
				}
			}
		} else {
			int currchan=imp.getChannel()-1;
			int currslice=imp.getSlice()-1;
			if(imp.getNFrames()==1) currslice=0;
			for(int i=0;i<frames;i++){
				Object pixels=stack.getPixels(i*slices*channels+currslice*channels+currchan+1);
				float roiavg=jstatistics.getstatistic(stat,pixels,width,height,mask,null);
				if(pixels instanceof float[]){
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((float[])pixels)[j]=roiavg;
						}
					}
				} else if(pixels instanceof short[]) {
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((short[])pixels)[j]=(short)roiavg;
						}
					}
				} else {
					for(int j=0;j<width*height;j++){
						if(!mask[j]){
							((byte[])pixels)[j]=(byte)roiavg;
						}
					}
				}
			}
		}
		imp.updateAndDraw();		
	}

}
