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
import jguis.*;
import jalgs.*;

public class subtract_line_from_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){
				titles[i]=imp.getTitle();
			} else {
				titles[i]="NA";
			}
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Plot Window",titles,titles[0]);
		gd.addChoice("Stack",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);

		ImageWindow iw=imp.getWindow();
		float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
		int length=yvals.length;
		ImageStack stack=imp2.getStack();
		int width=imp2.getWidth(); int height=imp2.getHeight();
		ImageStack retstack=new ImageStack(width,height);
		String name="Subtracted";
		if(stack.getSize()==length){
			for(int i=0;i<length;i++){
				float[] pixels=algutils.convert_arr_float(stack.getPixels(i+1));
				for(int j=0;j<width*height;j++){
					pixels[j]-=yvals[i];
				}
				retstack.addSlice("",pixels);
			}
			(new ImagePlus(name,retstack)).show();
		} else {
			int nframes=imp2.getNFrames();
			int nslices=imp2.getNSlices();
			int nchannels=imp2.getNChannels();
			if(length==nslices){
				//assume dividing along the z axis
				int counter=0;
				for(int i=0;i<nframes;i++){
					for(int j=0;j<nslices;j++){
						for(int k=0;k<nchannels;k++){
							float[] pixels=algutils.convert_arr_float(stack.getPixels(counter+1)); counter++;
							float sub=yvals[j];
							for(int m=0;m<width*height;m++) pixels[m]-=sub;
							retstack.addSlice("",pixels);
						}
					}
				}				
			} else if(length==nframes){
				//assume dividing along the t axis
				int counter=0;
				for(int i=0;i<nframes;i++){
					for(int j=0;j<nslices;j++){
						for(int k=0;k<nchannels;k++){
							float[] pixels=algutils.convert_arr_float(stack.getPixels(counter+1)); counter++;
							float sub=yvals[i];
							for(int m=0;m<width*height;m++) pixels[m]-=sub;
							retstack.addSlice("",pixels);
						}
					}
				}	
			} else {
				//assume dividing along the c axis
				int counter=0;
				for(int i=0;i<nframes;i++){
					for(int j=0;j<nslices;j++){
						for(int k=0;k<nchannels;k++){
							float[] pixels=algutils.convert_arr_float(stack.getPixels(counter+1)); counter++;
							float sub=yvals[k];
							for(int m=0;m<width*height;m++) pixels[m]-=sub;
							retstack.addSlice("",pixels);
						}
					}
				}	
			}
			jutils.create_hyperstack(name,retstack,imp2).show();
		}

	}

}
