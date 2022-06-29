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

public class channel_alignment_z_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		if(slices==1){slices=frames; frames=1;}
		ImageStack stack=imp.getStack();
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<channels;i++){
			gd.addNumericField("Ch"+(i+1)+"_shift (slices)",0.0,5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()) return;
		float[] xvals=new float[channels];
		for(int i=0;i<channels;i++) xvals[i]=(float)gd.getNextNumber();
		float xmin=xvals[0]; float xmax=xvals[0];
		for(int i=1;i<xvals.length;i++){
			if(xvals[i]<xmin){xmin=xvals[i];}
			if(xvals[i]>xmax){xmax=xvals[i];}
		}
		int intxmax=(int)xmax;
		int intxmin=(int)xmin;
		int newslices=slices+intxmax-intxmin+1;
		int typeindex;
		Object pixels=stack.getPixels(1);
		if(pixels instanceof float[]){
			typeindex=0;
		} else {
			if(pixels instanceof short[]){
				typeindex=1;
			} else {
				typeindex=2;
			}
		}
		ImageStack retstack=new ImageStack(width,height);
		if(typeindex==0){for(int i=0;i<frames*newslices*channels;i++) retstack.addSlice("",new float[width*height]);}
		if(typeindex==1){for(int i=0;i<frames*newslices*channels;i++) retstack.addSlice("",new short[width*height]);}
		if(typeindex==2){for(int i=0;i<frames*newslices*channels;i++) retstack.addSlice("",new byte[width*height]);}
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					int index5=j+(int)xvals[k]-intxmin+1;
					Object oldimage=jutils.get3DSlice(stack,i,j,k,frames,slices,channels);
					jutils.set3DSlice(retstack,algutils.clone_obj_array(oldimage),i,index5,k,frames,newslices,channels);
				}
				IJ.showProgress(i*slices+j,frames*slices);
			}
			IJ.showStatus("frame"+i);
		}
		jutils.create_hyperstack("Realigned",retstack,imp,frames,newslices,channels).show();
	}

}
