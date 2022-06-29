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

public class mask_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		GenericDialog gd=new GenericDialog("Choose Images");
		gd.addChoice("Stack",titles,titles[0]);
		gd.addChoice("Mask",titles,titles[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int index1=gd.getNextChoiceIndex();
		int index2=gd.getNextChoiceIndex();
		ImagePlus imp1=WindowManager.getImage(wList[index1]);
		ImagePlus imp2=WindowManager.getImage(wList[index2]);
		ImageStack stack=imp1.getStack();
		int width=imp1.getWidth();
		int height=imp1.getHeight();
		//int slices=stack.getSize();
		int nchans=imp1.getNChannels();
		int nslices=imp1.getNSlices();
		int nframes=imp1.getNFrames();
		ImageStack stack2=imp2.getStack();
		int slices2=stack2.getSize();
		byte[] mask=(byte[])stack2.getPixels(1);
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
					int pos=k+j*nchans+i*nslices*nchans+1;
					if(slices2==nframes*nslices*nchans) mask=(byte[])stack2.getPixels(pos); //different mask for each plane
					else if(slices2==nframes) mask=(byte[])stack2.getPixels(i+1); //same mask for all channels and slices
					else if(slices2==nslices) mask=(byte[])stack2.getPixels(j+1); //same mask for all frames and channels
					else if(slices2==nchans) mask=(byte[])stack2.getPixels(k+1); //same mask for all frames and slices
					else if(slices2==nframes*nslices) mask=(byte[])stack2.getPixels(j+i*nslices+1); //same mask for all channels
					else if(slices2==nframes*nchans) mask=(byte[])stack2.getPixels(k+i*nchans+1); //same mask for all slices
					Object inp=stack.getPixels(pos);
					if(inp instanceof float[]){
						float[] out=((float[])inp).clone();
						for(int l=0;l<width*height;l++) if(mask[l]==(byte)0) out[l]=0.0f;
						retstack.addSlice("",out);
					} else if(inp instanceof short[]) {
						short[] out=((short[])inp).clone();
						for(int l=0;l<width*height;l++) if(mask[l]==(byte)0) out[l]=(short)0;
						retstack.addSlice("",out);
					} else {
						byte[] out=((byte[])inp).clone();
						for(int l=0;l<width*height;l++) if(mask[l]==(byte)0) out[l]=(byte)0;
						retstack.addSlice("",out);
					}
				}
			}
		}
		jutils.create_hyperstack("Masked Image",retstack,imp1).show();
	}

}
