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
		int slices=stack.getSize();
		ImageStack stack2=imp2.getStack();
		int slices2=stack2.getSize();
		byte[] mask=(byte[])stack2.getPixels(1);
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<slices;i++){
			if(slices2>1){
				mask=(byte[])stack2.getPixels(i+1);
			}
			Object inp=stack.getPixels(i+1);
			if(inp instanceof float[]){
				float[] out=new float[width*height];
				for(int j=0;j<width*height;j++){
					int temp=mask[j]&0xff;
					if(temp!=0){out[j]=((float[])inp)[j];}
				}
				retstack.addSlice("",out);
			} else {
				short[] out=new short[width*height];
				for(int j=0;j<width*height;j++){
					int temp=mask[j]&0xff;
					if(temp!=0){out[j]=((short[])inp)[j];}
				}
				retstack.addSlice("",out);
			}
		}
		(new ImagePlus("Masked Stack",retstack)).show();
	}

}
