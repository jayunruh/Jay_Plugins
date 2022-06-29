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

public class overlay_image_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		GenericDialog gd=new GenericDialog("Choose Images");
		gd.addChoice("Image",titles,titles[0]);
		gd.addChoice("Overlay",titles,titles[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int index1=gd.getNextChoiceIndex();
		int index2=gd.getNextChoiceIndex();
		ImagePlus imp1=WindowManager.getImage(wList[index1]);
		ImagePlus imp2=WindowManager.getImage(wList[index2]);
		ImageStack stack=imp1.getStack();
		if(stack.getSize()>1){
			int size=stack.getSize();
			ImageStack stack2=imp2.getStack();
			ImageStack resstack=new ImageStack(imp1.getWidth(),imp1.getHeight());
			for(int k=0;k<size;k++){
				int[] pixels=(int[])stack.getProcessor(k+1).convertToRGB().getPixels();
				int pos=k;
				if(pos>=stack2.getSize()) pos=stack2.getSize()-1;
				int[] overlay=(int[])stack2.getProcessor(pos+1).convertToRGB().getPixels();
				int width1=imp1.getWidth();
				int height1=imp1.getHeight();
				int width2=imp2.getWidth();
				int height2=imp2.getHeight();
				int[] result=pixels.clone();
				for(int i=0;i<height2;i++){
					for(int j=0;j<width2;j++){
						int temp=overlay[j+i*width2]&0x00ffffff;
						if(temp!=0x00000000) result[j+i*width1]=overlay[j+i*width2];
					}
				}
				resstack.addSlice("",result);
			}
			new ImagePlus("Overlay",resstack).show();
		} else {
			int[] pixels=(int[])imp1.getProcessor().convertToRGB().getPixels();
			int[] overlay=(int[])imp2.getProcessor().convertToRGB().getPixels();
			int width1=imp1.getWidth();
			int height1=imp1.getHeight();
			int width2=imp2.getWidth();
			int height2=imp2.getHeight();
			int[] result=pixels.clone();
			for(int i=0;i<height2;i++){
				for(int j=0;j<width2;j++){
					int temp=overlay[j+i*width2]&0x00ffffff;
					if(temp!=0x00000000) result[j+i*width1]=overlay[j+i*width2];
				}
			}
			new ImagePlus("Overlay",new ColorProcessor(width1,height1,result)).show();
		}
	}

}
