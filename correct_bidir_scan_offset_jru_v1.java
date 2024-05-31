/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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

public class correct_bidir_scan_offset_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		Object[] stack2=jutils.stack2array(stack);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Shift (pixels)",2,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int shift=(int)gd.getNextNumber();
		int dtype=algutils.get_array_type(stack2[0]);
		for(int i=0;i<stack2.length;i++){
			for(int j=0;j<height;j+=2){
				Object line=algutils.get_image_row(stack2[i],width,height,j);
				for(int k=0;k<width;k++){
					algutils.setPixelVal(stack2[i],0.0f,k,j,width,height,dtype);
					int xpos=k-shift;
					if(xpos>=0 && xpos<width){
						float val=algutils.getArrVal(line,xpos,dtype);
						algutils.setPixelVal(stack2[i],val,k,j,width,height,dtype);
					}
				}
			}
		}
		imp.updateAndDraw();
	}

}
