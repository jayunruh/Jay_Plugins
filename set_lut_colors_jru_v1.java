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
import java.awt.image.*;

public class set_lut_colors_jru_v1 implements PlugIn {

	public void run(String arg) {
		CompositeImage imp=(CompositeImage)WindowManager.getCurrentImage();
		int mode=imp.getMode();
		ImageStack stack=imp.getStack();
		LUT[] luts=imp.getLuts();
		int[] colors=new int[luts.length];
		for(int i=0;i<luts.length;i++){
			colors[i]=jutils.get_LUT_color_index(luts[i]);
		}
		GenericDialog gd=new GenericDialog("Options");
		String[] colornames={"white","blue","green","red","magenta","cyan","yellow","orange"};
		for(int i=0;i<luts.length;i++) gd.addChoice("Color"+(i+1),colornames,colornames[colors[i]]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		for(int i=0;i<luts.length;i++){
			int index=gd.getNextChoiceIndex();
			if(index!=colors[i]){
				luts[i]=jutils.get_lut_for_color(jutils.colors[index]);
			}
		}
		//here we set a channel lookup table explicitly
		imp.setStack(stack);
		imp.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		((CompositeImage)imp).setLuts(luts);
		imp.setSlice(1);
	}

}
