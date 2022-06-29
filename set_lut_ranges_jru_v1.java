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

public class set_lut_ranges_jru_v1 implements PlugIn {

	public void run(String arg) {
		CompositeImage imp=(CompositeImage)WindowManager.getCurrentImage();
		int mode=imp.getMode();
		ImageStack stack=imp.getStack();
		int nch=imp.getNChannels();
		int currframe=imp.getFrame();
		int currslice=imp.getSlice();	
		int currch=imp.getChannel();
		double[][] currranges=new double[nch][2];
		for(int i=0;i<nch;i++){
			imp.setPositionWithoutUpdate(i+1,currslice,currframe);
			currranges[i][0]=imp.getProcessor().getMin();
			currranges[i][1]=imp.getProcessor().getMax();
		}
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<nch;i++){
			gd.addNumericField("Range_"+(i+1)+"_min",currranges[i][0],5,15,null);
			gd.addNumericField("Range_"+(i+1)+"_max",currranges[i][1],5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()){return;}
		for(int i=0;i<nch;i++){
			currranges[i][0]=gd.getNextNumber();
			currranges[i][1]=gd.getNextNumber();
		}
		for(int i=0;i<nch;i++){
			imp.setPositionWithoutUpdate(i+1,currslice,currframe);
			imp.setDisplayRange(currranges[i][0],currranges[i][1]);
		}
		imp.reset();
		imp.setPosition(currch,currslice,currframe);
		imp.updateChannelAndDraw();
	}

}
