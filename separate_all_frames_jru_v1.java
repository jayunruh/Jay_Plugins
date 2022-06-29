/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
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

public class separate_all_frames_jru_v1 implements PlugIn {
	//here we separate all of the time frames in a stack, naming them after their slice labels

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		int chans=imp.getNChannels();
		ImageStack stack=imp.getStack();
		String[] labels=new String[frames];
		for(int i=0;i<frames;i++){
			labels[i]=stack.getSliceLabel(i*chans*slices+1);
			if(labels[i]==null || labels[i].equals("")) labels[i]=""+(i+1);
		}
		Object[] pixels=jutils.stack2array(stack);
		double psize=jutils.get_psize(imp);
		imp.changes=false;
		imp.close();
		for(int i=0;i<frames;i++){
			ImageStack tempstack=new ImageStack(width,height);
			for(int j=0;j<slices*chans;j++){
				tempstack.addSlice("",pixels[i*slices*chans+j]);
			}
			ImagePlus imp2=new ImagePlus(labels[i],tempstack);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(chans,slices,1);
			jutils.set_psize(imp2,psize);
			if(chans>1) new CompositeImage(imp2,CompositeImage.COLOR).show();
			else imp2.show();
		}
	}

}
