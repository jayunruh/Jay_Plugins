/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.measure.*;

public class set_min_max_all_chan_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int currframe=imp.getFrame();
		int currslice=imp.getSlice();	
		int currch=imp.getChannel();
		int nch=imp.getNChannels();
		double min=imp.getProcessor().getMin();
		double max=imp.getProcessor().getMax();
		for(int i=0;i<nch;i++){
			imp.setPositionWithoutUpdate(i+1,currslice,currframe);
			imp.setDisplayRange(min,max);
		}
		((CompositeImage)imp).reset();
		imp.setPosition(currch,currslice,currframe);
		imp.updateChannelAndDraw();
	}

}
