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
import jalgs.jseg.*;

public class make_edge_connected_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		binary_processing bp=new binary_processing(width,height);
		for(int i=0;i<stack.getSize();i++){
			byte[] pixels=(byte[])stack.getProcessor(i+1).getPixels();
			bp.make_edge_connected(pixels);
			IJ.showProgress(i,stack.getSize());
		}
		imp.updateAndDraw();
	}

}
