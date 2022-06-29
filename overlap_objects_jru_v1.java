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

public class overlap_objects_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int width=imp.getWidth(); int height=imp.getHeight(); int slices=imp.getStackSize();
		findblobs3 fb=new findblobs3(width,height);
		for(int i=0;i<slices/2;i++){
			byte[] pix1=(byte[])stack.getPixels(2*i+1);
			byte[] pix2=(byte[])stack.getPixels(2*i+2);
			float[] obj1=fb.dofindblobs(pix1);
			float[] obj2=fb.dofindblobs(pix2);
			int[] stats=fb.overlap_objects(obj1,obj2);
			IJ.log("Slice "+(i+1));
			IJ.log("# Overlaping = "+stats[0]);
			IJ.log("# Nonoverlaping ch1 = "+stats[1]);
			IJ.log("# Nonoverlaping ch2 = "+stats[2]);
		}
	}

}
