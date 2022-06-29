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
import jalgs.jseg.*;

public class convert_to_mask_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int currslice=imp.getCurrentSlice();
		ImageStack stack2=imp.createEmptyStack();
		for(int i=0;i<stack.getSize();i++){
			String label=stack.getSliceLabel(i+1);
			Object pixels=stack.getPixels(i+1);
			byte[] mask=findblobs3.threshimage(pixels,0.5f);
			stack2.addSlice(label,mask);
		}
		stack2.setColorModel(LookUpTable.createGrayscaleColorModel(false));
		imp.setStack(stack2);
		imp.setSlice(currslice);
		imp.setCalibration(imp.getCalibration());
	}

}
