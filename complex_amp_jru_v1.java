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

public class complex_amp_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack is=imp.getStack();
		int size=is.getSize();
		ImageStack ampstack=new ImageStack(width,height);
		for(int j=0;j<size/2;j++){
			float[] real=(float[])is.getPixels(2*j+1);
			float[] im=(float[])is.getPixels(2*j+2);
			float[] amp=new float[width*height];
			for(int i=0;i<width*height;i++){
				amp[i]=(float)Math.sqrt(real[i]*real[i]+im[i]*im[i]);
			}
			ampstack.addSlice("",amp);
			IJ.showProgress(j,size/2);
		}
		ImagePlus imp2=new ImagePlus("Amp",ampstack);
		imp2.show();
	}

}
