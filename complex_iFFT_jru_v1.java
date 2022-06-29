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
import jalgs.jfft.*;

public class complex_iFFT_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack cstack=new ImageStack(width,height);
		po4cfft2D fft=new po4cfft2D(width,height);
		for(int j=0;j<size/2;j++){
			float[] rpixels=(float[])stack.getPixels(2*j+1);
			float[] ipixels=(float[])stack.getPixels(2*j+2);
			float[] real=new float[width*height];
			float[] im=new float[width*height];
			System.arraycopy(rpixels,0,real,0,width*height);
			System.arraycopy(ipixels,0,im,0,width*height);
			fft.docfft2D(real,im,true);
			cstack.addSlice("",(Object)real);
			cstack.addSlice("",(Object)im);
			IJ.showProgress(j,size/2);
		}
		ImagePlus imp2=new ImagePlus("iFFT",cstack);
		imp2.show();
	}

}
