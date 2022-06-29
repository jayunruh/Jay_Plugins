/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.io.*;
import jalgs.*;
import jguis.*;

public class test_tiff_writer_jru_v1 implements PlugIn, FrameInterface {
	int currframe;
	ImageStack stack;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		stack=imp.getStack();
		String name=imp.getTitle();
		SaveDialog sd=new SaveDialog("Save As",name,".tif");
		String dir=sd.getDirectory();
		String fname=sd.getFileName();
		if(fname==null || fname.length()==0) return;
		//(new Tiff_Writer()).saveAsTiff(dir+fname,imp);
		currframe=0;
		Tiff_Writer tw=new Tiff_Writer(imp,stack.getSize(),this);
		tw.saveAsTiffStack(dir+fname);
	}

	public Object getNextFrame(){
		currframe++;
		IJ.showProgress(currframe,stack.getSize());
		return stack.getPixels(currframe);
	}

}
