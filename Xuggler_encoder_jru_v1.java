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
import ij.io.*;
import java.io.*;
import jguis.*;

public class Xuggler_encoder_jru_v1 implements PlugIn, FrameInterface {
	int frameIterator;
	ImageStack stack;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		stack=imp.getStack();
		frameIterator=0;
		Xuggler_encoder e=new Xuggler_encoder(imp,stack.getSize(),this,5.0);
		SaveDialog sd=new SaveDialog("Save Movie File",arg,".mp4");
		String dir=sd.getDirectory();
		String name=sd.getFileName();
		if(name==null || name.length()==0){
			return;
		}
		e.saveAsMovie(dir+File.separator+name);
	}

	public Object getNextFrame(){
		frameIterator++;
		return stack.getProcessor(frameIterator).convertToRGB().getPixels();
	}

}
