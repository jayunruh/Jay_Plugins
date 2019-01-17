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
import ij.text.*;

public class show_values_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageProcessor ip=imp.getProcessor();
		float[] pixels=(float[])ip.getPixels();
		int width=imp.getWidth();
		int height=imp.getHeight();
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<height;i++){
			sb.append(""+pixels[i*width]);
			for(int j=1;j<width;j++){
				sb.append("\t "+pixels[j+i*width]);
			}
			sb.append("\n");
		}
		TextWindow tw=new TextWindow("Image Values",sb.toString(),200,400);
	}

}
