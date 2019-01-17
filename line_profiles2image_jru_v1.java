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
import ij.util.*;
import jguis.*;

public class line_profiles2image_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int width=yvals[0].length;
		int height=yvals.length;
		float[] carpet=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				carpet[j+i*width]=yvals[i][j];
			}
		}
		new ImagePlus("Profile Carpet",new FloatProcessor(width,height,carpet,null)).show();
	}

}
