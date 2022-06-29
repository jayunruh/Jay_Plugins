/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import java.lang.reflect.*;

public class get_window_class_info_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		IJ.log("Window class: "+iw.getClass().getName());
		Object[] refmethods=jutils.getReflectionMethods(iw);
		Method[] meth=(Method[])refmethods[0];
		IJ.log("Methods:");
		for(int i=0;i<meth.length;i++){
			IJ.log(meth[i].getName());	
		}
	}

}
