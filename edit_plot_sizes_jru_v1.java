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
import jguis.*;

public class edit_plot_sizes_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		Object plot=jutils.getReflectionField(iw,"p3");
		int fontsize=(Integer)jutils.getReflectionField(plot,"fontsize");
		int shapesize=(Integer)jutils.getReflectionField(plot,"shapesize");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Font Size",fontsize,0);
		gd.addNumericField("Shape Size",shapesize,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		fontsize=(int)gd.getNextNumber();
		shapesize=(int)gd.getNextNumber();
		jutils.runReflectionMethod(plot,"setFontSize",new Object[]{fontsize});
		jutils.runReflectionMethod(plot,"setShapeSize",new Object[]{shapesize});
		jutils.runPW4VoidMethod(iw,"updatePlot");
	}

}
