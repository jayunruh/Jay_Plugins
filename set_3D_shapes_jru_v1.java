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
import jguis.*;

public class set_3D_shapes_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		int[] shapes=(int[])jutils.runPW4VoidMethod(iw,"getShapes");
		int[] colors=(int[])jutils.runPW4VoidMethod(iw,"getColors");
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Shape",Plotter.shape_names,Plotter.shape_names[1]);
		gd.addCheckbox("Change Color",false);
		gd.addChoice("Color",Plot4.color_names,Plot4.color_names[0]);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int newshape=gd.getNextChoiceIndex();
		boolean changecolors=gd.getNextBoolean();
		int newcolor=gd.getNextChoiceIndex();
		for(int i=0;i<shapes.length;i++) shapes[i]=newshape;
		if(changecolors){
			for(int i=0;i<colors.length;i++) colors[i]=newcolor;
		}
		jutils.runPW4VoidMethod(iw,"updatePlot");
	}

}
