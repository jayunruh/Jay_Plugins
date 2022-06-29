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
import java.awt.Frame;
import ij.plugin.*;
import jguis.*;
import ij.text.*;
import java.util.*;

public class hist_2D_2_1D_gated_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		int[] indices=(int[])jutils.runReflectionMethod(iw,"getroiindices",null);
		if(indices==null){
			IJ.log("Select Roi First");
			return;
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Y_Axis_Histogram",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean yaxis=gd.getNextBoolean();
		float[] xvals=(float[])jutils.runPW4VoidMethod(iw,"getXValues");
		float[] yvals=(float[])jutils.runPW4VoidMethod(iw,"getYValues");
		String xlab=(String)jutils.runReflectionMethod(iw,"getxLabel",null);
		String ylab=(String)jutils.runReflectionMethod(iw,"getyLabel",null);
		float[] subset=new float[indices.length];
		for(int i=0;i<indices.length;i++){
			if(yaxis) subset[i]=yvals[indices[i]];
			else subset[i]=xvals[indices[i]];
		}
		String label=xlab;  if(yaxis) label=ylab;
		new PlotWindowHist("Histogram",label,"Occurrences",subset,3).draw();
	}

}
