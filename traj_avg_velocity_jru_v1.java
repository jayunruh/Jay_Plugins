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
import jalgs.*;

public class traj_avg_velocity_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Minimum_Length",5,0);
		gd.addNumericField("Step_Size",1,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int minpts=(int)gd.getNextNumber();
		int step=(int)gd.getNextNumber();
		if(minpts<=step) minpts=(step+1);
		int totpts=0;
		for(int i=0;i<npts.length;i++){
			if(npts[i]>=minpts) totpts+=(npts[i]-step);
		}
		float[] steps=new float[totpts];
		int counter=0;
		for(int i=0;i<npts.length;i++){
			if(npts[i]>=minpts){
				for(int j=0;j<(npts[i]-step);j++){
					steps[counter]=(float)Math.sqrt((xvals[i][j+step]-xvals[i][j])*(xvals[i][j+step]-xvals[i][j])+(yvals[i][j+step]-yvals[i][j])*(yvals[i][j+step]-yvals[i][j]));
					steps[counter]/=(float)step;
					counter++;
				}
			}
		}
		new PlotWindowHist("Step Size Histogram","Step Size","Frequency",steps,3).draw();
		float avgsize=jstatistics.getstatistic("Avg",steps,null);
		IJ.log("Avg Stepsize = "+avgsize);
	}

}
