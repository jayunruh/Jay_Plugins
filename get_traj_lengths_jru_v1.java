/*******************************************************************************
 * Copyright (c) 2019 Jay Unruh, Stowers Institute for Medical Research.
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

public class get_traj_lengths_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][] zvals=null;
		int[] npts=null;
		if(jutils.is3DPlot(iw)){
			npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
			zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];
		} else {
			npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		}
		for(int i=0;i<npts.length;i++){
			float length=getTrajLength(xvals[i],yvals[i],zvals[i],npts[i]);
			IJ.log("Trajectory"+(i+1)+" Length = "+length);
		}
	}

	public float getTrajLength(float[] xvals,float[] yvals,float[] zvals,int npts){
		float length=0.0f;
		for(int i=1;i<npts;i++){
			if(zvals!=null){
				length+=Math.sqrt((xvals[i]-xvals[i-1])*(xvals[i]-xvals[i-1])+(yvals[i]-yvals[i-1])*(yvals[i]-yvals[i-1])+(zvals[i]-zvals[i-1])*(zvals[i]-zvals[i-1]));
			} else {
				length+=Math.sqrt((xvals[i]-xvals[i-1])*(xvals[i]-xvals[i-1])+(yvals[i]-yvals[i-1])*(yvals[i]-yvals[i-1]));
			}
		}
		return length;
	}

}
