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
import jguis.*;

public class separate_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		String[] labels=(String[])jutils.runPW4VoidMethod(iw,"getAllLabels");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		for(int i=0;i<yvals.length;i++){
			float[] tempx=new float[npts[i]];
			float[] tempy=new float[npts[i]];
			System.arraycopy(xvals[i],0,tempx,0,npts[i]);
			System.arraycopy(yvals[i],0,tempy,0,npts[i]);
			new PlotWindow4(labels[0]+(i+1),labels[1],labels[2],tempx,tempy).draw();
		}
	}
}
