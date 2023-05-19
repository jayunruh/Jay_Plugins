/*******************************************************************************
 * Copyright (c) 2022 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;

public class avg_mirror_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin averages each series with it's mirror image
		//the mirror image is simple a reversed version of the data set
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][] outy=new float[yvals.length][yvals[0].length];
		for(int i=0;i<npts.length;i++){
			float[] yvals2=(float[])algutils.get_subarray(yvals[i],0,npts[i]);
			for(int j=0;j<yvals2.length;j++){
				outy[i][j]=0.5f*(yvals2[yvals2.length-j-1]+yvals2[j]);
			}
		}
		String xlab=(String)jutils.runPW4VoidMethod(iw,"getxLabel");
		String ylab=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		new PlotWindow4("Mirror_Avg",xlab,ylab,xvals,outy,npts).draw();
	}

}
