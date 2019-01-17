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

public class traj_rescale_x_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[] limits=(float[])jutils.runPW4VoidMethod(iw,"getLimits");
		GenericDialog gd=new GenericDialog("Options");
		float xmin=limits[0];
		gd.addNumericField("X_Min",xmin,5,10,null);
		float xmax=limits[1];
		gd.addNumericField("X_Max",xmax,5,10,null);
		gd.addCheckbox("Log_X",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		xmin=(float)gd.getNextNumber();
		xmax=(float)gd.getNextNumber();
		boolean logx=gd.getNextBoolean();
		if(!logx){
			for(int i=0;i<xvals.length;i++){
				float xinc=(xmax-xmin)/(float)(npts[i]-1);
				for(int j=0;j<npts[i];j++){
					xvals[i][j]=xmin+xinc*(float)j;
				}
			}
		} else {
			float logxmin=(float)Math.log(xmin);
			float logxmax=(float)Math.log(xmax);
			for(int i=0;i<xvals.length;i++){
				float logxinc=(logxmax-logxmin)/(float)(npts[i]-1);
				for(int j=0;j<npts[i];j++){
					float logval=logxmin+logxinc*(float)j;
					xvals[i][j]=(float)Math.exp(logval);
				}
			}
		}
		jutils.runPW4VoidMethod(iw,"xautoscale");
	}

}
