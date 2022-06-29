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

public class resample_plot_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float xinc=xvals[0][1]-xvals[0][0];
		float[] limits=(float[])jutils.runPW4VoidMethod(iw,"getLimits");
		String xlabel=(String)jutils.runPW4VoidMethod(iw,"getxLabel");
		String ylabel=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		float xstart=limits[0];
		float xend=limits[1];
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Increment",xinc,5,15,null);
		gd.addNumericField("X_Start",xstart,5,15,null);
		gd.addNumericField("X_End",xend,5,15,null);
		gd.addMessage("X Values must be evenly spaced");
		gd.addMessage("Out of Range values will be set to zero");
		gd.showDialog(); if(gd.wasCanceled()){return;}
		xinc=(float)gd.getNextNumber();
		xstart=(float)gd.getNextNumber();
		xend=(float)gd.getNextNumber();
		int newlength=1+(int)((xend-xstart)/xinc);
		float[] newxvals=new float[newlength];
		for(int i=0;i<newlength;i++){
			newxvals[i]=xstart+i*xinc;
		}
		PlotWindow4 pw=null;
		for(int i=0;i<xvals.length;i++){
			float tempxinc=xvals[i][1]-xvals[i][0];
			float scaledxinc=xinc/tempxinc;
			float scaledstart=(xstart-xvals[i][0])/tempxinc;
			//float scaledend=(float)npts[i]-1.0f+(xend-xvals[i][npts[i]-1])/tempxinc;
			float[] newtraj=resample_trajectory(scaledstart,scaledxinc,newlength,yvals[i],npts[i]);
			if(i==0){
				pw=new PlotWindow4(iw.getTitle(),xlabel,ylabel,newxvals,newtraj);
				pw.draw();
			} else {
				pw.addPoints(newxvals,newtraj,true);
			}
			IJ.showProgress(i,xvals.length);
		}
	}

	/*public float[] resample_trajectory(float start,float end,float inc,float[] traj,int trajlen){
		int npts=1+(int)((end-start)/inc);
		float[] interp=new float[npts];
		for(int i=0;i<npts;i++){
			float x=start+i*inc;
			interp[i]=interpolation.interp1D(traj,trajlen,x);
		}
		return interp;
	}*/

	public float[] resample_trajectory(float start,float inc,int npts,float[] traj,int trajlen){
		float[] interp=new float[npts];
		for(int i=0;i<npts;i++){
			float x=start+i*inc;
			interp[i]=interpolation.interp1D(traj,trajlen,x);
		}
		return interp;
	}

}
