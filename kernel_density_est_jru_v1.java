/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jfit.*;

public class kernel_density_est_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		if(jutils.isPW4(iw)) xvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=null;
		if(!iw.getClass().getName().equals("jguis.PlotWindowHist")) npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		String xlab=(String)jutils.runPW4VoidMethod(iw,"getxLabel");
		String ylab=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		float[] data=xvals[0];
		if(npts!=null && npts[0]!=data.length){
			float[] temp=(float[])algutils.get_subarray(data,0,npts[0]);
			data=temp;
		}
		//data=new float[]{1.0f,2.0f,3.0f};
		//for(int i=0;i<data.length;i++) IJ.log(""+data[i]);
		float bandwidth=fit_kde_gaus.get_kde_stdev(data);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Bandwidth",bandwidth,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		bandwidth=(float)gd.getNextNumber();
		//bandwidth=0.1f;
		IJ.log("bandwidth = "+bandwidth);
		float[][] kde=fit_kde_gaus.get_kde(data,bandwidth);
		PlotWindow4 pw=new PlotWindow4("KDE_Gaussian",xlab,ylab,kde[0],kde[1]);
		pw.draw();
		/*float[] tempx=pw.getXValues(0);
		float start=tempx[0];
		float dx=tempx[1]-tempx[0];
		int npts1=pw.getNpts()[0];
		float[][] allgaus=new float[data.length][];
		float[][] xvals2=new float[data.length][];
		gausfunc gf=new gausfunc();
		for(int i=0;i<data.length;i++){
			xvals2[i]=kde[0];
			allgaus[i]=gf.get_norm_func(2.0f*start-data[i],npts1,dx,bandwidth);
		}
		new PlotWindow4("All Kernels","x","y",xvals2,allgaus,null).draw();*/
	}

}
