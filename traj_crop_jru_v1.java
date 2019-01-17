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

public class traj_crop_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4Copy(iw);
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		float[][][] errs=pw.getErrors();
		int[] npts=pw.getNpts();
		int length=npts[0];
		GenericDialog gd=new GenericDialog("Options");
		float[] limits=pw.getLimits();
		float xmin=limits[0];
		float xmax=limits[1];
		gd.addNumericField("Start Value",xmin,5,10,null);
		gd.addNumericField("End Value",xmax,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		xmin=(float)gd.getNextNumber();
		xmax=(float)gd.getNextNumber();
		int counter1=-1;
		do{
			counter1++;
		}while(xvals[0][counter1]<xmin && counter1<(length-1));
		int counter2=counter1;
		do{
			counter2++;
		}while(xvals[0][counter2]<xmax && counter2<(length-1));
		int newlength=counter2-counter1+1;
		float[][] newxvals=new float[xvals.length][newlength];
		float[][] newyvals=new float[xvals.length][newlength];
		float[][][] newerrs=null;
		if(errs!=null) newerrs=new float[2][xvals.length][newlength];
		for(int j=0;j<xvals.length;j++){
			if((counter1+newlength)>npts[j]){
				newxvals[j]=new float[npts[j]-counter1];
				newyvals[j]=new float[npts[j]-counter1];
				if(errs!=null){
					newerrs[0][j]=new float[npts[j]-counter1];
					newerrs[1][j]=new float[npts[j]-counter1];
				}
			}
			for(int i=0;i<newlength;i++){
				if((counter1+i)<npts[j]){
					newxvals[j][i]=xvals[j][counter1+i];
					newyvals[j][i]=yvals[j][counter1+i];
					if(errs!=null){
						newerrs[0][j][i]=errs[0][j][counter1+i];
						newerrs[1][j][i]=errs[1][j][counter1+i];
					}
				}
			}
			pw.updateSeries(newxvals[j],newyvals[j],j,true);
		}
		pw.addErrors(newerrs);
	}

}
