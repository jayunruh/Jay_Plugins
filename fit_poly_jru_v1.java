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
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;

public class fit_poly_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin fits a trajectory to a polynomial line of specified order
		GenericDialog gd=new GenericDialog("Options");
		int order=2;
		gd.addNumericField("Polynomial Order?",order,0);
		boolean fixoff=false;
		gd.addCheckbox("Fix offset?",fixoff);
		float offset=0.0f;
		gd.addNumericField("Offset (if fixed)",offset,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		order=(int)gd.getNextNumber();
		fixoff=gd.getNextBoolean();
		offset=(float)gd.getNextNumber();

		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4SelCopy(iw);
		float[] xvals=pw.getXValues(0);
		float[] yvals=pw.getYValues(0);
		int[] colors=pw.getColors();
		colors[0]=0;

		float[] tempy=new float[yvals.length];
		if(fixoff){
			for(int i=0;i<yvals.length;i++){tempy[i]=yvals[i]-offset;}
		} else {
			for(int i=0;i<yvals.length;i++){tempy[i]=yvals[i];}
		}
		fitpoly fitclass=new fitpoly(order,xvals,fixoff);
		float[] params=fitclass.fitdatapoly(tempy,null);
		float[] fit=fitclass.getfit(params);
		for(int i=0;i<params.length;i++){
			if(fixoff){
				IJ.log("Coefficient "+(i+1)+" = "+params[i]);
			} else {
				IJ.log("Coefficient "+i+" = "+params[i]);
			}
		}
		pw.addPoints(xvals,fit,false);
	}

}
