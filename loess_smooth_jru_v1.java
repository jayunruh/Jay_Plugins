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
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import jalgs.jfit.*;

public class loess_smooth_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Span (odd)",11,0);
		String[] types={"linear","quadratic"};
		gd.addChoice("Type",types,types[1]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int span=(int)gd.getNextNumber();
		int type=gd.getNextChoiceIndex();
		int order=type+1;
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4SelCopy(iw);
		float[][] yvals=pw.getYValues();
		fit_loess fl=new fit_loess(yvals[0].length,span,order);
		float[] fit=fl.getfit(yvals[0]);
		float[] xvals=pw.getXValues()[0];
		pw.addPoints(xvals,fit,true);
	}

}
