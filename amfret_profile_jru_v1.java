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
import java.awt.Frame;
import java.awt.Polygon;
import ij.plugin.*;
import ij.plugin.frame.RoiManager;
import jguis.*;
import java.util.*;
import ij.text.*;

public class amfret_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin takes a set of window boundaries and returns fraction outside an roi gate for each of them
		//start by getting the histogram values and limits
		ImageWindow[] iw=jutils.selectPlotFamily(false,1);
		if(iw==null || iw.length<1) return;
		float[] xvals=(float[])jutils.runPW4VoidMethod(iw[0],"getXValues");
		float[] yvals=(float[])jutils.runPW4VoidMethod(iw[0],"getYValues");
		float[] limits=(float[])jutils.runPW4VoidMethod(iw[0],"getLimits");
		//get the gate roi
		RoiManager rman=RoiManager.getInstance();
		Roi roi=rman.getRoi(0);
		if(roi==null){IJ.error("Need gate roi in manager"); return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Bin_By",1,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int binby=(int)gd.getNextNumber();
		Object plot=jutils.runPW4VoidMethod(iw[0],"getPlot");
		FloatProcessor histimage=(FloatProcessor)jutils.runReflectionMethod(plot,"getHistImage",null);
		int histwidth=histimage.getWidth(); int histheight=histimage.getHeight();
		float[] histpix=(float[])histimage.getPixels();
		float[][] data=amfret_utils.getAmFRETProfile(histpix,histwidth,histheight,limits[2],limits[3],roi);
		float logxmin=(float)Math.log(limits[0]);
		float logxmax=(float)Math.log(limits[1]);
		float logxinc=(logxmax-logxmin)/(float)(histwidth-1);
		float[] xvals2=new float[histwidth];
		for(int i=0;i<histwidth;i++){
			float logval=logxmin+logxinc*(float)i;
			xvals2[i]=(float)Math.exp(logval);
		}
		float[][] populations={data[1],data[3]};
		for(int i=0;i<populations[0].length;i++){
			populations[0][i]/=16.0f; //correct for binning
			populations[1][i]/=16.0f; //correct for binning
		}
		float[][] amfretavgs={data[0],data[2]};
		//now optionally bin the data
		if(binby>1){
			float[][] newpops=new float[2][populations[0].length/binby];
			float[][] newamfretavgs=new float[2][populations[0].length/binby];
			float[] newxvals=new float[populations[0].length/binby];
			for(int i=0;i<populations[0].length;i+=binby){
				float amfretsum1=0.0f;
				float amfretsum2=0.0f;
				float popsum1=0.0f;
				float popsum2=0.0f;
				float xsum=0.0f;
				for(int j=0;j<binby;j++){
					amfretsum1+=amfretavgs[0][i+j]*populations[0][i+j];
					amfretsum2+=amfretavgs[1][i+j]*populations[1][i+j];
					popsum1+=populations[0][i+j];
					popsum2+=populations[1][i+j];
					xsum+=xvals2[i+j];
				}
				newpops[0][i/binby]=popsum1;
				newpops[1][i/binby]=popsum2;
				newamfretavgs[0][i/binby]=amfretsum1/popsum1;
				newamfretavgs[1][i/binby]=amfretsum2/popsum2;
				newxvals[i/binby]=xsum/(float)binby;
			}
			populations=newpops;
			amfretavgs=newamfretavgs;
			xvals2=newxvals;
		}
		float[][] xvals3={xvals2,xvals2};
		Plot4 popplot=new Plot4("Acceptor","Cell_Population",xvals3,populations,null);
		popplot.setLogAxes(true,false);
		new PlotWindow4("Populations",popplot).draw();

		Plot4 afplot=new Plot4("Acceptor","AmFRET_Avg",xvals3,amfretavgs,null);
		afplot.setLogAxes(true,false);
		new PlotWindow4("AmFRET_Avgs",afplot).draw();
	}

}
