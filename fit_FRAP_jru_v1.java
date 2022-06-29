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
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;

public class fit_FRAP_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	float xinc;
	float[] xvals2;
	int fitpts;

	public void run(String arg) {
		ImageWindow[] plots=jutils.selectPlots(false,1,new String[]{"FRAP_curve"});
		if(plots==null) return;
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(plots[0],"getYValues");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(plots[0],"getXValues");
		int sel=(Integer)jutils.runPW4VoidMethod(plots[0],"getSelected");
		int[] npts=(int[])jutils.runPW4VoidMethod(plots[0],"getNpts");
		if(sel<0) sel=0;
		xinc=xvals[sel][1]-xvals[sel][0];
		float[] weights=new float[npts[sel]];
		boolean haserrs=(Boolean)jutils.runPW4VoidMethod(plots[0],"getShowErrors");
		if(haserrs){
			//assume that error bars are sem values
			float[][][] errs=(float[][][])jutils.runPW4VoidMethod(plots[0],"getErrors");
			for(int i=0;i<npts[sel];i++) weights[i]=1.0f/(errs[0][0][i]*errs[0][0][i]);
		} else {
			//for(int i=0;i<npts[sel];i++) weights[i]=1.0f;
			weights=null;
		}
		int startpre=0;
		int endpre=3;
		int startrec=4;
		int endrec=npts[sel]-1;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Start_Pre_Bleach",startpre,0);
		gd.addNumericField("End_Pre_Bleach",endpre,0);
		gd.addNumericField("Start_Recovery",startrec,0);
		gd.addNumericField("End_Recovery",endrec,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		startpre=(int)gd.getNextNumber();
		endpre=(int)gd.getNextNumber();
		startrec=(int)gd.getNextNumber();
		endrec=(int)gd.getNextNumber();
		float[] pre=new float[endpre-startpre+1];
		for(int i=startpre;i<=endpre;i++) pre[i-startpre]=yvals[sel][i];
		fitpts=endrec-startrec+1;
		float[] rec=new float[fitpts];
		float[] weights2=null;
		if(weights!=null) weights2=new float[fitpts];
		float[][] errs2=null;
		if(weights!=null) errs2=new float[2][fitpts];
		xvals2=new float[fitpts];
		for(int i=startrec;i<=endrec;i++){
			rec[i-startrec]=yvals[sel][i];
			if(weights!=null) weights2[i-startrec]=weights[i];
			//xvals2[i-startrec]=xinc*(float)(i-startrec);
			xvals2[i-startrec]=xvals[sel][i]-xvals[sel][startrec];
			if(weights!=null) errs2[0][i-startrec]=1.0f/(float)Math.sqrt(weights[i]);
		}
		String[] labels={"Background","Amp1","Tau1","Amp2","Tau2"};
		float preavg=jstatistics.getstatistic("Avg",pre,null);
		float presem=jstatistics.getstatistic("StErr",pre,null);
		double[] params={rec[0],rec[rec.length-1]-rec[0],0.33f*(float)rec.length*xinc,0.0,(double)rec.length*xinc};
		double[][] constraints={{-10.0*Math.abs(params[0]),-10.0*Math.abs(params[1]),0.5*xinc,-10.0*Math.abs(params[1]),0.5*xinc},
					{10.0*Math.abs(params[0]),10.0*Math.abs(params[1]),10.0*xinc*(double)rec.length,10.0*Math.abs(params[1]),4.0*xinc*(double)rec.length}};
		int[] fixes={0,0,0,1,1};
		PlotWindow4 pw=new PlotWindow4("Recovery Fit","Time","Intensity",xvals2,rec);
		pw.draw();
		if(weights!=null) pw.addErrors(errs2);
		FitDialog fd=new FitDialog(pw,this,labels);
		fd.run_fit(params,weights2,constraints,fixes);
		IJ.log("Pre-FRAP = "+preavg+" SEM = "+presem);
	}

	public double[] fitfunc(double[] params){
		double[] fit=new double[fitpts];
		for(int i=0;i<fitpts;i++){
			//double t=(double)xinc*(double)i;
			double t=(double)xvals2[i];
			fit[i]=params[0]+params[1]+params[3]-params[1]*Math.exp(-t/params[2])-params[3]*Math.exp(-t/params[4]);
		}
		return fit;
	}

	public void showresults(String results){IJ.log(results);}

}
