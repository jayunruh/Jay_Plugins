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
import jalgs.jsim.*;
import jguis.*;

public class propagate_errors_jru_v1 implements PlugIn {
	String function;
	int iterations;
	double stdev, uplim, lowlim;
	boolean hist;

	public void run(String arg) {
		iterations=100000;
		lowlim=0.0;
		uplim=100000.0;
		function="pow(2.0,-P1)";
		stdev=0.0;
		double[] params=new double[5];
		double[] stdevs=new double[5];
		if(!showoptions(params,stdevs)){
			return;
		}
		int numfit=0;
		for(int i=0;i<5;i++){
			if(stdevs[i]>0.0){numfit++;}
		}
		//IJ.log("uplim = "+uplim+" , "+lowlim);
		rngs random=new rngs();
		double val=fitfunc(params);
		double avg=0.0;
		double avgsq=0.0;
		float[] traj=new float[iterations];
		for(int i=0;i<iterations;i++){
			double[] params2=new double[5];
			for(int j=0;j<numfit;j++){
				params2[j]=random.gasdev(params[j],stdevs[j]);
			}
			double temp=fitfunc(params2);
			if(temp>uplim) temp=uplim;
			if(temp<lowlim) temp=lowlim;
			avg+=temp/(double)iterations;
			avgsq+=(temp*temp)/(double)iterations;
			traj[i]=(float)temp;
			IJ.showProgress(i,iterations);
		}
		if(hist){(new PlotWindowHist("Histogram","function","frequency",traj,3)).draw();}
		float[] quartiles={25.0f,50.0f,75.0f};
		jstatistics.fpercentile(traj,quartiles);
		stdev=Math.sqrt(avgsq-avg*avg);
		IJ.log("Avg = "+val+" St. Dev. = "+stdev);
		IJ.log("Quartiles = "+quartiles[0]+" , "+quartiles[1]+" , "+quartiles[2]);
	}

	boolean showoptions(double[] params,double[] stdevs){
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<5;i++){
			gd.addNumericField("P"+(i+1),params[i],5,10,null);
			gd.addNumericField("Stdev"+(i+1),stdevs[i],5,10,null);
		}
		gd.addStringField("Equation",function,50);
		gd.addNumericField("Iterations",iterations,0,10,null);
		gd.addCheckbox("Histogram",hist);
		gd.addNumericField("Lower_Limit",0.0,5,15,null);
		gd.addNumericField("Upper_Limit",100000.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		for(int i=0;i<5;i++){
			params[i]=gd.getNextNumber();
			stdevs[i]=gd.getNextNumber();
		}
		function=gd.getNextString();
		iterations=(int)gd.getNextNumber();
		hist=gd.getNextBoolean();
		lowlim=gd.getNextNumber();
		uplim=gd.getNextNumber();
		return true;
	}

	public double fitfunc(double[] fitparams){
		String semicol=";";
		String macro="P1 = "+fitparams[0]+semicol+"\n"+
		"P2 = "+fitparams[1]+semicol+"\n"+
		"P3 = "+fitparams[2]+semicol+"\n"+
		"P4 = "+fitparams[3]+semicol+"\n"+
		"P5 = "+fitparams[4]+semicol+"\n"+
		"retval = "+function+semicol+"\n"+
		"return \"\"+retval"+semicol;
		//IJ.log(macro);
		Macro_Runner mr=new Macro_Runner();
		String retval=mr.runMacro(macro,null);
		return (double)Float.parseFloat(retval);
	}

}
