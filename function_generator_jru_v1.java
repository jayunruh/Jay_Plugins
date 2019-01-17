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
import jalgs.jsim.*;
import jguis.*;
import javax.script.*;

public class function_generator_jru_v1 implements PlugIn {
	String exdef,function;
	int noisetype;
	float gasstdev;
	rngs random;
	double xmin,xmax,dx;
	ScriptEngineManager manager;
	ScriptEngine engine;

	public void run(String arg) {
		//here we generate a function and add noise to it
		if(!def_function()){
			return;
		}
		manager=new ScriptEngineManager();
		engine=manager.getEngineByName("js");
		random=new rngs();
		int npts=1+(int)((xmax-xmin)/dx);
		float[] xvals=new float[npts];
		float[] yvals=new float[npts];
		for(int i=0;i<npts;i++){
			xvals[i]=(float)(xmin+(double)i*dx);
			yvals[i]=get_function(xvals[i]);
		}
		new PlotWindow4("Custom Function","x","y",xvals,yvals).draw();
	}

	boolean def_function(){
		GenericDialog gd=new GenericDialog("Function Options");
		gd.addStringField("Extra Definitions",exdef,50);
		gd.addStringField("Equation",function,50);
		String[] noisetypes={"None","Poisson","Gaussian"};
		gd.addChoice("Noise Type",noisetypes,noisetypes[0]);
		gd.addNumericField("Gaussian St Dev",1.0f,5,15,null);
		gd.addNumericField("X_min",0.0,5,15,null);
		gd.addNumericField("X_max",10.0,5,15,null);
		gd.addNumericField("Delta X",0.5,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		exdef=gd.getNextString();
		function=gd.getNextString();
		noisetype=gd.getNextChoiceIndex();
		gasstdev=(float)gd.getNextNumber();
		xmin=gd.getNextNumber();
		xmax=gd.getNextNumber();
		dx=gd.getNextNumber();
		if(function.startsWith("=")) function=function.substring(1);
		return true;
	}

	float get_function(double x){
		String script="x="+x+"; "+
		exdef+"; "+
		"retval="+function+";";
		Double temp=new Double(0.0);
		try{
			temp=(Double)engine.eval(script);
		}catch(Exception e){
			IJ.log(e.getMessage());
		}
		double val=temp.doubleValue();
		if(noisetype==0){
			return (float)val;
		} else {
			if(noisetype==1){
				return (float)random.poidev(val);
			} else {
				return (float)random.gasdev(val,(double)gasstdev);
			}
		}
	}

}
