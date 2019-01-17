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
import javax.script.*;
import jguis.*;
import ij.text.*;

public class fit_custom_function_jru_v1 implements PlugIn, NLLSfitinterface_v2, DialogListener {
	float[] tempx,tempdata,weights,errs;
	double[][] constraints;
	String function,exdef,weightfunction;
	boolean checkc2,redirect;
	int iterations,series;
	int hitcounter;
	float c2;
	//ScriptEngineManager manager;
	ScriptEngine engine;
	Compilable ce;
	CompiledScript cs;
	PlotWindow4 pw;
	TextWindow tw;

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		pw=jutils.getPW4SelCopy(iw);
		String title=pw.getTitle();
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		int length=yvals[0].length;
		if(pw.getShowErrors()) errs=pw.getErrors(0,false);
		int[] colors=pw.getColors();
		colors[0]=0;
		ScriptEngineManager manager=new ScriptEngineManager();
		engine=manager.getEngineByName("js");
		ce=(Compilable)engine;
		//hitcounter=0;

		c2=0.0f;
		iterations=0;
		checkc2=false;

		double[] stats=new double[3];
		tempx=new float[length]; tempdata=new float[length];
		System.arraycopy(xvals[0],0,tempx,0,length);
		System.arraycopy(yvals[0],0,tempdata,0,length);
		pw.addPoints(tempx,new float[tempx.length],false);
		series=pw.getNpts().length-1;
		double[] params=new double[10];
		int[] fixes={0,0,0,1,1,1,1,1,1,1};
		init_options(params,fixes);
		if(!init_functions()){return;}
		
		while(showoptions(params,fixes)){
			NLLSfit_v2 fitclass;
			if(checkc2){
				fitclass=new NLLSfit_v2(this,0);
			} else {
				fitclass=new NLLSfit_v2(this,0.0001,50,0.1);
			}
			float[] fit=fitclass.fitdata(params,fixes,constraints,yvals[0],weights,stats,true);
			pw.updateSeries(fit,series,false);
			c2=(float)stats[1];
			iterations=(int)stats[0];
		}

		IJ.log("Chi Squared = "+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		for(int i=0;i<10;i++){
			IJ.log("P"+(i+1)+" = "+(float)params[i]+" fixed = "+fixes[i]);
		}
		IJ.log("AIC = "+(float)stats[2]);
		//IJ.log("hits = "+hitcounter);
		set_options(params,fixes);
	}

	boolean init_functions(){
		GenericDialog gd=new GenericDialog("Fitting Options");
		gd.addStringField("Extra Definitions",exdef,50);
		gd.addCheckbox("Weight Using Plot Errors",false);
		gd.addStringField("Weighting Equation (y is for data)",weightfunction,50);
		gd.addStringField("Fit_Equation",function,50);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		exdef=gd.getNextString();
		boolean errweights=gd.getNextBoolean();
		weightfunction=gd.getNextString();
		function=gd.getNextString();
		//first initialize the weights
		weights=new float[tempdata.length];
		if(errweights || weightfunction.equals("") || weightfunction==null || weightfunction=="1.0"){
			if(errweights){
				for(int i=0;i<tempdata.length;i++) weights[i]=1.0f/(errs[i]*errs[i]);
			} else {
				for(int i=0;i<tempdata.length;i++) weights[i]=1.0f;
			}
		} else {
			for(int i=0;i<tempdata.length;i++){
				String script="y ="+tempdata[i]+"; "+
				"x ="+tempx[i]+"; "+
				"retval="+weightfunction+";";
				Double temp=new Double(0.0);
				try{
					temp=(Double)engine.eval(script);
				}catch(Exception e){
					IJ.log(e.getMessage());
				}
				if(!(temp.isInfinite() || temp.isNaN())){
					weights[i]=temp.floatValue();
				}
			}
		}
		//now compile the function script
		try{
			String script1=exdef+"; retval="+function+";";
			cs=ce.compile(script1);
		}catch(Exception e){
			IJ.log(e.toString());
			return false;
		}
		return true;
	}

	boolean showoptions(double[] params,int[] fixes){
		//GenericDialog gd=new NonBlockingGenericDialog("Options");
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Check Chi Squared",checkc2);
		for(int i=0;i<10;i++){
			gd.addNumericField("P"+(i+1),params[i],5,10,null);
			gd.addCheckbox("Fix?",(fixes[i]==1));
		}
		gd.addCheckbox("Get_Errors",false);
		gd.addCheckbox("Set_Constraints",false);
		gd.addNumericField("Iterations",iterations,0,10,null);
		gd.addNumericField("chi squared",c2,5,10,null);
		gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		checkc2=gd.getNextBoolean();
		for(int i=0;i<10;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){fixes[i]=1;}
			else{fixes[i]=0;}
		}
		boolean geterrors=gd.getNextBoolean();
		boolean setconstraints=gd.getNextBoolean();
		for(int i=0;i<10;i++){
			if(function.indexOf("P"+(i+1))<0){fixes[i]=1;}
		}
		if(geterrors){
			if(!get_errors(params,fixes)){return false;}
		}
		if(setconstraints) constraints=get_constraints(params);
		return true;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
		checkc2=gd.getNextBoolean();
		double[] params=new double[10];
		int[] fixes=new int[10];
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){fixes[i]=1;}
			else{fixes[i]=0;}
		}
		gd.getNextBoolean();
		gd.getNextNumber();
		gd.getNextNumber();
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,tempdata,null,stats,true);
		pw.updateSeries(fit,series,false);
		c2=(float)stats[1];
		return true;
	}

	public boolean get_errors(double[] params,int[] fixes){
		GenericDialog gd=new GenericDialog("Error Options");
		String[] methods={"Support Plane","Monte Carlo"};
		gd.addChoice("Method",methods,methods[0]);
		float conf=0.67f;
		gd.addNumericField("SP_Confidence Limit (%)",(int)(conf*100.0f),5,10,null);
		String[] labels={"P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};
		gd.addChoice("SP_Parameter",labels,labels[0]);
		double spacing=0.01;
		gd.addNumericField("SP_Chi^2_plot_spacing (% of value)?",spacing*100.0,2,10,null);
		int ntrials=100;
		gd.addNumericField("MC_#_Trials",ntrials,0);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		int methodindex=gd.getNextChoiceIndex();
		conf=0.01f*(float)gd.getNextNumber();
		int paramindex=gd.getNextChoiceIndex();
		spacing=0.01*gd.getNextNumber();
		ntrials=(int)gd.getNextNumber();
		if(methodindex==0){
			support_plane_errors_v2 erclass=new support_plane_errors_v2(this,0.0001,50,false,0.1);
			int errindex=paramindex;
			int nfit=0;
			for(int i=0;i<labels.length;i++){
				if(fixes[i]==0){
					nfit++;
				}
			}
			int npts=tempdata.length;
			int dofnum=npts-(nfit-1)-1;
			int dofden=npts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,(double)conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN && flim < 1.0){IJ.showMessage("Invalid Limiting F Value"); return false;}
			double truespacing=Math.abs(params[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(params,fixes,constraints,tempdata,weights,flim,truespacing,errindex);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			IJ.log("upper error = "+(c2plot[1][0]-params[errindex])+" lower error = "+(params[errindex]-c2plot[0][0]));
			int templength=c2plot[0].length;
			float[][] c2plotf=new float[2][templength-1];
			for(int i=0;i<(templength-1);i++){
				c2plotf[0][i]=(float)c2plot[0][i+1];
				c2plotf[1][i]=(float)c2plot[1][i+1];
			}
			new PlotWindow4("c2 plot",labels[errindex],"Chi^2",c2plotf[0],c2plotf[1]).draw();
		} else {
			StringBuffer sb=new StringBuffer();
			sb.append("Trial\t");
			for(int i=0;i<labels.length;i++){
				if(fixes[i]==0) sb.append(labels[i]+"\t");
			}
			sb.append("chi^2");
			tw=new TextWindow("Monte Carlo Results",sb.toString(),"",400,400);
			redirect=true;
			monte_carlo_errors_v2 erclass=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
			double[][] errors=erclass.geterrors(params,fixes,constraints,tempdata,weights,ntrials);
			sb=new StringBuffer();
			sb.append("StDev\t");
			for(int i=0;i<errors.length;i++){
				float[] ferr=new float[errors[0].length];
				for(int j=0;j<ferr.length;j++) ferr[j]=(float)errors[i][j];
				float stdev=jstatistics.getstatistic("StDev",ferr,null);
				sb.append(""+stdev);
				if(i<(errors.length-1)) sb.append("\t");
			}
			tw.append(sb.toString());
			redirect=false;
		}
		return true;
	}

	public double[][] get_constraints(double[] params){
		//here we populate the constraints
		GenericDialog gd=new GenericDialog("Constraints");
		for(int i=0;i<10;i++){
			if(constraints==null){
				gd.addNumericField("P"+(i+1)+"_upper",params[i],5,10,null);
				gd.addNumericField("P"+(i+1)+"_lower",params[i],5,10,null);
			} else {
				gd.addNumericField("P"+(i+1)+"_upper",constraints[1][i],5,10,null);
				gd.addNumericField("P"+(i+1)+"_lower",constraints[0][i],5,10,null);
			}
		}
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		double[][] constraints=new double[2][10];
		for(int i=0;i<10;i++){
			constraints[1][i]=gd.getNextNumber();
			constraints[0][i]=gd.getNextNumber();
		}
		return constraints;
	}

	public double[] fitfunc(double[] fitparams){
		Bindings b=engine.createBindings();
		for(int i=0;i<10;i++) b.put("P"+(i+1),fitparams[i]);
		/*String script1="P1="+fitparams[0]+"; "+
		"P2="+fitparams[1]+"; "+
		"P3="+fitparams[2]+"; "+
		"P4="+fitparams[3]+"; "+
		"P5="+fitparams[4]+"; "+
		"P6="+fitparams[5]+"; "+
		"P7="+fitparams[6]+"; "+
		"P8="+fitparams[7]+"; "+
		"P9="+fitparams[8]+"; "+
		"P10="+fitparams[9]+"; "+
		exdef+"; x=";
		String script2="; retval="+function+";";*/
		try{
			double[] temp=new double[tempx.length];
			for(int i=0;i<tempx.length;i++){
				//temp[i]=((Double)engine.eval(script1+(double)tempx[i]+script2)).doubleValue();
				b.put("x",tempx[i]);
				b.put("y",tempdata[i]);
				temp[i]=(Double)cs.eval(b);
			}
			return temp;
		}catch(Exception e){
			IJ.log(e.getMessage());
			return null;
		}
	}

	public void showresults(String results){
		if(redirect) tw.append(results);
		else IJ.log(results);
	}

	public void init_options(double[] params,int[] fixes){
		String[] initvals=jutils.get_plugin_options("fit_custom_function_jru_v1",23+21);
		if(initvals==null){
			function="P1*Math.exp(-x/P2)+P3";
			for(int i=0;i<10;i++) params[i]=0.0;
			for(int i=0;i<10;i++) fixes[i]=1;
			fixes[0]=0; fixes[1]=0; fixes[2]=0;
			weightfunction="1.0/y";
			exdef="";
			constraints=null;
		} else {
			function=initvals[0];
			int counter=1;
			for(int i=0;i<10;i++){
				params[i]=Double.parseDouble(initvals[counter]);
				counter++;
				fixes[i]=Integer.parseInt(initvals[counter]);
				counter++;
			}
			weightfunction=initvals[counter]; counter++;
			exdef=initvals[counter]; counter++;
			int hasconstraints=Integer.parseInt(initvals[counter]); counter++;
			if(hasconstraints==1){
				constraints=new double[2][10];
				for(int i=0;i<10;i++){
					constraints[0][i]=Double.parseDouble(initvals[counter]); counter++;
					constraints[1][i]=Double.parseDouble(initvals[counter]); counter++;
				}
			} else {
				constraints=null;
			}
		}
	}

	public void set_options(double[] params,int[] fixes){
		String[] initvals=new String[23+21];
		initvals[0]=function;
		int counter=1;
		for(int i=0;i<10;i++){
			initvals[counter]=""+params[i];
			counter++;
			initvals[counter]=""+fixes[i];
			counter++;
		}
		initvals[counter]=weightfunction; counter++;
		initvals[counter]=exdef; counter++;
		if(constraints!=null){
			initvals[counter]="1"; counter++;
			for(int i=0;i<10;i++){
				initvals[counter]=""+constraints[0][i]; counter++;
				initvals[counter]=""+constraints[1][i]; counter++;
			}
		} else {
			initvals[counter]="0"; counter++;
			for(int i=0;i<10;i++){
				initvals[counter]="0"; counter++;
				initvals[counter]="0"; counter++;
			}
		}
		jutils.set_plugin_options("fit_custom_function_jru_v1",initvals);
	}

}
