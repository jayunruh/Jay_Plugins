/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.text.TextWindow;
import jalgs.jdist;
import jalgs.jstatistics;
import jalgs.jfit.NLLSfit_v2;
import jalgs.jfit.NLLSfitinterface_v2;
import jalgs.jfit.monte_carlo_errors_v2;
import jalgs.jfit.support_plane_errors_v2;

import java.awt.AWTEvent;
import java.awt.GridLayout;

public class FitDialog_v3 implements DialogListener,NLLSfitinterface_v2{
	public PlotWindow4 pw;
	public NLLSfitinterface_v2 callclass;
	public NLLSfit_v2 fitclass;
	public String[] labels;
	public boolean manual;
	public int iterations;
	public boolean checkc2;
	public float[] xvals,yvals,weights;
	public double[][] constraints;
	public double[] params;
	public float[] errs;
	public int[] fixes;
	public float c2;
	private TextWindow tw;
	private boolean redirect;

	/******************************
	 * this class handles the procedures associated with fitting a plot, this version uses a grid layout genericdialog
	 * @param pw
	 * @param callclass
	 * @param labels
	 */
	public FitDialog_v3(PlotWindow4 pw,NLLSfitinterface_v2 callclass,String[] labels){
		this.pw=pw;
		this.callclass=callclass;
		this.labels=labels;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Manual_Fit",true);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		manual=gd.getNextBoolean();
		fitclass=new NLLSfit_v2(this,0.0001,50,0.1);
		if(pw!=null){
			xvals=pw.getXValues(0);
			yvals=pw.getYValues(0);
			pw.addPoints(xvals,new float[xvals.length],true);
		}
		redirect=false;
	}
	
	public FitDialog_v3(float[] ypts,NLLSfitinterface_v2 callclass,String[] labels){
		//this.pw=pw;
		this.callclass=callclass;
		this.labels=labels;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Manual_Fit",true);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		manual=gd.getNextBoolean();
		fitclass=new NLLSfit_v2(this,0.0001,50,0.1);
		this.yvals=ypts;
		redirect=false;
	}

	public void run_fit(double[] params,float[] weights,double[][] constraints,int[] fixes){
		this.weights=weights;
		this.params=params;
		this.constraints=constraints;
		this.fixes=fixes;
		double[] stats=new double[2];
		while(showoptions(params,fixes)){
			if(checkc2){
				fitclass.maxiter=0;
			}else{
				fitclass.maxiter=50;
			}
			float[] fit=fitclass.fitdata(params,fixes,constraints,yvals,weights,stats,true);
			c2=(float)stats[1];
			if(pw!=null) pw.updateSeries(fit,1,true);
			iterations=(int)stats[0];
			if(!manual){
				break;
			}
		}

		IJ.log("Chi Squared = "+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(labels[i]+" : "+params[i]+" : Fixed : "+(fixes[i]==1));
		}
	}

	public void run_fit(double[] params,float[] weights,double[][] constraints,int[] fixes,TextWindow outtable){
		this.weights=weights;
		this.params=params;
		this.constraints=constraints;
		this.fixes=fixes;
		double[] stats=new double[2];
		while(showoptions(params,fixes)){
			if(checkc2){
				fitclass.maxiter=0;
			}else{
				fitclass.maxiter=50;
			}
			float[] fit=fitclass.fitdata(params,fixes,constraints,yvals,weights,stats,true);
			c2=(float)stats[1];
			if(pw!=null) pw.updateSeries(fit,1,true);
			iterations=(int)stats[0];
			if(!manual){
				break;
			}
		}

		StringBuffer sb=new StringBuffer();
		if(pw!=null) sb.append(pw.getTitle());
		else sb.append("no plot");
		IJ.log("Chi Squared = "+(float)stats[1]);
		sb.append("\t"+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		sb.append("\t"+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(labels[i]+" : "+params[i]+" : Fixed : "+(fixes[i]==1));
			sb.append("\t"+(float)params[i]);
		}
		if(outtable!=null){
			//IJ.log(sb.toString());
			outtable.append(sb.toString());
		}
	}

	public TextWindow make_outtable(String title){
		String collabels="title\tc^2\tIter\t"+table_tools.print_string_array(labels);
		TextWindow outtable=new TextWindow(title,collabels,"",400,200);
		return outtable;
	}
	
	public void append_outtable_params(String title,String label,double[] params) {
		//this assumes the outtable has been created
		TextWindow tw=jutils.selectTable(title);
		tw.append(label+"\t"+c2+"\t"+iterations+"\t"+table_tools.print_double_array(params));
	}
	
	public void append_outtable_errs(String title,String label) {
		//this assumes the outtable has been created and monte carlo errors have been acquired
		//note that the errs array has nparams+1 values with the last one being the chi squared stdev
		if(errs!=null) {
			TextWindow tw=jutils.selectTable(title);
			StringBuffer sb=new StringBuffer();
			sb.append(label+"_errs"+"\t"+errs[labels.length]+"\t"+0);
			for(int i=0;i<labels.length;i++) {
				sb.append("\t"+errs[i]);
			}
			tw.append(sb.toString());
		}
	}

	public static TextWindow make_outtable(String title,String[] labels1){
		String collabels="title\tc^2\tIter\t"+table_tools.print_string_array(labels1);
		TextWindow outtable=new TextWindow(title,collabels,"",400,200);
		return outtable;
	}

	public boolean showoptions(double[] params,int[] fixes){
		GenericDialog gd=new GenericDialog("Starting Fit Parameters");
		gd.setLayout(new GridLayout(params.length+5, 2));
		gd.addCheckbox("Check Chi Squared",checkc2); gd.addMessage(" "); gd.addMessage(" ");
		for(int i=0;i<params.length;i++){
			gd.addNumericField(labels[i],params[i],5,10,null); gd.addCheckbox("Fix"+(i+1)+"?",(fixes[i]==1));
		}
		gd.addCheckbox("Get_Errors",false); gd.addCheckbox("Edit_Constraints",false); gd.addMessage(" ");
		gd.addNumericField("Iterations",iterations,0); gd.addMessage(" ");
		gd.addNumericField("Chi Squared",c2,5,10,null); gd.addMessage(" ");
		gd.addDialogListener(this);
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
		checkc2=gd.getNextBoolean();
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){
				fixes[i]=1;
			}else{
				fixes[i]=0;
			}
		}
		if(gd.getNextBoolean()){
			if(!get_errors(params,fixes)){
				return false;
			}
		}
		if(gd.getNextBoolean()){
			if(!showconstraintsdialog()){
				return false;
			}
		}
		return true;
	}
	
	private boolean showconstraintsdialog(){
		GenericDialog gd=new GenericDialog("Constraints");
		int nparams=labels.length;
		gd.setLayout(new GridLayout(nparams+2,4));
		gd.addMessage(" "); gd.addMessage("Lower_Limit"); gd.addMessage(" "); gd.addMessage("Upper_Limit"); 
		for(int i=0;i<nparams;i++){
			gd.addNumericField(labels[i],constraints[0][i],5,15,null); gd.addNumericField("upper",constraints[1][i],5,15,null); 
		}
		gd.showDialog(); if(gd.wasCanceled()) return false;
		for(int i=0;i<nparams;i++){
			constraints[0][i]=gd.getNextNumber();
			constraints[1][i]=gd.getNextNumber();
		}
		return true;
	}

	public boolean get_errors(double[] params,int[] fixes){
		GenericDialog gd=new GenericDialog("Error Options");
		String[] methods={"Support Plane","Monte Carlo"};
		gd.addChoice("Method",methods,methods[1]);
		float conf=0.67f;
		gd.addNumericField("SP_Confidence Limit (%)",(int)(conf*100.0f),5,10,null);
		gd.addChoice("SP_Parameter",labels,labels[0]);
		double spacing=0.01;
		gd.addNumericField("SP_Chi^2_plot_spacing (% of value)?",spacing*100.0,2,10,null);
		int ntrials=100;
		gd.addNumericField("MC_#_Trials",ntrials,0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
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
			int npts=yvals.length;
			int dofnum=npts-(nfit-1)-1;
			int dofden=npts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN&&flim<1.0){
				IJ.showMessage("Invalid Limiting F Value");
				return false;
			}
			double truespacing=Math.abs(params[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(params,fixes,constraints,yvals,weights,flim,truespacing,errindex);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			IJ.log("upper error = "+(c2plot[1][0]-params[errindex])+" lower error = "+(params[errindex]-c2plot[0][0]));
			int templength=c2plot[0].length;
			float[][] c2plotf=new float[2][templength-1];
			for(int i=0;i<(templength-1);i++){
				c2plotf[0][i]=(float)c2plot[0][i+1];
				c2plotf[1][i]=(float)c2plot[1][i+1];
			}
			new PlotWindow4("c2 plot",labels[errindex],"Chi^2",c2plotf[0],c2plotf[1]).draw();
		}else{
			StringBuffer sb=new StringBuffer();
			sb.append("Trial\t");
			int[] errindices=new int[labels.length+1];
			int cnt=0;
			for(int i=0;i<labels.length;i++){
				if(fixes[i]==0) {
					sb.append(labels[i]+"\t");
					errindices[cnt]=i;
					cnt++;
				}
			}
			sb.append("chi^2");
			errindices[cnt]=labels.length;
			tw=new TextWindow("Monte Carlo Results",sb.toString(),"",400,400);
			redirect=true;
			monte_carlo_errors_v2 erclass=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
			double[][] errors=erclass.geterrors(params,fixes,constraints,yvals,weights,ntrials);
			sb=new StringBuffer();
			sb.append("StDev\t");
			errs=new float[labels.length+1];
			for(int i=0;i<errors.length;i++){
				float[] ferr=new float[errors[0].length];
				for(int j=0;j<ferr.length;j++)
					ferr[j]=(float)errors[i][j];
				float stdev=jstatistics.getstatistic("StDev",ferr,null);
				errs[errindices[i]]=stdev;
				sb.append(""+stdev);
				if(i<(errors.length-1))
					sb.append("\t");
			}
			tw.append(sb.toString());
			redirect=false;
		}
		return true;
	}

	public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		checkc2=gd.getNextBoolean();
		double[] params=new double[labels.length];
		int[] fixes=new int[labels.length];
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()) fixes[i]=1;
			else fixes[i]=0;
		}
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,yvals,weights,stats,true);
		if(pw!=null) pw.updateSeries(fit,1,true);
		c2=(float)stats[1];
		return true;
	}

	public double[] fitfunc(double[] params){
		return callclass.fitfunc(params);
	}

	public void showresults(String results){
		if(redirect)
			tw.append(results);
		else
			callclass.showresults(results);
	}
}
