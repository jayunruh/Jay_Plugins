/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.text.TextPanel;
import ij.text.TextWindow;
import jalgs.jdist;
import jalgs.jstatistics;
import jalgs.jfit.NLLSfit_v2;
import jalgs.jfit.NLLSfitinterface_v2;
import jalgs.jfit.monte_carlo_errors_v2;
import jalgs.jfit.support_plane_errors_v2;

import javax.swing.JTable;
import javax.swing.event.TableModelEvent;

public class FitDialog_v2 implements NLLSfitinterface_v2,TableDialogListener{
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
	public int[] fixes;
	public float c2;
	private TextWindow tw;
	private boolean redirect;

	// this class handles the procedures associated with fitting
	// this version uses a table dialog for more fitting parameters
	public FitDialog_v2(PlotWindow4 pw,NLLSfitinterface_v2 callclass,String[] labels){
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
		xvals=pw.getXValues(0);
		yvals=pw.getYValues(0);
		pw.addPoints(xvals,new float[xvals.length],false);
		redirect=false;
	}

	public void run_fit(double[] params,float[] weights,double[][] constraints,int[] fixes){
		this.weights=weights;
		this.params=params;
		this.constraints=constraints;
		this.fixes=fixes;
		double[] stats=new double[2];
		while(showoptions(this.params,this.fixes)){
			if(checkc2){
				fitclass.maxiter=0;
			}else{
				fitclass.maxiter=50;
			}
			float[] fit=fitclass.fitdata(this.params,this.fixes,this.constraints,yvals,weights,stats,true);
			c2=(float)stats[1];
			pw.updateSeries(fit,1,true);
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
		while(showoptions(this.params,this.fixes)){
			if(checkc2){
				fitclass.maxiter=0;
			}else{
				fitclass.maxiter=50;
			}
			float[] fit=fitclass.fitdata(this.params,this.fixes,this.constraints,yvals,weights,stats,true);
			c2=(float)stats[1];
			pw.updateSeries(fit,1,true);
			iterations=(int)stats[0];
			if(!manual){
				break;
			}
			IJ.selectWindow(pw.getTitle());
		}

		StringBuffer sb=new StringBuffer();
		sb.append(pw.getTitle());
		IJ.log("Chi Squared = "+(float)stats[1]);
		sb.append("\t"+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		sb.append("\t"+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(labels[i]+" : "+params[i]+" : Fixed : "+(fixes[i]==1));
			sb.append("\t"+(float)params[i]);
		}
		if(outtable!=null)
			outtable.append(sb.toString());
	}

	public TextWindow make_outtable(String title){
		String collabels="title\tc^2\tIter\t"+table_tools.print_string_array(labels);
		TextWindow outtable=new TextWindow(title,collabels,"",400,200);
		return outtable;
	}

	public static TextWindow make_outtable(String title,String[] labels1){
		String collabels="title\tc^2\tIter\t"+table_tools.print_string_array(labels1);
		TextWindow outtable=new TextWindow(title,collabels,"",400,200);
		return outtable;
	}
	
	/***************
	 * here we add columns only if the new fit has more parameters
	 * @param tw
	 * @param labels1
	 */
	public static void adapt_outtable(TextWindow tw,String[] labels1){
		TextPanel tp=tw.getTextPanel();
		String[] oldlabels=table_tools.getcollabels(tp);
		if(labels1.length>(oldlabels.length-3)){
			String collabels="title\tc^2\tIter\t"+table_tools.print_string_array(labels1);
			table_tools.change_table_labels(tp,collabels);
		}
	}

	public boolean showoptions(double[] params,int[] fixes){
		updatePlot(params);
		int nparams=params.length;
		Object[][] tabledata=new Object[nparams+5][3];
		String[] columnlabels={"Parameters","Values","Fix?"};
		tabledata[0][0]="Check chi^2?";
		tabledata[0][1]=new Boolean(checkc2);
		for(int i=0;i<nparams;i++){
			tabledata[i+1][0]=labels[i];
			tabledata[i+1][1]=new Double(params[i]);
			tabledata[i+1][2]=new Boolean(fixes[i]==1);
		}
		tabledata[nparams+1][0]="Get Errors";
		tabledata[nparams+1][1]=new Boolean(false);
		tabledata[nparams+2][0]="Edit Constraints";
		tabledata[nparams+2][1]=new Boolean(false);
		tabledata[nparams+3][0]="Iterations";
		tabledata[nparams+3][1]=new Integer(iterations);
		tabledata[nparams+4][0]="Chi Squared";
		tabledata[nparams+4][1]=new Double(c2);
		TableDialog2 td=new TableDialog2(null,null,"Fit Parameters",columnlabels,tabledata,null);
		td.addTableDialogListener(this);
		Object[][] retvals=TableDialog2.showDialog(td);
		//Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Fit Parameters",columnlabels,tabledata,null);
		if(retvals==null) return false;
		checkc2=((Boolean)retvals[0][1]).booleanValue();
		for(int i=0;i<nparams;i++){
			params[i]=((Double)retvals[i+1][1]).doubleValue();
			fixes[i]=((Boolean)retvals[i+1][2]).booleanValue()?1:0;
		}
		boolean geterrs=((Boolean)retvals[nparams+1][1]).booleanValue();
		if(geterrs){
			if(!get_errors(params,fixes)) return false;
		}
		boolean showconstraints=((Boolean)retvals[nparams+2][1]).booleanValue();
		if(showconstraints){
			if(!showconstraintsdialog()) return false;
		}
		return true;
	}
	
	private boolean showconstraintsdialog(){
		int nparams=labels.length;
		Object[][] tabledata=new Object[nparams][3];
		String[] columnlabels={"Parameters","Lower Limit","Upper Limit"};
		for(int i=0;i<nparams;i++){
			tabledata[i][0]=labels[i];
			tabledata[i][1]=constraints[0][i];
			tabledata[i][2]=constraints[1][i];
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Constraints",columnlabels,tabledata,null);
		if(retvals==null){
			return false;
		}
		for(int i=0;i<nparams;i++){
			constraints[0][i]=((Double)retvals[i][1]).doubleValue();
			constraints[1][i]=((Double)retvals[i][2]).doubleValue();
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
			for(int i=0;i<labels.length;i++){
				if(fixes[i]==0)
					sb.append(labels[i]+"\t");
			}
			sb.append("chi^2");
			tw=new TextWindow("Monte Carlo Results",sb.toString(),"",400,400);
			redirect=true;
			monte_carlo_errors_v2 erclass=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
			double[][] errors=erclass.geterrors(params,fixes,constraints,yvals,weights,ntrials);
			sb=new StringBuffer();
			sb.append("StDev\t");
			for(int i=0;i<errors.length;i++){
				float[] ferr=new float[errors[0].length];
				for(int j=0;j<ferr.length;j++)
					ferr[j]=(float)errors[i][j];
				float stdev=jstatistics.getstatistic("StDev",ferr,null);
				sb.append(""+stdev);
				if(i<(errors.length-1))
					sb.append("\t");
			}
			tw.append(sb.toString());
			redirect=false;
		}
		return true;
	}

	/*public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		checkc2=gd.getNextBoolean();
		double[] params=new double[labels.length];
		int[] fixes=new int[labels.length];
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){
				fixes[i]=1;
			}else{
				fixes[i]=0;
			}
		}
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,yvals,weights,stats,true);
		pw.updateSeries(fit,1,true);
		c2=(float)stats[1];
		return true;
	}*/

	public double[] fitfunc(double[] params){
		return callclass.fitfunc(params);
	}

	public void showresults(String results){
		if(redirect)
			tw.append(results);
		else
			callclass.showresults(results);
	}

	public void tableDataChanged(JTable table,TableModelEvent e,Object[][] tabledata){
		checkc2=((Boolean)tabledata[0][1]).booleanValue();
		double[] params=new double[labels.length];
		int[] fixes=new int[labels.length];
		for(int i=0;i<params.length;i++){
			params[i]=((Double)tabledata[i+1][1]).doubleValue();
			fixes[i]=((Boolean)tabledata[i+1][2]).booleanValue()?1:0;
		}
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,yvals,weights,stats,true);
		pw.updateSeries(fit,1,false);
		c2=(float)stats[1];
		return;
	}
	
	public void updatePlot(double[] params){
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,yvals,weights,stats,true);
		pw.updateSeries(fit,1,false);
	}
}
