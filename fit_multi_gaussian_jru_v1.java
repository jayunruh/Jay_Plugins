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
import ij.text.*;

public class fit_multi_gaussian_jru_v1 implements PlugIn, NLLSfitinterface, DialogListener {
	boolean checkc2;
	float[] tempx,tempdata;
	float c2;
	int iterations,series;
	NLLSfit fitclass;
	String[] paramnames={"baseline","xc1","stdev1","amp1","xc2","stdev2","amp2","xc3","stdev3","amp3","xc4","stdev4","amp4"};
	PlotWindow4 pw;

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		pw=jutils.getPW4SelCopy(iw);
		String title=pw.getTitle();
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		int[] colors=pw.getColors();
		colors[0]=0;
		int length=yvals[0].length;
		float min=yvals[0][0];
		int maxloc=0;
		float max=yvals[0][0];
		for(int i=1;i<length;i++){
			if(yvals[0][i]<min){min=yvals[0][i];}
			if(yvals[0][i]>max){max=yvals[0][i]; maxloc=i;}
		}
		float halfmax=(max+min)/2.0f;
		float fwhm=0.0f;
		for(int i=(maxloc+1);i<length;i++){
			if(yvals[0][i]<halfmax){
				fwhm+=0.5f*(xvals[0][i]-xvals[0][maxloc]);
				break;
			}
			if(i==(length-1)){
				fwhm+=0.5f*(xvals[0][i]-xvals[0][maxloc]);
			}
		}
		for(int i=(maxloc-1);i>=0;i--){
			if(yvals[0][i]<halfmax){
				fwhm+=0.5f*(xvals[0][maxloc]-xvals[0][i]);
				break;
			}
			if(i==0){
				fwhm+=0.5f*(xvals[0][maxloc]-xvals[0][i]);
			}
		}

		//parameters are baseline,xc1,stdev1,amp1,xc2,stdev2,amp2,xc3,stdev3,amp3,xc4,stdev4,amp4
		double[] params={(double)min,(double)xvals[0][maxloc],(double)fwhm/2.35,(double)max,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		int[] fixes={0,0,0,0,1,1,1,1,1,1,1,1,1};
		c2=0.0f;
		iterations=0;
		checkc2=false;

		double[] stats=new double[2];
		double[][] constraints=new double[2][13];
		constraints[0][0]=-0.5*(double)max; constraints[1][0]=(double)max;
		constraints[0][1]=(double)xvals[0][0]; constraints[1][1]=(double)xvals[0][length-1];
		constraints[0][2]=0.1*(double)(xvals[0][1]-xvals[0][0]); constraints[1][2]=(double)(xvals[0][length-1]-xvals[0][0]);
		constraints[0][3]=0.0; constraints[1][3]=10.0*(float)max;
		for(int i=1;i<4;i++){
			for(int j=0;j<3;j++){
				constraints[0][j+i*3+1]=constraints[0][j+1];
				constraints[1][j+i*3+1]=constraints[1][j+1];
			}
		}

		pw.addPoints(xvals[0],new float[length],false);
		series=pw.getNpts().length-1;

		tempx=new float[length];
		System.arraycopy(xvals[0],0,tempx,0,length);
		tempdata=new float[length];
		System.arraycopy(yvals[0],0,tempdata,0,length);
		fitclass=new NLLSfit(this,0.0001,50,0.1);

		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Manual_Fit",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean manual=gd.getNextBoolean();

		while(showoptions(params,fixes)){
			if(checkc2){
				fitclass.maxiter=0;
			} else {
				fitclass.maxiter=50;
			}
			float[] fit=fitclass.fitdata(params,fixes,constraints,tempdata,null,stats,true);
			pw.updateSeries(fit,series,false);
			c2=(float)stats[1];
			iterations=(int)stats[0];
			if(!manual){break;}
		}

		TextWindow outtable=jutils.selectTable("MultiGauss Fits");
		if(outtable==null){
			outtable=FitDialog.make_outtable("MultiGauss Fits",paramnames);
		}
		StringBuffer sb=new StringBuffer();
		sb.append(pw.getTitle());
		IJ.log("Chi Squared = "+(float)stats[1]); sb.append("\t"+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]); sb.append("\t"+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(paramnames[i]+" : "+params[i]+" : Fixed : "+(fixes[i]==1));  sb.append("\t"+(float)params[i]);
		}
		outtable.append(sb.toString());
	}

	public boolean showoptions(double[] params,int[] fixes){
		GenericDialog gd=new GenericDialog("Starting Fit Parameters");
		gd.addCheckbox("Check Chi Squared",checkc2);
		for(int i=0;i<params.length;i++){
			gd.addNumericField(paramnames[i],params[i],5,10,null);
			gd.addCheckbox("Fix"+(i+1)+"?",(fixes[i]==1));
		}
		gd.addNumericField("Iterations",iterations,0);
		gd.addNumericField("Chi Squared",c2,5,10,null);
		gd.addDialogListener(this);
		gd.showDialog();
		if(gd.wasCanceled()){return false;}
		checkc2=gd.getNextBoolean();
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){fixes[i]=1;}
			else{fixes[i]=0;}
		}
		return true;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
		checkc2=gd.getNextBoolean();
		double[] params=new double[paramnames.length];
		int[] fixes=new int[paramnames.length];
		for(int i=0;i<params.length;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){fixes[i]=1;}
			else{fixes[i]=0;}
		}
		NLLSfit fitclass=new NLLSfit(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,tempdata,null,stats,true);
		pw.updateSeries(fit,series,false);
		c2=(float)stats[1];
		return true;
	}

	public double fitfunc(double[] params,int indvar){
		//parameters are baseline,xc1,stdev1,amp1,xc2,stdev2,amp2,xc3,stdev3,amp3,xc4,stdev4,amp4
		double retval=params[0];
		for(int i=0;i<4;i++){
			if(params[3+i*3]>0.0){
				double xshift2=(params[1+i*3]-tempx[indvar])*(params[1+i*3]-tempx[indvar]);
				retval+=params[3+i*3]*Math.exp(-xshift2/(2.0*params[2+i*3]*params[2+i*3]));
			}
		}
		return retval;
	}

	public void showresults(String results){
		IJ.log(results);
	}
}
