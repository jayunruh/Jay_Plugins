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
import ij.plugin.frame.RoiManager;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import ij.text.*;
import java.util.*;

public class fit_multi_gaussian_jru_v3 implements PlugIn, NLLSfitinterface_v2, DialogListener {
	boolean checkc2;
	float[] tempx,tempdata;
	float c2;
	int iterations,series;
	NLLSfit_v2 fitclass;
	String[] paramnames={"baseline","xc1","stdev1","amp1","xc2","stdev2","amp2","xc3","stdev3","amp3","xc4","stdev4","amp4"};
	PlotWindow4 pw;

	public void run(String arg) {
		//this plugin takes multiple spot selections and uses them as a starting point for multigaussian fitting
		ImageWindow iw=WindowManager.getCurrentWindow();
		pw=jutils.getPW4SelCopy(iw);
		String title=pw.getTitle();
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		int[] colors=pw.getColors();
		colors[0]=0;
		int length=yvals[0].length;
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.showMessage("No Peaks are Selected");
			return;
		}
		//get all of the peak positions on the x axis
		Roi[] rois=rman.getRoisAsArray();
		int npeaks=rois.length;
		if(npeaks>4) npeaks=4;
		float[] coords=new float[npeaks];
		Plot4 plot=pw.getPlot();
		for(int i=0;i<npeaks;i++){
			Rectangle r=rois[i].getBounds();
			coords[i]=plot.getPlotCoords(r.x,r.y)[0];
		}

		tempx=new float[length];
		System.arraycopy(xvals[0],0,tempx,0,length);
		tempdata=new float[length];
		System.arraycopy(yvals[0],0,tempdata,0,length);

		//now grid search over the standard deviations assuming they are all the same
		float xinc=tempx[1]-tempx[0];
		float minstdev=xinc;
		float maxstdev=(0.5f*xinc*(float)tempx.length)/(float)npeaks;
		double[][] beststdev=get_best_stdev(minstdev,maxstdev,0.25f*xinc,coords);
		//now initialize the parameters
		//parameters are baseline,xc1,stdev1,amp1,xc2,stdev2,amp2,xc3,stdev3,amp3,xc4,stdev4,amp4
		double[] params=new double[13];
		int[] fixes=new int[13]; Arrays.fill(fixes,1);
		float max=jstatistics.getstatistic("Max",tempdata,null);
		double[] stats=new double[2];
		double[][] constraints=new double[2][13];
		params[0]=beststdev[0][0];
		fixes[0]=0;
		constraints[0][0]=-0.5*(double)max; constraints[1][0]=(double)max;
		for(int i=0;i<npeaks;i++){
			params[3*i+1]=coords[i];
			params[3*i+2]=beststdev[1][0];
			params[3*i+3]=beststdev[0][i+1];
			fixes[3*i+1]=0;
			fixes[3*i+2]=0;
			fixes[3*i+3]=0;
			constraints[0][3*i+1]=coords[i]-beststdev[1][0];
			constraints[1][3*i+1]=coords[i]+beststdev[1][0];
			constraints[0][3*i+2]=0.1f*beststdev[1][0];
			constraints[1][3*i+2]=10.0f*beststdev[1][0];
			constraints[0][3*i+3]=0.1f*params[3*i+3];
			constraints[1][3*i+3]=10.0f*params[3*i+3];
		}
		c2=0.0f;
		iterations=0;
		checkc2=false;

		pw.addPoints(xvals[0],new float[length],false);
		series=pw.getNpts().length-1;

		fitclass=new NLLSfit_v2(this,0.0001,50,0.1);

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
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,null,tempdata,null,stats,true);
		pw.updateSeries(fit,series,false);
		c2=(float)stats[1];
		return true;
	}

	public double[] fitfunc(double[] params){
		//parameters are baseline,xc1,stdev1,amp1,xc2,stdev2,amp2,xc3,stdev3,amp3,xc4,stdev4,amp4
		double[] retarray=new double[tempx.length];
		for(int k=0;k<tempx.length;k++){
			retarray[k]=params[0];
			for(int i=0;i<4;i++){
				if(params[3+i*3]>0.0){
					double xshift2=(params[1+i*3]-tempx[k])*(params[1+i*3]-tempx[k]);
					retarray[k]+=params[3+i*3]*Math.exp(-xshift2/(2.0*params[2+i*3]*params[2+i*3]));
				}
			}
		}
		return retarray;
	}

	public double[][] get_leastsquares_amps(float[] centers,double stdev){
		//least squares gives baseline and amps
		//return params are base,amp1,amp2,...
		double[][] indvars=new double[centers.length+1][];
		indvars[0]=new double[tempx.length];
		Arrays.fill(indvars[0],1.0);
		for(int i=0;i<centers.length;i++){
			double[] tempparams=new double[13];
			tempparams[1]=(double)centers[i]; tempparams[2]=stdev; tempparams[3]=1.0f;
			indvars[i+1]=fitfunc(tempparams);
		}
		linleastsquares lls=new linleastsquares(indvars);
		double[][] output=new double[2][];
		output[0]=lls.fitdata(tempdata,null);
		output[1]=new double[1];
		output[1][0]=lls.get_c2(output[0],tempdata,null);
		return output;
	}

	public double[][] get_best_stdev(float minstdev,float maxstdev,float stdevinc,float[] centers){
		double[][] temp=get_leastsquares_amps(centers,minstdev);
		double minc2=temp[1][0];
		double best=(double)minstdev;
		for(float i=minstdev+stdevinc;i<maxstdev;i+=stdevinc){
			temp=get_leastsquares_amps(centers,i);
			double c2=temp[1][0];
			if(c2<minc2){
				minc2=c2;
				best=i;
			}
		}
		temp=get_leastsquares_amps(centers,best);
		temp[1]=new double[]{best,minc2};
		return temp;
	}

	public void showresults(String results){
		IJ.log(results);
	}
}
