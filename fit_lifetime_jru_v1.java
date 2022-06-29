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

public class fit_lifetime_jru_v1 implements PlugIn, NLLSfitinterface {
	boolean wraparound,checkc2;
	float c2;
	double[] currfit,currparams,irf,irfshift;
	double currshift,nschan;
	int iterations,startfit,endfit,length;
	float[] tempdata;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		wraparound=false;
		gd.addCheckbox("Wrap Around Data",wraparound);
		String[] fitalgs={"NLLS","Simplex"};
		gd.addChoice("Fit Algorithm",fitalgs,fitalgs[0]);
		gd.addCheckbox("Manual_Fit",true);
		gd.addCheckbox("Table_Output",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		wraparound=gd.getNextBoolean();
		int algindex=gd.getNextChoiceIndex();
		boolean manual=gd.getNextBoolean();
		boolean tableout=gd.getNextBoolean();

		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4Copy(iw);
		String title=pw.getTitle();
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		nschan=(double)(xvals[0][1]-xvals[0][0]);
		length=yvals[0].length;
		startfit=0; endfit=length-1;
		irf=new double[length];
		float[] data;
		if(yvals.length==1){
			//need to read in the irf separately
			int[] wList = WindowManager.getIDList();
			String[] titles = new String[wList.length];
			for(int i=0;i<wList.length;i++){
				ImagePlus imp = WindowManager.getImage(wList[i]);
				ImageWindow iw1=imp.getWindow();
				if(iw1.getClass().getName().equals("jguis.PlotWindow4")){titles[i]=imp.getTitle();}
				else{titles[i]="NA";}
			}
			GenericDialog gd2 = new GenericDialog("Options");
			gd2.addChoice("Scatter Function",titles,titles[0]);
			gd2.showDialog();
			if(gd2.wasCanceled()){return;}
			int index1 = gd2.getNextChoiceIndex();
			ImagePlus imp = WindowManager.getImage(wList[index1]);
			ImageWindow iw1=imp.getWindow();
			float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw1,"getYValues");
			for(int i=0;i<length;i++){irf[i]=(double)yvals2[0][i];}
			data=yvals[0];
		} else {
			data=yvals[1];
			for(int i=0;i<length;i++){irf[i]=(double)yvals[0][i];}
		}
		
		//parameters are shift,baseline,a1,t1,...
		double[] params={0.0,0.0,1.0,3.0,0.0,1.0,0.0,0.2,0.0,0.2};
		int[] fixes={1,1,0,0,1,1,1,1,1,1};
		c2=0.0f;
		iterations=0;
		checkc2=false;

		double[] stats=new double[2];
		double[][] constraints=new double[2][10];
		constraints[0][0]=-10.0; constraints[1][0]=10.0;
		constraints[0][1]=-1000.0; constraints[1][1]=10000.0;
		constraints[0][2]=-100.0; constraints[1][2]=1000.0;
		constraints[0][3]=0.01; constraints[1][3]=1000.0;
		constraints[0][4]=-100.0; constraints[1][4]=1000.0;
		constraints[0][5]=0.01; constraints[1][5]=1000.0;
		constraints[0][6]=-100.0; constraints[1][6]=1000.0;
		constraints[0][7]=0.01; constraints[1][7]=1000.0;
		constraints[0][8]=-100.0; constraints[1][8]=1000.0;
		constraints[0][9]=0.01; constraints[1][9]=1000.0;

		pw.addPoints(xvals[0],new float[length],false);
		int series=pw.getNpts().length-1;

		while(showoptions(params,fixes)){
			tempdata=new float[endfit-startfit+1];
			float[] weights=new float[endfit-startfit+1];
			for(int i=startfit;i<=endfit;i++){
				tempdata[i-startfit]=data[i];
				if(tempdata[i-startfit]>0.0f){
					weights[i-startfit]=1.0f/tempdata[i-startfit];
				} else {
					weights[i-startfit]=1.0f;
				}
			}
			currparams=new double[10];
			for(int i=0;i<10;i++){
				currparams[i]=params[i];
			}
			currshift=params[0]+1.0;
			updatefit(params);
			if(algindex==0){
				NLLSfit fitclass;
				if(checkc2){
					fitclass=new NLLSfit(this,0);
				} else {
					fitclass=new NLLSfit(this,0.0001,50,0.1);
				}
				float[] fit=fitclass.fitdata(params,fixes,constraints,tempdata,weights,stats,true);
				float[] fullfit=new float[length];
				for(int i=0;i<startfit;i++){fullfit[i]=0.0f;}
				for(int i=startfit;i<=endfit;i++){fullfit[i]=fit[i-startfit];}
				for(int i=endfit+1;i<length;i++){fullfit[i]=0.0f;}
				pw.updateSeries(fullfit,series,false);
				c2=(float)stats[1];
				iterations=(int)stats[0];
			} else {
				simplexfit fitclass;
				if(checkc2){
					fitclass=new simplexfit(this,0,0);
				} else {
					fitclass=new simplexfit(this,100,20);
				}
				float[] fit=fitclass.fitdata(params,fixes,constraints,tempdata,weights,stats,true);
				float[] fullfit=new float[length];
				for(int i=0;i<startfit;i++){fullfit[i]=0.0f;}
				for(int i=startfit;i<=endfit;i++){fullfit[i]=fit[i-startfit];}
				for(int i=endfit+1;i<length;i++){fullfit[i]=0.0f;}
				pw.updateSeries(fullfit,series,false);
				c2=(float)stats[1];
				iterations=(int)stats[0];
			}
			if(!manual){break;}
		}

		if(tableout){
			TextWindow tw=jutils.selectTable("Lifetime_Output");
			if(tw==null) tw=new TextWindow("Lifetime_Output","c2\titer\tshift\tbaseline\tamp1\ttau1\tamp2\ttau2\tamp3\ttau3\tamp4\ttau4","",400,200);
			tw.append(""+(float)stats[1]+"\t"+(int)stats[0]+"\t"+(float)params[0]+"\t"+(float)params[1]+"\t"+(float)params[2]+"\t"+(float)params[3]+"\t"+(float)params[4]+"\t"+(float)params[5]+"\t"+(float)params[6]+"\t"+(float)params[7]+"\t"+(float)params[8]+"\t"+(float)params[9]+"\n");
		} else {
			IJ.log("Chi Squared = "+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		IJ.log("shift = "+(float)params[0]+" fixed = "+fixes[0]);
		IJ.log("baseline = "+(float)params[1]+" fixed = "+fixes[1]);
		IJ.log("amp1 = "+(float)params[2]+" fixed = "+fixes[2]);
		IJ.log("tau1 = "+(float)params[3]+" fixed = "+fixes[3]);
		IJ.log("amp2 = "+(float)params[4]+" fixed = "+fixes[4]);
		IJ.log("tau2 = "+(float)params[5]+" fixed = "+fixes[5]);
		IJ.log("amp3 = "+(float)params[6]+" fixed = "+fixes[6]);
		IJ.log("tau3 = "+(float)params[7]+" fixed = "+fixes[7]);
		IJ.log("amp4 = "+(float)params[8]+" fixed = "+fixes[8]);
		IJ.log("tau4 = "+(float)params[9]+" fixed = "+fixes[9]);
		}
	}

	public boolean showoptions(double[] params,int[] fixes){
		GenericDialog gd=new GenericDialog("Starting Fit Parameters");
		gd.addCheckbox("Check Chi Squared",checkc2);
		gd.addNumericField("Start Fit",(float)startfit,0,10,null);
		gd.addNumericField("End Fit",(float)endfit,0,10,null);
		gd.addNumericField("Shift",(float)params[0],5,10,null);
		gd.addCheckbox("Fix1?",(fixes[0]==1));
		gd.addNumericField("Baseline",(float)params[1],5,10,null);
		gd.addCheckbox("Fix2?",(fixes[1]==1));
		gd.addNumericField("amp_1",(float)params[2],5,10,null);
		gd.addCheckbox("Fix3?",(fixes[2]==1));
		gd.addNumericField("tau_1",(float)params[3],5,10,null);
		gd.addCheckbox("Fix4?",(fixes[3]==1));
		gd.addNumericField("amp_2",(float)params[4],5,10,null);
		gd.addCheckbox("Fix5?",(fixes[4]==1));
		gd.addNumericField("tau_2",(float)params[5],5,10,null);
		gd.addCheckbox("Fix6?",(fixes[5]==1));
		gd.addNumericField("amp_3",(float)params[6],5,10,null);
		gd.addCheckbox("Fix7?",(fixes[6]==1));
		gd.addNumericField("tau_3",(float)params[7],5,10,null);
		gd.addCheckbox("Fix8?",(fixes[7]==1));
		gd.addNumericField("amp_4",(float)params[8],5,10,null);
		gd.addCheckbox("Fix9?",(fixes[8]==1));
		gd.addNumericField("tau_4",(float)params[9],5,10,null);
		gd.addCheckbox("Fix10?",(fixes[9]==1));
		gd.addNumericField("Iterations",iterations,0);
		gd.addNumericField("Chi Squared",c2,7,10,null);
		gd.showDialog();
		if(gd.wasCanceled()){return false;}
		checkc2=gd.getNextBoolean();
		startfit=(int)gd.getNextNumber();
		endfit=(int)gd.getNextNumber();
		for(int i=0;i<10;i++){
			params[i]=gd.getNextNumber();
			if(gd.getNextBoolean()){fixes[i]=1;}
			else{fixes[i]=0;}
		}
		return true;
	}

	public double fitfunc(double[] params,int indvar){
		//params are shift,baseline,a1,t1,a2,t2,a3,t3,a4,t4;
		for(int i=0;i<10;i++){
			if(params[i]!=currparams[i]){
				updatefit(params);
				for(int j=i;j<10;j++){
					currparams[j]=params[j];
				}
				break;
			}
		}
		return params[1]+nschan*currfit[indvar+startfit];
	}

	private void updatefit(double[] params){
		//start by updating the irfshift array if necessary
		if(params[0]!=currshift){
			updateirf(params[0]);
			currshift=params[0];
		}
		if(wraparound){
			double[] model=new double[2*length];
			for(int i=0;i<2*length;i++){
				double time=nschan*(double)i;
				for(int j=0;j<4;j++){
					model[i]+=params[j*2+2]*Math.exp(-time/params[j*2+3]);
				}
			}
			currfit=convolute(model,irfshift);
		} else {
			double[] model=new double[length];
			for(int i=0;i<length;i++){
				double time=nschan*(double)i;
				for(int j=0;j<4;j++){
					model[i]+=params[j*2+2]*Math.exp(-time/params[j*2+3]);
				}
			}
			currfit=convolute(model,irfshift);
		}
	}

	private void updateirf(double shift){
		irfshift=new double[length];
		if(shift>0.0){
			int intshift=(int)shift;
			double shiftrem=shift-(double)intshift;
			for(int i=1;i<(length-intshift);i++){
				irfshift[i+intshift]=irf[i]-(irf[i]-irf[i-1])*shiftrem;
			}
		} else {
			if(shift<0.0){
				int intshift=(int)Math.abs(shift);
				double shiftrem=-(double)intshift-shift;
				for(int i=intshift;i<(length-1);i++){
					irfshift[i-intshift]=irf[i]+(irf[i+1]-irf[i])*shiftrem;
				}
			} else {
				for(int i=1;i<(length-1);i++){
					irfshift[i]=irf[i];
				}
			}
		}
	}

	public void showresults(String results){
		IJ.log(results);
	}

	double[] convolute(double[] model,double[] response){
		int rlength=response.length;
		int mlength=model.length;
		if(mlength>2*rlength){mlength=2*rlength;}
		double[] conv=new double[rlength];
		for(int i=0;i<rlength;i++){
			for(int j=i;j<mlength;j++){
				if(j<rlength){
					conv[j] += response[i]*model[j-i];
				} else {
					if(wraparound){
						conv[j-rlength]+=response[i]*model[j-i];
					}
				}
			}
		}
		return conv;
	}
}
