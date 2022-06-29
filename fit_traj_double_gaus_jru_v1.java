/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import jalgs.*;
import jalgs.jfit.*;
import ij.text.*;
import jalgs.jsim.*;

public class fit_traj_double_gaus_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	public float[] tempx;
	public gausfunc gf;
	public boolean output;
	public rngs random;

	public void run(String arg) {
		output=true;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Stdev_Guess",50.0f,5,15,null);
		gd.addNumericField("Stdev_ratio_guess",1.0f,5,15,null);
		gd.addNumericField("Min_Stdev",25.0f, 5,15,null);
		gd.addNumericField("Max_Stdev",100.0f,5,15,null);
		gd.addCheckbox("Get_Errors",true);
		gd.addCheckbox("Single_Gaus_Fit",false);
		gd.addCheckbox("Output_Ratio_Parameters",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float stdevguess=(float)gd.getNextNumber();
		float srguess=(float)gd.getNextNumber();
		float minstdev=(float)gd.getNextNumber();
		float maxstdev=(float)gd.getNextNumber();
		boolean errs=gd.getNextBoolean();
		boolean single=gd.getNextBoolean();
		boolean outratio=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4SelCopy(iw);
		tempx=pw.getXValues()[0];
		float[] yvals=pw.getYValues()[0];
		gf=new gausfunc();
		random=new rngs();
		//first find the maximum
		float maxpos=jstatistics.getstatistic("MaxPos",yvals,null);
		float x1=tempx[(int)maxpos];
		float min=jstatistics.getstatistic("Min",yvals,null);
		float max=jstatistics.getstatistic("Max",yvals,null);
		//now find the highest value more than 1.5 guessstdev away from the first
		float x2=tempx[0]-1.0f;
		float maxaway=min;
		for(int i=0;i<tempx.length;i++){
			if((float)Math.abs(x1-tempx[i])>2.0f*stdevguess){
				if(yvals[i]>maxaway){maxaway=yvals[i]; x2=tempx[i];}
			}
		}
		//now that we have an x2 estimate, we can start searching
		double xcguess=0.5*(x2+x1);
		double dguess=Math.abs(x2-x1);
		IJ.log("xcguess = \t"+xcguess+"\t dguess = \t"+dguess+"\t amp1 = \t"+(max-min)+" \t amp2 = \t"+(maxaway-min));
		//search over xc, d, and stdev
		//should have relatively tight constraints on stdevratio
		//parameters are baseline, amp, ampratio,xc,d,stdev,stdevratio
		double[] params={min,max-min,1.0,xcguess,dguess,stdevguess,1.0};
		if(single) params=new double[]{min,max-min,0.0,xcguess,0.0,stdevguess,0.0};
		//IJ.log("guess params = "+table_tools.print_float_array(algutils.convert_arr_float(params),1));
		int[] fixes={0,0,0,1,1,1,0};
		if(single) fixes=new int[]{0,0,1,1,1,1,1};
		double dstart=dguess*0.25;
		double dend=dstart+2.0*(dguess-dstart);
		double xcstart=xcguess-dguess*0.25;
		double xcend=xcstart+2.0*(xcguess-xcstart);
		if(single){
			xcstart=xcguess-stdevguess*2.0;
			xcend=xcstart+stdevguess*4.0;
		}
		double stdevstart=minstdev;
		double stdevend=maxstdev;
		double[][] constraints=new double[2][7];
		constraints[0][0]=params[0]-params[1]; constraints[1][0]=params[0]+params[1];
		constraints[0][1]=0.0; constraints[1][1]=10.0*params[1];
		constraints[0][2]=0.2; constraints[1][2]=5.0;
		constraints[0][3]=xcstart; constraints[1][3]=xcend;
		constraints[0][4]=dstart; constraints[1][4]=dend;
		constraints[0][5]=stdevstart; constraints[1][5]=stdevend;
		constraints[0][6]=0.5; constraints[1][6]=2.0;
		if(single){
			dstart=dguess; dend=dguess;
		}
		double dd=0.05*dguess;
		double dxc=0.05*dguess;
		if(single) dxc=0.1*stdevguess;
		double dstdev=0.05*stdevguess;
		double mind=dstart; double minxc=xcstart; double minstdev2=stdevstart;
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,10,0.0);
		double[] stats=new double[2];
		fitclass.fitdata(params,fixes,constraints,yvals,null,stats,false);
		double c2min=stats[1];
		double[] oparams=params.clone();
		for(double d=dstart;d<=dend;d+=dd){
			for(double xc=xcstart;xc<=xcend;xc+=dxc){
				for(double stdev=stdevstart;stdev<=stdevend;stdev+=dstdev){
					params=oparams.clone();
					params[3]=xc; params[4]=d; params[5]=stdev;
					fitclass.fitdata(params,fixes,constraints,yvals,null,stats,false);
					if(stats[1]<c2min){
						mind=params[4]; minxc=params[3]; minstdev2=params[5]; c2min=stats[1];
					}
				}
			}
		}
		IJ.log("mind \t"+mind+"\t minxc \t"+minxc+"\t minstdev2 \t"+minstdev2);
		params=oparams.clone();
		params[3]=minxc; params[4]=mind; params[5]=minstdev2;
		if(single) params[4]=0.0;
		fixes=new int[]{0,0,0,0,0,0,0};
		if(single) fixes=new int[]{0,0,1,0,1,0,1};
		float[] fit=fitclass.fitdata(params,fixes,constraints,yvals,null,stats,false);
		pw.addPoints(tempx,fit,true);
		pw.updatePlot();
		String[] paramnames={"baseline","amp","ampratio","xc","d","stdev","stdevratio"};
		if(!outratio){
			paramnames[1]="amp1"; paramnames[2]="xc1"; paramnames[3]="stdev1";
			paramnames[4]="amp2"; paramnames[5]="xc2"; paramnames[6]="stdev2";
		}
		TextWindow outtable=jutils.selectTable("Double Gauss Fits");
		if(outtable==null){
			outtable=FitDialog.make_outtable("Double Gauss Fits",paramnames);
		}
		float[] err=null;
		if(errs){
			output=false;
			monte_carlo_errors_v2 mcerrs=new monte_carlo_errors_v2(this,0.0001,10,false,0.1);
			double[][] temp1=mcerrs.geterrors(params,fixes,constraints,yvals,null,1000);
			float[][] temp=algutils.convert_arr_float(temp1);
			err=new float[temp.length];
			for(int i=0;i<err.length;i++) err[i]=jstatistics.getstatistic("StDev",temp[i],null);
		}
		//now convert the parameters away from ratio parameters if desired
		double[][] convparams={params,algutils.convert_arr_double(err)};
		if(!outratio){
			convparams=convertParams(params,algutils.convert_arr_double(err));
			fixes=new int[]{0,0,0,0,0,0,0};
		}
		StringBuffer sb=new StringBuffer();
		sb.append(pw.getTitle());
		IJ.log("Chi Squared = "+(float)stats[1]); sb.append("\t"+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]); sb.append("\t"+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(paramnames[i]+" = "+convparams[0][i]);  sb.append("\t"+(float)convparams[0][i]);
		}
		outtable.append(sb.toString());
		if(errs){
			StringBuffer sb2=new StringBuffer();
			sb2.append(pw.getTitle()+"_errs");
			sb2.append("\t"+(float)err[err.length-1]);
			sb2.append("\t"+0.0f);
			int counter=0;
			for(int i=0;i<params.length;i++){
				if(fixes[i]==0){
					sb2.append("\t"+(float)convparams[1][counter]);
					counter++;
				} else {
					sb2.append("\t0.0 ");
				}
			}
			outtable.append(sb2.toString());
		}
	}

	public double[] fitfunc(double[] params){
		//params are 0base, 1ampavg, 2ampratio, 3xc, 4d, 5sdavg, 6sdratio
		double a2=2.0*params[1]/(1.0+params[2]);
		double a1=2.0*params[1]-a2;
		double s2=2.0*params[5]/(1.0+params[6]);
		double s1=2.0*params[5]-s2;
		double xc1=params[3]-0.5*params[4];
		double xc2=params[3]+0.5*params[4];
		double[] fit=new double[tempx.length];
		for(int i=0;i<fit.length;i++){
			double xd1=xc1-(double)tempx[i];
			double xd2=xc2-(double)tempx[i];
			fit[i]=params[0]+a1*gf.getinterpgaus(xd1,s1)+a2*gf.getinterpgaus(xd2,s2);
		}
		return fit;
	}

	public void showresults(String results){
		if(output) IJ.log(results);
	}

	public double[][] convertParams(double[] params,double[] sems){
		//this converts peak pair parameters (using ratios) to individual parameters and their sem values
		//output parameters are 0base,1a1,2xc1,3sd1,4a2,5xc2,6sd2
		if(params[2]==0.0){
			double[] newparams={params[0],params[1],params[3],params[5],0.0,0.0,0.0};
			double[] newsems=new double[newparams.length];
			if(sems!=null){
				newsems[0]=sems[0]; newsems[1]=sems[1]; newsems[2]=sems[2]; newsems[3]=sems[3];
			}
			return new double[][]{newparams,newsems};
		}
		double[] sems2=new double[params.length];
		if(sems!=null){
			sems2=sems.clone();
		}
		double[] ampconverted=avgRatio2Vals(params[1],params[2],sems[1],sems[2]);
		double[] sdconverted=avgRatio2Vals(params[5],params[6],sems[5],sems[6]);
		double xc1=params[3]-0.5*params[4];
		double xc2=xc1+params[4];
		double xc1sem=0.0;
		double xc2sem=0.0;
		if(sems!=null){
			xc1sem=Math.sqrt(sems[3]*sems[3]+0.25*sems[4]*sems[4]);
			xc2sem=Math.sqrt(xc1sem*xc1sem+sems[4]*sems[4]);
		}
		double[] newparams={params[0],ampconverted[0],xc1,sdconverted[0],ampconverted[1],xc2,sdconverted[1]};
		double[] newsems={sems[0],ampconverted[2],xc1sem,sdconverted[2],ampconverted[3],xc2sem,sdconverted[3]};
		return new double[][]{newparams,newsems};
	}

	public double[] avgRatio2Vals(double avg,double ratio,double avgsem,double ratiosem){
		//converts average and ratio parameters to individual paramaters
		//sem values are propagated by monte carlo simulations
		double val2=2.0*avg/(1.0+ratio);
		double val1=2.0*avg-val2;
		double val1sem=0.0;
		double val2sem=0.0;
		if(avgsem>0.0 && ratiosem>0.0){
			float[][] sims=new float[2][10000];
			for(int i=0;i<10000;i++){
				double avgval=random.gasdev(avg,avgsem);
				double ratioval=random.gasdev(ratio,ratiosem);
				sims[1][i]=(float)(2.0*avgval/(1.0+ratioval));
				sims[0][i]=(float)(2.0*avgval-sims[1][i]);
			}
			val1sem=(double)jstatistics.getstatistic("StDev",sims[0],null);
			val2sem=(double)jstatistics.getstatistic("StDev",sims[1],null);
		}
		return new double[]{val1,val2,val1sem,val2sem};
	}
}
