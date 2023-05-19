/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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

public class fit_stretched_exp_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	//fits to x<xshift?b:b+A*(1-exp(ln(0.5)*(((x-xshift)/(EC50-xshift))^alpha))))
	int model;
	float[] xvals,yvals;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_EC50_",1.00,5,15,null);
		gd.addNumericField("Max_EC50_",500.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		double minec50=gd.getNextNumber();
		double maxec50=gd.getNextNumber();

		String[][] paramsnames={{"baseline","amp","EC50","alpha","xshift"}};
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4SelCopy(iw);
		FitDialog_v3 fd=new FitDialog_v3(pw,this,paramsnames[0]);
		xvals=pw.getXValues()[0];
		yvals=pw.getYValues()[0];
		float[] errs=pw.getErrors(0,false);
		float[] weights=null;
		if(errs!=null){
			weights=new float[errs.length];
			for(int i=0;i<errs.length;i++) weights[i]=1.0f/(errs[i]*errs[i]);
		}
		double[] params=null;
		double[] guesses=initEC50(yvals,minec50,maxec50);
		params=new double[]{0.0,1.0,guesses[0],guesses[1],0.0};
		double[][] constraints=new double[2][5];
		constraints[0][0]=params[0]-params[1]; constraints[1][0]=params[0]+params[1];
		constraints[0][1]=0.1*params[1]; constraints[1][1]=10.0*params[1];
		constraints[0][2]=minec50; constraints[1][2]=maxec50;
		constraints[0][3]=0.1; constraints[1][3]=10.0;
		constraints[0][4]=0.0; constraints[1][4]=maxec50;
		//double[][] constraints=null;
		int[] fixes={1,1,0,0,0};
		fd.run_fit(params,weights,constraints,fixes);
		TextWindow outtable=jutils.selectTable("Stretched Exp Fits");
		if(outtable==null){
			outtable=fd.make_outtable("Stretched Exp Fits");
		}
		fd.append_outtable_params("Stretched Exp Fits",pw.getTitle(),params);
		fd.append_outtable_errs("Stretched Exp Fits",pw.getTitle());
	}

	public double getc2(double[] fit){
		double c2=0.0;
		for(int i=0;i<fit.length;i++){
			c2+=(fit[i]-(double)yvals[i])*(fit[i]-(double)yvals[i]);
		}
		c2/=(double)(fit.length-3);
		return c2;
	}

	public double[] initEC50(float[] data,double min,double max){
		double minc2=-1.0;
		double minec50=min;
		double minamp=0.0;
		double minoff=0.0;
		double minalpha=2.0;
		//added to check alpha values 2, 4, and 6
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,2.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=2.0;}
		}
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,4.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=4.0;}
		}
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,6.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=6.0;}
		}
		return new double[]{minec50,minalpha,minc2,minamp,minoff};
	}

	public double[] fitfunc(double[] params){
		//params are 0baseline,1amp,2ec50,3alpha,4xshift
		double[] func=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			double x=(double)xvals[i];
			if(x<params[4]) func[i]=params[0];
			else {
				double temp=(x-params[4])/(params[2]-params[4]);
				func[i]=params[0]+params[1]*(1.0-Math.exp(-0.693*Math.pow(temp,params[3])));
			}
		}
		return func;
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
