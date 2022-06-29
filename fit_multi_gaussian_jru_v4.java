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

public class fit_multi_gaussian_jru_v4 implements PlugIn, NLLSfitinterface_v2 {
	boolean checkc2;
	float[] tempx,tempdata;
	float c2;
	int iterations,series,npeaks;
	NLLSfit_v2 fitclass;
	String[] paramnames;
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
		if(rman==null || rman.getCount()<1){ //if there are no rois, assume roi is at maximum and we are fitting a single gaussian
			IJ.log("No Peaks are Selected, Assuming Single Peak");
			//return;
			if(rman==null) rman=new RoiManager();
			int maxpos=(int)jstatistics.getstatistic("MaxPos",yvals[0],null);
			int[] maxcoords=pw.getPlot().getCoordsPosition(xvals[0][maxpos],yvals[0][maxpos]);
			rman.addRoi(new PointRoi(maxcoords[0],135));
		}
		//get all of the peak positions on the x axis
		Roi[] rois=rman.getRoisAsArray();
		npeaks=rois.length;
		if(npeaks>8) npeaks=8;
		float[] coords=new float[npeaks];
		Plot4 plot=pw.getPlot();
		String[] paramnames=new String[3*npeaks+1];
		paramnames[0]="baseline";
		for(int i=0;i<npeaks;i++){
			paramnames[3*i+1]="xc"+(i+1);
			paramnames[3*i+2]="stdev"+(i+1);
			paramnames[3*i+3]="amp"+(i+1);
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
		double[] params=new double[1+3*npeaks];
		int[] fixes=new int[1+3*npeaks]; Arrays.fill(fixes,1);
		float max=jstatistics.getstatistic("Max",tempdata,null);
		double[] stats=new double[2];
		double[][] constraints=new double[2][1+3*npeaks];
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

		FitDialog_v3 fd=new FitDialog_v3(pw,this,paramnames);
		TextWindow outtable=jutils.selectTable("MultiGauss Fits");
		if(outtable==null){
			outtable=fd.make_outtable("MultiGauss Fits",paramnames);
		}
		//WindowManager.putBehind();
		//WindowManager.setCurrentWindow(pw);
		//WindowManager.toFront(pw);
		IJ.selectWindow(title);
		fd.run_fit(params,null,constraints,fixes);
		fd.append_outtable_params("MultiGauss Fits",pw.getTitle(),params);
		fd.append_outtable_errs("MultiGauss Fits",pw.getTitle());
	}

	public double[] fitfunc(double[] params){
		//parameters are baseline,xc1,stdev1,amp1,xc2,stdev2,amp2,xc3,stdev3,amp3,xc4,stdev4,amp4
		double[] retarray=new double[tempx.length];
		for(int k=0;k<tempx.length;k++){
			retarray[k]=params[0];
			for(int i=0;i<npeaks;i++){
				if(params[3+i*3]!=0.0){
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
			double[] tempparams=new double[1+npeaks*3];
			tempparams[1]=(double)centers[i]; tempparams[2]=stdev; tempparams[3]=1.0f;
			indvars[i+1]=fitfunc(tempparams);
		}
		linleastsquares lls=new linleastsquares(indvars);
		double[][] output=new double[2][];
		output[0]=lls.fitdata(tempdata,null);
		output[1]=new double[1];
		output[1][0]=lls.get_c2(output[0],tempdata,null);
		//new PlotWindow4("test indvars","x","y",algutils.convert_arr_float(indvars),null).draw();
		//IJ.log(""+output[0][0]+" , "+indvars[1][0]);
		return output;
	}

	public double[][] get_best_stdev(float minstdev,float maxstdev,float stdevinc,float[] centers){
		double minc2=Float.MAX_VALUE;
		double best=(double)minstdev;
		for(float i=minstdev;i<maxstdev;i+=stdevinc){
			double[][] temp=get_leastsquares_amps(centers,i);
			double c2=temp[1][0];
			//if(i==minstdev) break;
			if(c2!=Float.NaN && c2<minc2){
				minc2=c2;
				best=i;
			}
			//IJ.log(""+i+" , "+c2);
		}
		double[][] temp=get_leastsquares_amps(centers,best);
		temp[1]=new double[]{best,minc2};
		return temp;
	}

	public void showresults(String results){
		IJ.log(results);
	}
}
