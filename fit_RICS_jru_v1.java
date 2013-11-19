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
import ij.plugin.filter.*;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;

public class fit_RICS_jru_v1 implements PlugInFilter, NLLSfitinterface {
	ImagePlus imp;
	double[] params;
	int[] fixes;
	double[] stats;
	boolean c2test;
	double pi = 3.14159265;
	int xc,yc,xpts,xskip,ypts,maxiter,photons,psfindex;
	float[] g0val;
	double pixelsize,pixeltime,linetime;
	String[] psfoptions={"3D Gaussian","2D Gaussian"};

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL;
	}

	public void run(ImageProcessor ip) {
		float[] pixels = (float[])ip.getPixels();
		int fullwidth = imp.getWidth();
		int fullheight = imp.getHeight();
		
		//initialize the variables
		xc=fullwidth/2; yc=fullheight/2; xpts=16; ypts=16; xskip=1; maxiter=10; c2test=false;
		photons=1; pixelsize=imp.getCalibration().pixelWidth; pixeltime=0.00001; linetime=0.00256;
		GenericDialog gd2=new GenericDialog("Options");
		psfindex=0;
		gd2.addChoice("PSF?",psfoptions,psfoptions[psfindex]);
		gd2.addNumericField("X center",xc,0);
		gd2.addNumericField("Y center",yc,0);
		gd2.addNumericField("Pts in x to fit",xpts,0);
		gd2.addNumericField("Pts in x to skip",xskip,0);
		gd2.addNumericField("Pts in y to fit",ypts,0);
		gd2.addNumericField("Photons",photons,0);
		gd2.addNumericField("Pixel Size (um)",(float)(pixelsize),5,10,null);
		gd2.addNumericField("Pixel Time (us)",(float)(pixeltime*1000000.0),5,10,null);
		gd2.addNumericField("Line Time (ms)",(float)(linetime*1000.0),5,10,null);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		psfindex=gd2.getNextChoiceIndex();
		xc=(int)gd2.getNextNumber();
		yc=(int)gd2.getNextNumber();
		xpts=(int)gd2.getNextNumber();
		xskip=(int)gd2.getNextNumber();
		ypts=(int)gd2.getNextNumber();
		photons=(int)gd2.getNextNumber();
		pixelsize=(gd2.getNextNumber());
		pixeltime=(gd2.getNextNumber())/1000000.0;
		linetime=(gd2.getNextNumber())/1000.0;
		params=new double[7];
		params[0]=0.3; params[1]=5.0; params[3]=1.0;
		params[4]=40.0; params[6]=5.0;
		fixes=new int[7];
		fixes[0]=fixes[1]=1;
		fixes[5]=fixes[6]=1;
		stats=new double[2];
		double[][] constraints={{0.1,0.0,-10.0,0.0,0.001,0.0,0.001},{2.0,20.0,10.0,1000.0,1000.0,1000.0,1000.0}};
		float[] ac=new float[xpts*ypts];
		float[][] acdisp=new float[xpts][ypts];
		float[][] fitdisp=new float[xpts][ypts];
		float[] xaxis=new float[xpts];
		float[] yaxis=new float[ypts];
		for(int i=0;i<ypts;i++){
			yaxis[i]=(float)pixelsize*(float)i;
			for(int j=0;j<xpts;j++){
				if(i==0){xaxis[j]=(float)pixelsize*(float)j;}
				ac[j+i*xpts]=pixels[(i+yc)*fullwidth+j+xc];
				acdisp[j][i]=ac[j+i*xpts];
			}
		}
		PlotWindow3D pw=new PlotWindow3D("RICS fit","x (\u00B5m)","y (\u00B5m)","G(\u03BE,\u03C8)",xaxis,yaxis,acdisp);
		pw.draw();
		pw.addPoints(xaxis,yaxis,fitdisp,true);
		g0val=new float[xskip];
		for(int i=0;i<xskip;i++){
			g0val[i]=ac[i];
		}

		while (showDialog())
		{
			NLLSfit nf;
			if(c2test){
				nf=new NLLSfit(this,0.0001,0,0.1);				
			}
			else{
				nf=new NLLSfit(this,0.0001,maxiter,0.1);
			}
			float[] fit=nf.fitdata(params,fixes,constraints,ac,null,stats,false);
			for(int i=0;i<ypts;i++){
				for(int j=0;j<xpts;j++){
					fitdisp[j][i]=fit[j+i*xpts];
				}
			}
			pw.updateSeries(xaxis,yaxis,fitdisp,1,true);
		}
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog("Fitting Parameters");
		gd.addCheckbox("Test Chi Squared?",c2test);
		gd.addNumericField("Max Iterations?",maxiter,0);
		gd.addNumericField("w0 (um)",params[0],10,15,"");
		gd.addCheckbox("fix?",(fixes[0]==1));
		gd.addNumericField("z0/w0",params[1],10,15,"");
		gd.addCheckbox("fix?",(fixes[1]==1));
		gd.addNumericField("Baseline",params[2],10,15,"");
		gd.addCheckbox("fix?",(fixes[2]==1));
		gd.addNumericField("G(0) 1",params[3],10,15,"");
		gd.addCheckbox("fix?",(fixes[3]==1));
		gd.addNumericField("D1 (um^2/sec)",params[4],10,15,"");
		gd.addCheckbox("fix?",(fixes[4]==1));
		gd.addNumericField("G(0) 2",params[5],10,15,"");
		gd.addCheckbox("fix?",(fixes[5]==1));
		gd.addNumericField("D2 (um^2/sec)",params[6],10,15,"");
		gd.addCheckbox("fix?",(fixes[6]==1));
		gd.addNumericField("Iterations Completed",(int)stats[0],0);
		gd.addNumericField("Chi Squared",(float)stats[1],5,15,"");
		gd.showDialog();
		if(gd.wasCanceled()){
			String temp;
			IJ.log("Fit Results");
			IJ.log("Points in x to fit "+xpts);
			IJ.log("Points in x to skip "+xskip);
			IJ.log("Points in y to fit "+ypts);
			temp = (fixes[0]==1) ? "true" : "false";
			IJ.log("Fix w0? "+temp);
			IJ.log("w0 (um) = "+(float)params[0]);
			temp = (fixes[1]==1) ? "true" : "false";
			IJ.log("Fix z0/w0? "+temp);
			IJ.log("z0/w0 = "+(float)params[1]);
			temp = (fixes[2]==1) ? "true" : "false";
			IJ.log("Fix Baseline? "+temp);
			IJ.log("Baseline = "+(float)params[2]);
			temp = (fixes[3]==1) ? "true" : "false";
			IJ.log("Fix G(0) 1? "+temp);
			IJ.log("G(0) 1 = "+(float)params[3]);
			temp = (fixes[4]==1) ? "true" : "false";
			IJ.log("Fix D1? "+temp);
			IJ.log("D1 = "+(float)params[4]);
			temp = (fixes[5]==1) ? "true" : "false";
			IJ.log("Fix G(0) 2? "+temp);
			IJ.log("G(0) 2 = "+(float)params[5]);
			temp = (fixes[6]==1) ? "true" : "false";
			IJ.log("Fix D2? "+temp);
			IJ.log("D2 = "+(float)params[6]);
			return false;
		}
		c2test=gd.getNextBoolean();
		maxiter = (int)gd.getNextNumber();
		for(int i=0;i<7;i++){
			params[i]=gd.getNextNumber();
			fixes[i]=gd.getNextBoolean() ? 1 : 0;
		}
		return true;
	}
	
	public double fitfunc(double[] params,int indvar)
	{
		//the params list is w0,z0/w0,baseline,g01,D1,g02,D2
		if(indvar<xskip){return g0val[indvar];}
		else{
			int j=(int)(indvar%xpts);
			int i=(int)((indvar-j)/xpts);
			double mult=4.0*(double)photons;
			double xpixels=(double)j;
			double xdistance=xpixels*pixelsize;
			double ypixels=(double)i;
			double ydistance=ypixels*pixelsize;
			double tau=pixeltime*xpixels+linetime*ypixels;
			double sqr_distance=xdistance*xdistance+ydistance*ydistance;
			double w02=params[0]*params[0];
			double z02=params[0]*params[1]*params[0]*params[1];
			double dumdouble1=params[3]/(1.0+((mult*params[4]*tau)/w02));
			if(psfindex==0){
				dumdouble1/=Math.sqrt(1.0+((mult*params[4]*tau)/z02));
			}
			//note that the exponent differs from digman et al by a factor of 2
			dumdouble1*=Math.exp(((double)photons)*(((-sqr_distance)/w02)/(1.0+((mult*params[4]*tau)/w02))));
			double dumdouble=dumdouble1;
			dumdouble1=params[5]/(1.0+((mult*params[6]*tau)/w02));
			if(psfindex==0){
				dumdouble1/=Math.sqrt(1.0+((mult*params[6]*tau)/z02));
			}
			dumdouble1*=Math.exp(((double)photons)*(((-sqr_distance)/w02)/(1.0+((mult*params[6]*tau)/w02))));
			return dumdouble1+dumdouble+params[2];
		}
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
