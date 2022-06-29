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
//import ij.plugin.filter.*;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;

public class fit_ICS_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	ImagePlus imp;
	double[] params;
	int[] fixes;
	double[] stats;
	boolean g0skip,c2test;
	double pi = 3.14159265;
	int xc,yc,xpts,ypts,maxiter;
	float g0val;
	NLLSfit_v2 fitclass;

	/*public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL;
	}*/

	//public void run(ImageProcessor ip) {
	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageProcessor ip=imp.getProcessor();
		float[] pixels = (float[])ip.getPixels();
		int fullwidth = imp.getWidth();
		int fullheight = imp.getHeight();

		//initialize the variables
		xc=fullwidth/2; yc=fullheight/2; xpts=16; ypts=16; g0skip=true; maxiter=10; c2test=false;
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Skip_G0",true);
		gd2.addNumericField("X center",xc,0);
		gd2.addNumericField("Y center",yc,0);
		gd2.addNumericField("Pts in x to fit",xpts,0);
		gd2.addNumericField("Pts in y to fit",ypts,0);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		g0skip=gd2.getNextBoolean();
		xc=(int)gd2.getNextNumber();
		yc=(int)gd2.getNextNumber();
		xpts=(int)gd2.getNextNumber();
		ypts=(int)gd2.getNextNumber();
		params=new double[5];
		params[1]=1.0; params[2]=3.0; params[4]=5.0;
		fixes=new int[5];
		fixes[3]=fixes[4]=1;
		stats=new double[2];
		double[][] constraints={{-10.0,0.0,0.1,0.0,0.1},{10.0,1000.0,100.0,1000.0,100.0}};
		float[] ac=new float[xpts*ypts];
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				ac[j+i*xpts]=pixels[(i+yc)*fullwidth+j+xc];
			}
		}
		g0val=ac[0];

		ImagePlus[] imps=new ImagePlus[3];
		imps[0] = new ImagePlus("Autocorr",new FloatProcessor(xpts,ypts,ac,null));
		imps[0].show();
		Rectangle r0=imps[0].getWindow().getBounds();
		imps[1] = new ImagePlus("Fit",new FloatProcessor(xpts,ypts,new float[xpts*ypts],null));
		imps[1].show();
		Rectangle r1=imps[1].getWindow().getBounds();
		imps[1].getWindow().setBounds(r0.x,r0.y+r0.height,r1.width,r1.height);
		imps[1].updateAndDraw();
		imps[2] = new ImagePlus("Residuals",new FloatProcessor(xpts,ypts,new float[xpts*ypts],null));
		imps[2].show();
		Rectangle r2=imps[2].getWindow().getBounds();
		imps[2].getWindow().setBounds(r0.x,r0.y+2*r0.height,r2.width,r2.height);
		imps[2].updateAndDraw();

		String[] paramnames={"Basline","G0_1","StDev_Particle1","G0_2","StDev_Particle2"};
		FitDialog_image fd=new FitDialog_image(imps,this,paramnames);
		fd.run_fit(params,null,constraints,fixes);

		/*while (showDialog())
		{
			NLLSfit nf;
			if(c2test){
				nf=new NLLSfit(this,0);				
			}
			else{
				nf=new NLLSfit(this,maxiter);
			}
			float[] fit=nf.fitdata(params,fixes,constraints,ac,null,stats,false);
			float[] resid=new float[xpts*ypts];
			for(int i=0;i<xpts*ypts;i++){
				resid[i]=ac[i]-fit[i];
			}
			imp2.setProcessor(null,new FloatProcessor(xpts,ypts,resid,null));
			if(imp2.getWindow()==null){imp2.show();}
			imp3.setProcessor(null,new FloatProcessor(xpts,ypts,fit,null));
			if(imp3.getWindow()==null){imp3.show();}
			imp4.setProcessor(null,new FloatProcessor(xpts,ypts,ac,null));
			if(imp4.getWindow()==null){imp4.show();}
		}*/
	}

	/*private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog("Fitting Parameters");
		gd.addCheckbox("Skip G(0)?",g0skip);
		gd.addCheckbox("Test Chi Squared?",c2test);
		gd.addNumericField("Max Iterations?",maxiter,0);
		gd.addNumericField("Baseline",params[0],10,15,"");
		gd.addCheckbox("fix?",(fixes[0]==1));
		gd.addNumericField("G(0) 1",params[1],10,15,"");
		gd.addCheckbox("fix?",(fixes[1]==1));
		gd.addNumericField("StDev Particle 1",params[2],10,15,"");
		gd.addCheckbox("fix?",(fixes[2]==1));
		gd.addNumericField("G(0) 2",params[3],10,15,"");
		gd.addCheckbox("fix?",(fixes[3]==1));
		gd.addNumericField("StDev Particle 2",params[4],10,15,"");
		gd.addCheckbox("fix?",(fixes[4]==1));
		gd.addNumericField("Iterations Completed",(int)stats[0],0);
		gd.addNumericField("Chi Squared",(float)stats[1],5,15,"");
		gd.showDialog();
		if(gd.wasCanceled()){
			String temp;
			IJ.log("Fit Results");
			temp = g0skip ? "true" : "false";
			IJ.log("Skip G(0) ? "+temp);
			IJ.log("Points in x to fit "+xpts);
			IJ.log("Points in y to fit "+ypts);
			temp = (fixes[0]==1) ? "true" : "false";
			IJ.log("Fix Baseline? "+temp);
			IJ.log("Baseline = "+(float)params[0]);
			temp = (fixes[1]==1) ? "true" : "false";
			IJ.log("Fix G(0) 1? "+temp);
			IJ.log("G(0) 1 = "+(float)params[1]);
			temp = (fixes[2]==1) ? "true" : "false";
			IJ.log("Fix StDev1? "+temp);
			IJ.log("StDev Particle 1 = "+(float)params[2]);
			temp = (fixes[3]==1) ? "true" : "false";
			IJ.log("Fix G(0) 2? "+temp);
			IJ.log("G(0) 2 = "+(float)params[3]);
			temp = (fixes[4]==1) ? "true" : "false";
			IJ.log("Fix StDev2? "+temp);
			IJ.log("StDev Particle 2 = "+(float)params[4]);
			return false;
		}
		g0skip=gd.getNextBoolean();
		c2test=gd.getNextBoolean();
		maxiter = (int)gd.getNextNumber();
		params[0] = gd.getNextNumber();
		fixes[0]=gd.getNextBoolean() ? 1 : 0;
		params[1]= gd.getNextNumber();
		fixes[1]=gd.getNextBoolean() ? 1 : 0;
		params[2] = gd.getNextNumber();
		fixes[2]=gd.getNextBoolean() ? 1 : 0;
		params[3] = gd.getNextNumber();
		fixes[3]=gd.getNextBoolean() ? 1 : 0;
		params[4] = gd.getNextNumber();
		fixes[4]=gd.getNextBoolean() ? 1 : 0;
		return true;
	}*/
	
	/*public double fitfunc(double[] params,int indvar)
	{
		//the params list is baseline,g01,stdev1,g02,stdev2
		if(indvar==0 && g0skip){return g0val;}
		else{
			int j=(int)(indvar%xpts);
			int i=(int)((indvar-j)/xpts);
			double dumdouble,dumdouble1;
			dumdouble1 = Math.exp(((-0.25)*(double)(j*j))/(params[2]*params[2]));
			dumdouble1 *= Math.exp(((-0.25)*(double)(i*i))/(params[2]*params[2]));
			dumdouble1 *= params[1];
			dumdouble=dumdouble1;
			dumdouble1 = Math.exp(((-0.25)*(double)(j*j))/(params[4]*params[4]));
			dumdouble1 *= Math.exp(((-0.25)*(double)(i*i))/(params[4]*params[4]));
			dumdouble1 *= params[3];
			return dumdouble1+dumdouble+params[0];
		}
	}*/

	public double[] fitfunc(double[] params){
		//the params list is baseline,g01,stdev1,g02,stdev2
		double[] func=new double[xpts*ypts];
		if(g0skip){func[0]=g0val;}
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				if(j==0 && i==0 && g0skip) continue;
				double dumdouble,dumdouble1;
				dumdouble1 = Math.exp(((-0.25)*(double)(j*j))/(params[2]*params[2]));
				dumdouble1 *= Math.exp(((-0.25)*(double)(i*i))/(params[2]*params[2]));
				dumdouble1 *= params[1];
				dumdouble=dumdouble1;
				dumdouble1 = Math.exp(((-0.25)*(double)(j*j))/(params[4]*params[4]));
				dumdouble1 *= Math.exp(((-0.25)*(double)(i*i))/(params[4]*params[4]));
				dumdouble1 *= params[3];
				func[j+i*xpts]=dumdouble1+dumdouble+params[0];
			}
		}
		return func;
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
