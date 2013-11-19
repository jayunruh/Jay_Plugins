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

public class fit_roi_gaussian_jru_v1 implements PlugIn, NLLSfitinterface {
	int xpts;

	public void run(String arg) {

		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageProcessor ip=imp.getProcessor();
		Rectangle r=ip.getRoi();
		FloatProcessor fp=null;
		if(ip instanceof FloatProcessor){
			fp=(FloatProcessor)ip;
		} else {
			fp=(FloatProcessor)ip.convertToFloat();
		}
		float[] pixels=(float[])fp.getPixels();
		float[] roi=new float[r.width*r.height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				roi[j-r.x+(i-r.y)*r.width]=pixels[j+width*i];
			}
		}
		float baseline=roi[0];
		for(int i=0;i<r.width*r.height;i++){
			roi[i]-=baseline;
		}
		double avgx=0.0;
		double avgy=0.0;
		double sum=0.0;
		for(int i=0;i<r.height;i++){
			for(int j=0;j<r.width;j++){
				double temp=roi[j+i*r.width];
				avgx+=temp*(double)j;
				avgy+=temp*(double)i;
				sum+=temp;
			}
		}
		avgx/=sum;
		int intavgx=(int)avgx;
		if(intavgx>=r.width){intavgx=r.width-1;}
		if(intavgx<0){intavgx=0;}
		avgy/=sum;
		int intavgy=(int)avgy;
		if(intavgy>=r.height){intavgy=r.height-1;}
		if(intavgy<0){intavgy=0;}
		float predamp=roi[intavgx+intavgy*r.width];
		float halfmax=predamp/2.0f;
		float predstdev=0.0f;
		for(int i=intavgx;i<r.width;i++){
			if(roi[i+intavgy*r.width]<=halfmax){
				predstdev+=(float)(i-intavgx);
				break;
			}
		}
		for(int i=intavgx;i>=0;i--){
			if(roi[i+intavgy*r.width]<=halfmax){
				predstdev+=(float)(intavgx-i);
				break;
			}
		}
		for(int i=intavgy;i<r.height;i++){
			if(roi[intavgx+i*r.width]<=halfmax){
				predstdev+=(float)(i-intavgy);
				break;
			}
		}
		for(int i=intavgy;i>=0;i--){
			if(roi[intavgx+i*r.width]<=halfmax){
				predstdev+=(float)(intavgy-i);
				break;
			}
		}
		predstdev/=4.0f;
		//IJ.log("pred stdev = "+(float)predstdev+" , amp = "+(float)predamp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("x center",avgx,5,10,null);
		gd.addNumericField("y center",avgy,5,10,null);
		gd.addNumericField("stdev",predstdev,5,10,null);
		gd.addNumericField("amp",predamp,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		avgx=(float)gd.getNextNumber();
		avgy=(float)gd.getNextNumber();
		predstdev=(float)gd.getNextNumber();
		predamp=(float)gd.getNextNumber();
		xpts=r.width;
		NLLSfit fitclass=new NLLSfit(this,0.0001,50,0.1);
		double[] params={0.0,predamp,predstdev,avgx,avgy};
		double[][] constraints=new double[2][5];
		constraints[0][0]=-predamp; constraints[1][0]=predamp;
		constraints[0][1]=0.2*predamp; constraints[1][1]=5.0*predamp;
		constraints[0][2]=0.05; constraints[1][2]=(double)xpts;
		constraints[0][3]=0.0; constraints[1][3]=(double)r.width;
		constraints[0][4]=0.0; constraints[1][4]=(double)r.height;
		double[] stats=new double[2];
		int[] fixes=new int[5];
		boolean[] fitmask=new boolean[r.width*r.height];
		double S=0.0;
		float[] fit=fitclass.fitintensitydata(params,fixes,constraints,roi,stats,true,S,fitmask);
		IJ.log("baseline = "+(float)(params[0]+baseline));
		IJ.log("amplitude = "+(float)params[1]);
		IJ.log("stdev = "+(float)params[2]);
		IJ.log("x center = "+(float)params[3]);
		IJ.log("y center = "+(float)params[4]);
		new ImagePlus("fit",new FloatProcessor(r.width,r.height,fit,null)).show();		
	}

	public double fitfunc(double[] params,int indvar)
	{
		//the params list is baseline,amplitude,stdev,xc,yc
		int j=(int)(indvar%xpts);
		int i=(int)((indvar-j)/xpts);
		double xdiff=((double)j-params[3]);
		double ydiff=((double)i-params[4]);
		double dumdouble;
		dumdouble = Math.exp(((-0.5)*(xdiff*xdiff))/(params[2]*params[2]));
		dumdouble *= Math.exp(((-0.5)*(ydiff*ydiff))/(params[2]*params[2]));
		dumdouble *= params[1];
		return dumdouble+params[0];
	}

	public void showresults(String results){
		IJ.log(results);
	}


}
