/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class fitpoly{
	public linleastsquares lls; // linear least squares class used for the fit
	public float[][] indvars;
	public int npts;
	
	/*****************
	 * here we assume that x starts at 0 and increases linearly
	 * @param order
	 * @param npts
	 * @param fixoff
	 */
	public fitpoly(int order,int npts,boolean fixoff) {
		float[] xvals=new float[npts];
		for(int i=0;i<npts;i++) xvals[i]=(float)i;
		init(order,xvals,fixoff);
	}

	/**********
	 * This class fits a set of data to a polynomial of order "order" Once the
	 * object is constructed it can be used to fit any data with the same number
	 * of data points and the same fit order
	 */
	public fitpoly(int order,float[] xvals,boolean fixoff){
		init(order,xvals,fixoff);
	}
	
	private void init(int order,float[] xvals,boolean fixoff) {
		npts=xvals.length;
		if(fixoff){
			indvars=new float[order][npts];
			for(int i=1;i<=order;i++){
				for(int j=0;j<npts;j++){
					indvars[i-1][j]=(float)Math.pow(xvals[j],i);
				}
			}
			lls=new linleastsquares(indvars);
		}else{
			indvars=new float[order+1][npts];
			for(int j=0;j<npts;j++){
				indvars[0][j]=1.0f;
			}
			for(int i=1;i<=order;i++){
				for(int j=0;j<npts;j++){
					indvars[i][j]=(float)Math.pow(xvals[j],i);
				}
			}
			lls=new linleastsquares(indvars);
		}
	}

	public double[] fitdatapolyd(float[] data,float[] weights){
		if(weights==null){
			return lls.fitdata(data,null);
		}else{
			return lls.fitdata(data,weights);
		}
	}

	public float[] fitdatapoly(float[] data,float[] weights){
		double[] coef=fitdatapolyd(data,weights);
		float[] coef2=new float[coef.length];
		for(int i=0;i<coef.length;i++){
			coef2[i]=(float)coef[i];
		}
		return coef2;
	}

	public float[] getfit(double[] coefficients){
		float[] fit=new float[npts];
		for(int i=0;i<npts;i++){
			for(int j=0;j<coefficients.length;j++){
				fit[i]+=coefficients[j]*indvars[j][i];
			}
		}
		return fit;
	}
	
	public float getfitpt(double[] coefficients,int pt){
		float temp=0.0f;
		for(int i=0;i<coefficients.length;i++){
			temp+=coefficients[i]*indvars[i][pt];
		}
		return temp;
	}

	public float[] getfit(float[] coefficients){
		float[] fit=new float[npts];
		for(int i=0;i<npts;i++){
			for(int j=0;j<coefficients.length;j++){
				fit[i]+=coefficients[j]*indvars[j][i];
			}
		}
		return fit;
	}
	
	public float getfitpt(float[] coefficients,int pt){
		float temp=0.0f;
		for(int i=0;i<coefficients.length;i++){
			temp+=coefficients[i]*indvars[i][pt];
		}
		return temp;
	}

}
