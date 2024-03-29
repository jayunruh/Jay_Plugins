/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.algutils;
import jalgs.matrixsolve;

import java.util.Arrays;
import java.util.List;

public class linleastsquares{
	public float[][] indvars;
	public double[][] dindvars;
	public int nindvars,npts,startfit,endfit;
	public boolean isfloat;

	/*
	 * This class fits a set of data using linear least squares The fit must be
	 * a linear combination of the independent variable arrays passed to the
	 * constructor e.g.
	 * fit(x)=indvars[0][x]*c[0]+indvars[1][x]*c[1]+indvars[2][x]*c[2]+... Once
	 * the object is constructed it can be used to fit any data with the same
	 * number of data points and the same independent variables for example,
	 * this can be used for linear unmixing Copyright Jay Unruh Stowers
	 * Institute for Medical Research 4/25/08
	 */

	public linleastsquares(){
		// only use this constructor for the get_amp_offset method
	}
	
	public linleastsquares(List<List<Float>> indvars1){
		nindvars=indvars1.size();
		indvars=new float[nindvars][];
		for(int i=0;i<nindvars;i++){
			indvars[i]=algutils.convert_arr_float(indvars1.get(i));
		}
		npts=indvars[0].length;
		startfit=0;
		endfit=npts-1;
		isfloat=true;
	}

	public linleastsquares(float[][] indvars1){
		// here we provide the component independent variable arrays
		nindvars=indvars1.length;
		npts=indvars1[0].length;
		indvars=indvars1;
		startfit=0;
		endfit=npts-1;
		isfloat=true;
	}

	public linleastsquares(double[][] indvars1){
		// here we provide the component independent variable arrays
		nindvars=indvars1.length;
		npts=indvars1[0].length;
		dindvars=indvars1;
		startfit=0;
		endfit=npts-1;
		isfloat=false;
	}

	public linleastsquares(double[][] indvars1,boolean baseline,int startfit1,int endfit1){
		// here we provide the component independent variable arrays
		// this version can add a baseline array and allows for a fit range
		if(!baseline){
			nindvars=indvars1.length;
			npts=indvars1[0].length;
			dindvars=indvars1;
		}else{
			nindvars=indvars1.length+1;
			npts=indvars1[0].length;
			dindvars=new double[nindvars+1][];
			dindvars[0]=new double[npts];
			Arrays.fill(dindvars[0],1.0);
			for(int i=1;i<nindvars;i++){
				dindvars[i]=indvars1[i-1];
			}
		}
		startfit=startfit1;
		endfit=endfit1;
		isfloat=false;
	}

	public linleastsquares(float[][] indvars1,boolean baseline,int startfit1,int endfit1){
		// here we provide the component independent variable arrays
		// this version adds a baseline array
		if(!baseline){
			nindvars=indvars1.length;
			npts=indvars1[0].length;
			indvars=indvars1;
		}else{
			nindvars=indvars1.length+1;
			npts=indvars1[0].length;
			indvars=new float[nindvars+1][];
			indvars[0]=new float[npts];
			Arrays.fill(indvars[0],1.0f);
			for(int i=1;i<nindvars;i++){
				indvars[i]=indvars1[i-1];
			}
		}
		startfit=startfit1;
		endfit=endfit1;
		isfloat=true;
	}
	
	public double[] fitdata(List<Float> data,List<Float> weights){
		return fitdata(algutils.convert_arr_float(data),algutils.convert_arr_float(weights));
	}

	public double[] fitdata(float[] data,float[] weights){
		double[][] jacobian=new double[nindvars][nindvars];
		double[] jvector=new double[nindvars];
		if(isfloat){
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(weights!=null){
							jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k]*weights[k];
						}else{
							jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(weights!=null){
						jvector[i]+=(double)weights[k]*(double)indvars[i][k]*data[k];
					}else{
						jvector[i]+=(double)data[k]*(double)indvars[i][k];
					}
				}
			}
		}else{
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(weights!=null){
							jacobian[i][j]+=dindvars[i][k]*dindvars[j][k]*weights[k];
						}else{
							jacobian[i][j]+=dindvars[i][k]*dindvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(weights!=null){
						jvector[i]+=weights[k]*dindvars[i][k]*data[k];
					}else{
						jvector[i]+=data[k]*dindvars[i][k];
					}
				}
			}
		}
		double[] outcoef=new double[nindvars];
		(new matrixsolve()).gjsolve(jacobian,jvector,outcoef,nindvars);
		return outcoef;
	}
	
	/*************************
	 * here we get the fit coefficients as well as their standard errors
	 * @param data
	 * @return
	 */
	public double[][] getfiterrors(float[] data,float[] weights){
		double[][] jacobian=new double[nindvars][nindvars];
		double[] jvector=new double[nindvars];
		if(isfloat){
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(weights!=null){
							jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k]*weights[k];
						}else{
							jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(weights!=null){
						jvector[i]+=(double)weights[k]*(double)indvars[i][k]*data[k];
					}else{
						jvector[i]+=(double)data[k]*(double)indvars[i][k];
					}
				}
			}
		}else{
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(weights!=null){
							jacobian[i][j]+=dindvars[i][k]*dindvars[j][k]*weights[k];
						}else{
							jacobian[i][j]+=dindvars[i][k]*dindvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(weights!=null){
						jvector[i]+=weights[k]*dindvars[i][k]*data[k];
					}else{
						jvector[i]+=data[k]*dindvars[i][k];
					}
				}
			}
		}
		//the coefficient variances are long the diagonal of the inversion matrix
		double[][] inv=(new matrixsolve()).gjinv2(jacobian,nindvars);
		double[] outcoef=matrixsolve.vec_mult(inv,jvector);
		double[] se=new double[nindvars];
		double c2=get_c2(outcoef,data,weights);
		//int length=endfit-startfit+1;
		//c2*=(double)(length-nindvars)/(double)(length-1);
		for(int i=0;i<nindvars;i++) se[i]=Math.sqrt(c2*inv[i][i]);
		return new double[][]{outcoef,se};
	}

	public double[] fitintensitydata(float[] data){
		// here we assume that variance = 1/intensity
		double[][] jacobian=new double[nindvars][nindvars];
		double[] jvector=new double[nindvars];
		if(isfloat){
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(data[k]>0.0f){
							jacobian[i][j]+=((double)indvars[i][k]*(double)indvars[j][k])/data[k];
						}else{
							jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(data[k]>0.0f){
						jvector[i]+=indvars[i][k];
					}else{
						jvector[i]+=(double)data[k]*(double)indvars[i][k];
					}
				}
			}
		}else{
			for(int i=0;i<nindvars;i++){
				for(int j=0;j<=i;j++){
					for(int k=startfit;k<=endfit;k++){
						if(data[k]>0.0f){
							jacobian[i][j]+=(dindvars[i][k]*dindvars[j][k])/data[k];
						}else{
							jacobian[i][j]+=dindvars[i][k]*dindvars[j][k];
						}
					}
					if(j!=i){
						jacobian[j][i]=jacobian[i][j];
					}
				}
				for(int k=startfit;k<=endfit;k++){
					if(data[k]>0.0f){
						jvector[i]+=dindvars[i][k];
					}else{
						jvector[i]+=data[k]*dindvars[i][k];
					}
				}
			}
		}
		double[] outcoef=new double[nindvars];
		(new matrixsolve()).gjsolve(jacobian,jvector,outcoef,nindvars);
		return outcoef;
	}

	/*********
	 * gets the slope for a linear fit with the first point being 0 and increasing by 1
	 * @param data
	 * @return
	 */
	public float get_slope(float[] data){
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		for(int i=0;i<data.length;i++){
			sumx2+=i*i;
			sumx+=i;
			sumy+=data[i];
			sumxy+=data[i]*i;
		}
		double dlength=data.length;
		if(sumx2>0.0){
			double divider=dlength*sumx2-sumx*sumx;
			return (float)((dlength*sumxy-sumx*sumy)/divider);
		}else{
			return 0.0f;
		}
	}
	
	/**********************
	 * gets the slope and offset for a linear fit as with get_slope
	 * @param data
	 * @return
	 */
	public float[] get_slope_offset(float[] data){
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		for(int i=0;i<data.length;i++){
			sumx2+=(double)i*(double)i;
			sumx+=i;
			sumy+=data[i];
			sumxy+=data[i]*(double)i;
		}
		double dlength=data.length;
		double off,amp;
		if(sumx2>0.0){
			double divider=dlength*sumx2-sumx*sumx;
			off=(sumx2*sumy-sumx*sumxy)/divider;
			amp=(dlength*sumxy-sumx*sumy)/divider;
		}else{
			amp=0.0;
			off=sumy/dlength;
		}
		float[] fitparams={(float)amp,(float)off};
		return fitparams;
	}

	public float[] get_slope_se(float[] data){
		// here we get the slope and the standard error in the slope
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		for(int i=0;i<data.length;i++){
			sumx2+=i*i;
			sumx+=i;
			sumy+=data[i];
			sumxy+=data[i]*i;
		}
		double dlength=data.length;
		if(sumx2>0.0){
			double divider=dlength*sumx2-sumx*sumx;
			double slope=(dlength*sumxy-sumx*sumy)/divider;
			double avgx=sumx/dlength;
			double avgy=sumy/dlength;
			double offset=avgy-slope*avgx;
			double c2=0.0f;
			for(int i=0;i<data.length;i++)
				c2+=(data[i]-offset-slope*i)*(data[i]-offset-slope*i);
			double se=Math.sqrt(c2/((dlength-2.0)*(sumx2-sumx*sumx/dlength)));
			return new float[]{(float)slope,(float)se};
		}else{
			return null;
		}
	}

	/*****************************************
	 * this function gets the amplitude and offset for any function (fit[i]=amp*function[i]+offset)
	 * @param function
	 * @param data: the data to fit
	 * @param offset: whether or not to fit the offset
	 * @return
	 */
	public float[] get_amp_offset(float[] function,float[] data,boolean offset){
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		for(int i=0;i<data.length;i++){
			sumx2+=function[i]*function[i];
			sumx+=function[i];
			sumy+=data[i];
			sumxy+=data[i]*function[i];
		}
		if(offset){
			double dlength=data.length;
			double off,amp;
			if(sumx2>0.0){
				double divider=dlength*sumx2-sumx*sumx;
				off=(sumx2*sumy-sumx*sumxy)/divider;
				amp=(dlength*sumxy-sumx*sumy)/divider;
			}else{
				amp=0.0;
				off=sumy/dlength;
			}
			float[] fitparams={(float)amp,(float)off};
			return fitparams;
		}else{
			double amp;
			if(sumx2>0.0){
				amp=sumxy/sumx2;
			}else{
				amp=0.0;
			}
			float[] fitparams={(float)amp};
			return fitparams;
		}
	}

	public double[] get_amp_offset(double[] function,float[] data,boolean offset){
		return get_amp_offset(function,data,offset,0,data.length-1);
	}
	
	public double[] get_amp_offset(double[] function,float[] data,boolean offset,int start,int end) {
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		for(int i=start;i<=end;i++){
			sumx2+=function[i]*function[i];
			sumx+=function[i];
			sumy+=data[i];
			sumxy+=data[i]*function[i];
		}
		if(offset){
			double dlength=end-start+1;
			double off,amp;
			if(sumx2>0.0){
				double divider=dlength*sumx2-sumx*sumx;
				off=(sumx2*sumy-sumx*sumxy)/divider;
				amp=(dlength*sumxy-sumx*sumy)/divider;
			}else{
				amp=0.0;
				off=sumy/dlength;
			}
			double[] fitparams={amp,off};
			return fitparams;
		}else{
			double amp;
			if(sumx2>0.0){
				amp=sumxy/sumx2;
			}else{
				amp=0.0;
			}
			double[] fitparams={amp};
			return fitparams;
		}
	}
	
	/****************
	 * this is like loess smooth except with a user defined function--find the amp and offset at windowsize around each point
	 * @param function: the function
	 * @param data
	 * @param offset
	 * @param windowsize
	 * @return an array of coefficients (first one or two places) and the fit
	 */
	public double[][] get_running_amp_offset(double[] function,float[] data,boolean offset,int windowsize){
		int halfwsize=windowsize/2;
		double sumx2=0.0,sumxy=0.0,sumx=0.0,sumy=0.0;
		double[][] profile=new double[data.length][];
		for(int j=0;j<windowsize;j++) {
			sumx2+=function[j]*function[j];
			sumx+=function[j];
			sumy+=data[j];
			sumxy+=data[j]*function[j];
		}
		for(int i=halfwsize;i<(data.length-halfwsize);i++) {
			if(i>halfwsize) {
				int prev=i-halfwsize-1;
				int next=i-halfwsize+windowsize-1;
				sumx2-=function[prev]*function[prev]; sumx2+=function[next]*function[next];
				sumx-=function[prev]; sumx+=function[next];
				sumy-=data[prev]; sumy+=data[next];
				sumxy-=data[prev]*function[prev]; sumxy+=data[next]*function[next];
			}
			if(offset){
				double dlength=windowsize;
				double off,amp;
				if(sumx2>0.0){
					double divider=dlength*sumx2-sumx*sumx;
					off=(sumx2*sumy-sumx*sumxy)/divider;
					amp=(dlength*sumxy-sumx*sumy)/divider;
				}else{
					amp=0.0;
					off=sumy/dlength;
				}
				profile[i]=new double[]{amp,off,amp*function[i]+off};
			}else{
				double amp;
				if(sumx2>0.0){
					amp=sumxy/sumx2;
				}else{
					amp=0.0;
				}
				profile[i]=new double[]{amp,amp*function[i]};
			}
		}
		//now fill in the start and end with nearest values
		for(int i=0;i<halfwsize;i++) profile[i]=profile[halfwsize].clone();
		for(int i=(data.length-halfwsize);i<data.length;i++) profile[i]=profile[data.length-halfwsize-1].clone();
		return profile;
	}
	
	public double get_amp_offset_c2(double[] function,float[] data,double[] coef){
		double c2=0.0;
		if(coef.length>1){
			for(int i=0;i<function.length;i++){
				float resid=((float)coef[0]*(float)function[i]+(float)coef[1])-data[i];
				c2+=resid*resid;
			}
			c2/=(function.length-2);
		}else{
			for(int i=0;i<function.length;i++){
				float resid=(float)coef[0]*(float)function[i]-data[i];
				c2+=resid*resid;
			}
			c2/=(function.length-1);
		}
		return c2;
	}

	public double get_amp_offset_c2(float[] function,float[] data,float[] coef){
		double c2=0.0;
		if(coef.length>1){
			for(int i=0;i<function.length;i++){
				float resid=(coef[0]*function[i]+coef[1])-data[i];
				c2+=resid*resid;
			}
			c2/=(function.length-2);
		}else{
			for(int i=0;i<function.length;i++){
				float resid=coef[0]*function[i]-data[i];
				c2+=resid*resid;
			}
			c2/=(function.length-1);
		}
		return c2;
	}

	public double[] get_fit(double[] coef){
		double[] fit=new double[npts];
		if(isfloat){
			for(int i=0;i<npts;i++){
				for(int j=0;j<nindvars;j++){
					fit[i]+=coef[j]*indvars[j][i];
				}
			}
		}else{
			for(int i=0;i<npts;i++){
				for(int j=0;j<nindvars;j++){
					fit[i]+=coef[j]*dindvars[j][i];
				}
			}
		}
		return fit;
	}

	public float[] get_ffit(double[] coef){
		float[] fit=new float[npts];
		if(isfloat){
			for(int i=0;i<npts;i++){
				for(int j=0;j<nindvars;j++){
					fit[i]+=(float)(coef[j]*indvars[j][i]);
				}
			}
		}else{
			for(int i=0;i<npts;i++){
				for(int j=0;j<nindvars;j++){
					fit[i]+=(float)(coef[j]*dindvars[j][i]);
				}
			}
		}
		return fit;
	}

	public double get_c2(double[] coef,float[] data,float[] weights){
		int length=endfit-startfit+1;
		double tempc2=0.0;
		double[] fit=get_fit(coef);
		if(weights!=null){
			for(int i=startfit;i<=endfit;i++){
				tempc2+=(fit[i]-data[i])*(fit[i]-data[i])*weights[i];
			}
		}else{
			for(int i=startfit;i<=endfit;i++){
				tempc2+=(fit[i]-data[i])*(fit[i]-data[i]);
			}
		}
		return tempc2/(length-nindvars);
	}

	public double get_int_c2(double[] coef,float[] data){
		int length=endfit-startfit+1;
		double tempc2=0.0;
		double[] resid=get_int_resid(coef,data);
		for(int i=startfit;i<=endfit;i++){
			tempc2+=resid[i];
		}
		return tempc2/(length-nindvars);
	}

	public double[] get_resid(double[] coef,float[] data,float[] weights){
		int length=data.length;
		double[] resid=get_fit(coef);
		if(weights!=null){
			for(int i=0;i<length;i++){
				resid[i]-=data[i];
				resid[i]*=Math.sqrt(weights[i]);
			}
		}else{
			for(int i=0;i<length;i++){
				resid[i]-=data[i];
			}
		}
		return resid;
	}

	public float[] get_fresid(double[] coef,float[] data,float[] weights){
		int length=data.length;
		double[] fit=get_fit(coef);
		float[] resid=new float[length];
		if(weights!=null){
			for(int i=0;i<length;i++){
				resid[i]=(float)((fit[i]-data[i])*Math.sqrt(weights[i]));
			}
		}else{
			for(int i=0;i<length;i++){
				resid[i]=(float)(fit[i]-data[i]);
			}
		}
		return resid;
	}

	public double[] get_int_resid(double[] coef,float[] data){
		int length=data.length;
		double[] resid=get_fit(coef);
		for(int i=0;i<length;i++){
			resid[i]-=data[i];
			if(data[i]>0.0f){
				resid[i]*=Math.sqrt(1.0/data[i]);
			}
		}
		return resid;
	}

	public float[] get_int_fresid(double[] coef,float[] data){
		int length=data.length;
		double[] fit=get_fit(coef);
		float[] resid=new float[length];
		for(int i=0;i<length;i++){
			if(data[i]>0.0f){
				resid[i]=(float)((fit[i]-data[i])*Math.sqrt(1.0/data[i]));
			}else{
				resid[i]=(float)(fit[i]-data[i]);
			}
		}
		return resid;
	}

}
