/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.matrixsolve;

public class NLLSfit_v3{
	NLLSfitinterface_v3 fitclass;
	public double toler;
	public int maxiter;
	double dx=0.00001;
	public double lambda;

	/*
	 * This is a generic nonlinear least squares fitting class Any calling class
	 * must either implement the NLLSfitinterface or have access to a class that
	 * implements that interface the interface defined fitfunc method is used to
	 * let this class know what the fitting function is Copyright Jay Unruh
	 * Stowers Institute for Medical Research 4/25/08
	 */

	/************************
	 * this is the main constructor
	 * @param fitinterface: some class implementing the fit interface: usually the calling class
	 * @param toler1: the fractional convergence tolerance, usually 0.0001
	 * @param maxiter1: the max iterations, usually 50
	 * @param lambda1: the levenberg marquardt multiplier, usually 0.1, 0.0 is pure gaus-newton optimization
	 */
	public NLLSfit_v3(NLLSfitinterface_v3 fitinterface,double toler1,int maxiter1,double lambda1){
		// lambda of 0.1 seems to work pretty well
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=lambda1;
	}

	public NLLSfit_v3(NLLSfitinterface_v3 fitinterface,double toler1,int maxiter1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=0.0;
	}

	public NLLSfit_v3(NLLSfitinterface_v3 fitinterface,int maxiter1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		toler=0.0001;
		lambda=0.0;
	}

	public NLLSfit_v3(NLLSfitinterface_v3 fitinterface,double toler1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=50;
		lambda=0.0;
	}

	public NLLSfit_v3(NLLSfitinterface_v3 fitinterface){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=50;
		lambda=0.0;
	}

	/***********************************
	 * this is the main fitting class
	 * @param params: the initial parameter values
	 * @param fixes1: an int array containing 0 for unfixed parameters and 1 for fixed parameters
	 * @param data: the data array
	 * @param weights1: the weights array
	 * @param stats: a double array of length 2 which will contain number of iterations and chi squared on return
	 * @param output: a boolean to specify whether or not data is written to the showresults method of the fit interface
	 * @return
	 */
	public float[] fitdata(double[] params,int[] fixes1,float[] data,float[] weights1,double[] stats,boolean output){
		//TODO--need to implement getting derivatives through interface and add a static derivative function to help out
		// this function fits the array data to an arbitrary function denoted by
		// the NLLSfitinterface
		// the returned array is the fit
		// output dictates whether or not a running tally of the iterations and
		// their chisquareds are desired
		// stats returns the number of iterations and the chi squared
		// the weights are set to 1 if weights1 is null
		// if fixes[i] is zero, the ith parameter is fit, if not, the ith
		// parameter is fixed
		// constraints is a 2 by nparams array containing the upper and lower
		// bounds for each parameter
		int nparams=params.length;
		int npts=data.length;
		int fitparams=nparams;
		double currlambda=lambda;
		int[] fixes=new int[nparams];
		if(fixes1!=null){
			for(int i=0;i<nparams;i++){
				fitparams-=fixes1[i];
				fixes[i]=fixes1[i];
			}
		}
		double[] weights=new double[npts];
		if(weights1==null){
			for(int i=0;i<npts;i++){
				weights[i]=1.0;
			}
		}else{
			for(int i=0;i<npts;i++){
				weights[i]=weights1[i];
			}
		}
		double chisquared=0.0f;
		double c2old=0.0f;
		int iterations=0;
		float[] fit=new float[npts];
		if(maxiter==0){
			double[] dfit=fitclass.fitfunc(params);
			for(int i=0;i<fit.length;i++){
				fit[i]=(float)dfit[i];
			}
			chisquared=calculate_c2_fit(dfit,fitparams,data,weights);
		}else{
			double[] dparams=new double[fitparams];
			double[] dfit=fitclass.fitfunc(params);
			chisquared=calculate_c2_fit(dfit,fitparams,data,weights);
			double tempdouble=0.0;
			int overcount=0;
			do{
				c2old=chisquared;
				int counter=1;
				//calculate the shifted parameters for the numerical derivatives
				double[][] ders=fitclass.derivfunc(params,fixes,dfit);
				//here we calculate the jacobian with derivatives
				double[][] jacobian=new double[fitparams][fitparams];
				double[] jvector=new double[fitparams];
				int[] map=new int[fitparams];
				int counter1=0;
				for(int i=0;i<nparams;i++){
					if(fixes[i]==0){map[counter1]=i; counter1++;}
				}
				for(int i=0;i<fitparams;i++){
					for(int j=0;j<=i;j++){
						for(int k=0;k<npts;k++){
							jacobian[i][j]+=ders[k][map[i]]*ders[k][map[j]]*weights[k];
						}
						if(i!=j){
							jacobian[j][i]=jacobian[i][j];
						}
					}

					for(int k=0;k<npts;k++){
						jvector[i]+=ders[k][map[i]]*(data[k]-dfit[k])*weights[k];
					}
				}
				//implement the levenberg marquardt lagrangian multiplier
				for(int k=0;k<fitparams;k++){
					jacobian[k][k]*=(1.0+currlambda);
				}
				//solve the matrix equation
				(new matrixsolve()).gjsolve(jacobian,jvector,dparams,fitparams);
				//now implement the updates and constraints
				//parameters are truncated at the constraint boundary
				counter=0;
				for(int i=0;i<nparams;i++){
					if(fixes[i]==0){
						params[i]+=dparams[counter];
						counter++;
					}
				}
				fitclass.applyconstraints(params,fixes);
				dfit=fitclass.fitfunc(params);
				chisquared=calculate_c2_fit(dfit,fitparams,data,weights);
				iterations++;
				if(output){
					fitclass.showresults("iteration "+iterations+" c2 = "+chisquared);
				}
				if(iterations==maxiter){
					break;
				}
				tempdouble=(c2old-chisquared)/chisquared;
				if(tempdouble>=0.0&&currlambda>=0.0){
					currlambda/=10.0;
					overcount=0;
				}else{
					currlambda*=10.0;
					if(tempdouble<0.0) overcount++;
					if(overcount>3) break;
				}
			}while(tempdouble>toler||c2old<chisquared);
			stats[0]=iterations;
			for(int i=0;i<npts;i++){
				fit[i]=(float)dfit[i];
			}
		}
		stats[1]=chisquared;
		if(stats.length>2){
			stats[2]=calculate_aic_c2(chisquared,fitparams,data.length);
		}
		return fit;
	}
	
	/************
	 * this function provides the derivatives by numerical integration
	 * 
	 */
	public static double[][] derivfunc2(double[] params,int[] fixes,double[] fit,NLLSfitinterface_v3 fitclass,double dx){
		double[][] shiftfit=new double[params.length][];
		double[][] ders=new double[fit.length][params.length];
		for(int i=0;i<params.length;i++){
			if(fixes[i]==0){
				double[] temp=params.clone();
				temp[i]+=dx;
				shiftfit[i]=fitclass.fitfunc(temp);
			}
		}
		for(int i=0;i<fit.length;i++){
			for(int j=0;j<params.length;j++){
				if(fixes[j]==0){
					ders[i][j]=(shiftfit[j][i]-fit[i])/dx;
				}
			}
		}
		return ders;
	}

	public double calculate_aic_c2(double c2,int numfit,int size1){
		double size=size1;
		double dof=numfit;
		return size*Math.log(c2*(size-dof)/size)+2.0*size*dof/(size-dof-1);
	}
	
	public double calculate_aic_params(double[] params,int numfit,float[] data,double[] weights){
		double c2=calculate_c2_params(params,numfit,data,weights);
		return calculate_aic_c2(c2,numfit,data.length);
	}

	public double calculate_c2_params(double[] params,int numfit,float[] data,double[] weights){
		// this function calculates chisquared using the parameter array
		double[] tempfit=fitclass.fitfunc(params);
		return calculate_c2_fit(tempfit,numfit,data,weights);
	}

	public double calculate_c2_fit(double[] fit,int numfit,float[] data,double[] weights){
		// this function calculates chisquared using a fit array
		int length=data.length;
		double tempc2=0.0;
		for(int i=0;i<length;i++){
			tempc2+=(fit[i]-data[i])*(fit[i]-data[i])*weights[i];
		}
		return tempc2/(length-numfit);
	}
	
	public double calculate_aic_fit(double[] fit,int numfit,float[] data,double[] weights){
		double c2=calculate_c2_fit(fit,numfit,data,weights);
		return calculate_aic_c2(c2,numfit,data.length);
	}

	public double calculate_c2_params(double[] params,int numfit,float[] data,double[] weights,boolean[] fitmask){
		// this function calculates chisquared using the parameter array
		double[] tempfit=fitclass.fitfunc(params);
		return calculate_c2_fit(tempfit,numfit,data,weights,fitmask);
	}

	public double calculate_c2_fit(double[] fit,int numfit,float[] data,double[] weights,boolean[] fitmask){
		// this function calculates chisquared using a fit array
		int length=data.length;
		double tempc2=0.0;
		int nonfit=0;
		for(int i=0;i<length;i++){
			if(!fitmask[i]){
				tempc2+=(fit[i]-data[i])*(fit[i]-data[i])*weights[i];
			}else{
				nonfit++;
			}
		}
		return tempc2/(length-nonfit-numfit);
	}
}
