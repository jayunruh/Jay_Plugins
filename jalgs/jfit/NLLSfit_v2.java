/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.matrixsolve;

public class NLLSfit_v2{
	NLLSfitinterface_v2 fitclass;
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
	public NLLSfit_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,double lambda1){
		// lambda of 0.1 seems to work pretty well
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=lambda1;
	}

	public NLLSfit_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=0.0;
	}

	public NLLSfit_v2(NLLSfitinterface_v2 fitinterface,int maxiter1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		toler=0.0001;
		lambda=0.0;
	}

	public NLLSfit_v2(NLLSfitinterface_v2 fitinterface,double toler1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=50;
		lambda=0.0;
	}

	public NLLSfit_v2(NLLSfitinterface_v2 fitinterface){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=50;
		lambda=0.0;
	}

	/***********************************
	 * this is the main fitting function
	 * @param params: the initial parameter values
	 * @param fixes1: an int array containing 0 for unfixed parameters and 1 for fixed parameters
	 * @param constraints: an nparams x 2 2D array with the constraints: {{lower0,lower1,...},{upper0,upper1,...}}
	 * @param data: the data array
	 * @param weights1: the weights array
	 * @param stats: a double array of length 2 which will contain number of iterations and chi squared on return
	 * @param output: a boolean to specify whether or not data is written to the showresults method of the fit interface
	 * @return
	 */
	public float[] fitdata(double[] params,int[] fixes1,double[][] constraints,float[] data,float[] weights1,double[] stats,boolean output){
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
			Object[] shiftparams=new Object[fitparams+1];
			Object[] shiftfit=new Object[fitparams+1];
			double[] dparams=new double[fitparams];
			shiftparams[0]=params;
			shiftfit[0]=fitclass.fitfunc((double[])shiftparams[0]);
			chisquared=calculate_c2_fit((double[])shiftfit[0],fitparams,data,weights);
			double tempdouble=0.0;
			do{
				c2old=chisquared;
				int counter=1;
				//calculate the shifted parameters for the numerical derivatives
				for(int i=0;i<nparams;i++){
					if(fixes[i]==0){
						shiftparams[counter]=new double[nparams];
						for(int j=0;j<nparams;j++){
							((double[])shiftparams[counter])[j]=((double[])shiftparams[0])[j];
						}
						((double[])shiftparams[counter])[i]+=dx;
						counter++;
					}
				}
				//now calculate the shifted functions
				for(int i=1;i<=fitparams;i++){
					shiftfit[i]=fitclass.fitfunc((double[])shiftparams[i]);
				}
				//here we calculate the jacobian with derivatives
				double[][] jacobian=new double[fitparams][fitparams];
				double[] jvector=new double[fitparams];
				for(int i=0;i<fitparams;i++){
					for(int j=0;j<=i;j++){
						for(int k=0;k<npts;k++){
							jacobian[i][j]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*((((double[])shiftfit[j+1])[k]-((double[])shiftfit[0])[k])/dx)*weights[k];
						}
						if(i!=j){
							jacobian[j][i]=jacobian[i][j];
						}
					}

					for(int k=0;k<npts;k++){
						jvector[i]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*(data[k]-((double[])shiftfit[0])[k])*weights[k];
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
						((double[])shiftparams[0])[i]+=dparams[counter];
						if(constraints!=null){
							if(((double[])shiftparams[0])[i]<constraints[0][i]){
								((double[])shiftparams[0])[i]=constraints[0][i];
							}
							if(((double[])shiftparams[0])[i]>constraints[1][i]){
								((double[])shiftparams[0])[i]=constraints[1][i];
							}
						}
						counter++;
					}
				}
				shiftfit[0]=fitclass.fitfunc((double[])shiftparams[0]);
				chisquared=calculate_c2_fit((double[])shiftfit[0],fitparams,data,weights);
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
				}else{
					currlambda*=10.0;
				}
			}while(tempdouble>toler||c2old<chisquared);
			stats[0]=iterations;
			for(int i=0;i<npts;i++){
				fit[i]=(float)((double[])shiftfit[0])[i];
			}
		}
		stats[1]=chisquared;
		if(stats.length>2){
			stats[2]=calculate_aic_c2(chisquared,fitparams,data.length);
		}
		return fit;
	}

	public float[] fitintensitydata(double[] params,int[] fixes1,double[][] constraints,float[] data,double[] stats,boolean output,double S,boolean[] fitmask){
		// this function fits the array data to an arbitrary function denoted by
		// the NLLSfitinterface
		// the returned array is the fit
		// output dictates whether or not a running tally of the iterations and
		// their chisquareds are desired
		// stats returns the number of iterations and the chi squared
		// if fixes[i] is zero, the ith parameter is fit, if not, the ith
		// parameter is fixed
		// constraints is a 2 by nparams array containing the upper and lower
		// bounds for each parameter
		// this version assumes that var=S*Intensity
		// fitmask is true for pixels that are not to be counted
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
		for(int i=0;i<npts;i++){
			if(data[i]>0.0f&&S>0.0f){
				weights[i]=1.0/(S*data[i]);
			}else{
				weights[i]=1.0;
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
			chisquared=calculate_c2_params(params,fitparams,data,weights,fitmask);
		}else{
			Object[] shiftparams=new Object[fitparams+1];
			Object[] shiftfit=new Object[fitparams+1];
			double[] dparams=new double[fitparams];
			shiftparams[0]=params;
			shiftfit[0]=fitclass.fitfunc((double[])shiftparams[0]);
			chisquared=calculate_c2_fit((double[])shiftfit[0],fitparams,data,weights,fitmask);
			double tempdouble=0.0;
			do{
				c2old=chisquared;
				int counter=1;
				for(int i=0;i<nparams;i++){
					if(fixes[i]==0){
						shiftparams[counter]=new double[nparams];
						for(int j=0;j<nparams;j++){
							((double[])shiftparams[counter])[j]=((double[])shiftparams[0])[j];
						}
						((double[])shiftparams[counter])[i]+=dx;
						counter++;
					}
				}
				for(int i=1;i<=fitparams;i++){
					shiftfit[i]=fitclass.fitfunc((double[])shiftparams[i]);
				}
				double[][] jacobian=new double[fitparams][fitparams];
				double[] jvector=new double[fitparams];
				for(int i=0;i<fitparams;i++){
					for(int j=0;j<=i;j++){
						for(int k=0;k<npts;k++){
							if(!fitmask[k]){
								jacobian[i][j]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*((((double[])shiftfit[j+1])[k]-((double[])shiftfit[0])[k])/dx)*weights[k];
							}
						}
						if(i!=j){
							jacobian[j][i]=jacobian[i][j];
						}
					}

					for(int k=0;k<npts;k++){
						if(!fitmask[k]){
							jvector[i]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*(data[k]-((double[])shiftfit[0])[k])*weights[k];
						}
					}
				}
				for(int k=0;k<fitparams;k++){
					jacobian[k][k]*=(1.0+currlambda);
				}
				(new matrixsolve()).gjsolve(jacobian,jvector,dparams,fitparams);
				counter=0;
				for(int i=0;i<nparams;i++){
					if(fixes[i]==0){
						((double[])shiftparams[0])[i]+=dparams[counter];
						if(constraints!=null){
							if(((double[])shiftparams[0])[i]<constraints[0][counter]){
								((double[])shiftparams[0])[i]=constraints[0][counter];
							}
							if(((double[])shiftparams[0])[i]>constraints[1][counter]){
								((double[])shiftparams[0])[i]=constraints[1][counter];
							}
						}
						counter++;
					}
				}
				shiftfit[0]=fitclass.fitfunc((double[])shiftparams[0]);
				chisquared=calculate_c2_fit((double[])shiftfit[0],fitparams,data,weights,fitmask);
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
				}else{
					currlambda*=10.0;
				}
			}while(tempdouble>toler||c2old<chisquared);
			stats[0]=iterations;
			for(int i=0;i<npts;i++){
				fit[i]=(float)((double[])shiftfit[0])[i];
			}
		}
		stats[1]=chisquared;
		if(stats.length>2){
			int npts2=data.length;
			if(fitmask!=null){
				for(int i=0;i<fitmask.length;i++) if(fitmask[i]) npts2--;
			}
			stats[2]=calculate_aic_c2(chisquared,fitparams,npts2);
		}
		return fit;
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
