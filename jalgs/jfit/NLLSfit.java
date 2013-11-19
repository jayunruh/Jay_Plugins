/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.matrixsolve;

public class NLLSfit{
	NLLSfitinterface fitclass;
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

	public NLLSfit(NLLSfitinterface fitinterface,double toler1,int maxiter1,double lambda1){
		// lambda of 0.1 seems to work pretty well
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=lambda1;
	}

	public NLLSfit(NLLSfitinterface fitinterface,double toler1,int maxiter1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		lambda=0.0;
	}

	public NLLSfit(NLLSfitinterface fitinterface,int maxiter1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		toler=0.0001;
		lambda=0.0;
	}

	public NLLSfit(NLLSfitinterface fitinterface,double toler1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=50;
		lambda=0.0;
	}

	public NLLSfit(NLLSfitinterface fitinterface){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=50;
		lambda=0.0;
	}

	public float[] fitdata(double[] params,int[] fixes1,double[][] constraints,float[] data,float[] weights1,double[] stats,boolean output){
		// this function fits the array data to an arbitrary function denoted by
		// the NLLSfitinterface
		// the returned array is the fit
		// output dictates whether or not a running tally of the iterations and
		// their chisquareds are desired
		// stats returns the chisquared and the number of variables
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
				weights[i]=(double)weights1[i];
			}
		}
		double chisquared=0.0f;
		double c2old=0.0f;
		int iterations=0;
		float[] fit=new float[npts];
		if(maxiter==0){
			for(int i=0;i<npts;i++){
				fit[i]=(float)fitclass.fitfunc(params,i);
			}
			chisquared=calculate_c2_params(params,fitparams,data,weights);
		}else{
			Object[] shiftparams=new Object[fitparams+1];
			Object[] shiftfit=new Object[fitparams+1];
			double[] dparams=new double[fitparams];
			shiftparams[0]=params;
			shiftfit[0]=new double[npts];
			for(int i=0;i<npts;i++){
				((double[])shiftfit[0])[i]=(double)fitclass.fitfunc((double[])shiftparams[0],i);
			}
			chisquared=calculate_c2_fit((double[])shiftfit[0],fitparams,data,weights);
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
					shiftfit[i]=new double[npts];
					for(int j=0;j<npts;j++){
						((double[])shiftfit[i])[j]=(double)fitclass.fitfunc((double[])shiftparams[i],j);
					}
				}
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
						jvector[i]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*((double)data[k]-((double[])shiftfit[0])[k])*weights[k];
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
				for(int i=0;i<npts;i++){
					((double[])shiftfit[0])[i]=(double)fitclass.fitfunc((double[])shiftparams[0],i);
				}
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
			stats[0]=(double)iterations;
			for(int i=0;i<npts;i++){
				fit[i]=(float)((double[])shiftfit[0])[i];
			}
		}
		stats[1]=chisquared;
		return fit;
	}

	public float[] fitintensitydata(double[] params,int[] fixes1,double[][] constraints,float[] data,double[] stats,boolean output,double S,boolean[] fitmask){
		// this function fits the array data to an arbitrary function denoted by
		// the NLLSfitinterface
		// the returned array is the fit
		// output dictates whether or not a running tally of the iterations and
		// their chisquareds are desired
		// stats returns the chisquared and the number of variables
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
				weights[i]=1.0/(S*(double)data[i]);
			}else{
				weights[i]=1.0;
			}
		}
		double chisquared=0.0f;
		double c2old=0.0f;
		int iterations=0;
		float[] fit=new float[npts];
		if(maxiter==0){
			for(int i=0;i<npts;i++){
				fit[i]=(float)fitclass.fitfunc(params,i);
			}
			chisquared=calculate_c2_params(params,fitparams,data,weights,fitmask);
		}else{
			Object[] shiftparams=new Object[fitparams+1];
			Object[] shiftfit=new Object[fitparams+1];
			double[] dparams=new double[fitparams];
			shiftparams[0]=params;
			shiftfit[0]=new double[npts];
			for(int i=0;i<npts;i++){
				((double[])shiftfit[0])[i]=(double)fitclass.fitfunc((double[])shiftparams[0],i);
			}
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
					shiftfit[i]=new double[npts];
					for(int j=0;j<npts;j++){
						((double[])shiftfit[i])[j]=(double)fitclass.fitfunc((double[])shiftparams[i],j);
					}
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
							jvector[i]+=((((double[])shiftfit[i+1])[k]-((double[])shiftfit[0])[k])/dx)*((double)data[k]-((double[])shiftfit[0])[k])*weights[k];
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
				for(int i=0;i<npts;i++){
					((double[])shiftfit[0])[i]=(double)fitclass.fitfunc((double[])shiftparams[0],i);
				}
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
			stats[0]=(double)iterations;
			for(int i=0;i<npts;i++){
				fit[i]=(float)((double[])shiftfit[0])[i];
			}
		}
		stats[1]=chisquared;
		return fit;
	}

	public double calculate_c2_params(double[] params,int numfit,float[] data,double[] weights){
		// this function calculates chisquared using the parameter array
		int length=data.length;
		double tempc2=0.0;
		for(int i=0;i<length;i++){
			double tempdouble=fitclass.fitfunc(params,i);
			tempc2+=(tempdouble-(double)data[i])*(tempdouble-(double)data[i])*weights[i];
		}
		return tempc2/((double)(length-numfit));
	}

	public double calculate_c2_fit(double[] fit,int numfit,float[] data,double[] weights){
		// this function calculates chisquared using a fit array
		int length=data.length;
		double tempc2=0.0;
		for(int i=0;i<length;i++){
			tempc2+=(fit[i]-(double)data[i])*(fit[i]-(double)data[i])*weights[i];
		}
		return tempc2/((double)(length-numfit));
	}

	public double calculate_c2_params(double[] params,int numfit,float[] data,double[] weights,boolean[] fitmask){
		// this function calculates chisquared using the parameter array
		int length=data.length;
		int nonfit=0;
		double tempc2=0.0;
		for(int i=0;i<length;i++){
			if(!fitmask[i]){
				double tempdouble=fitclass.fitfunc(params,i);
				tempc2+=(tempdouble-(double)data[i])*(tempdouble-(double)data[i])*weights[i];
			}else{
				nonfit++;
			}
		}
		return tempc2/((double)(length-nonfit-numfit));
	}

	public double calculate_c2_fit(double[] fit,int numfit,float[] data,double[] weights,boolean[] fitmask){
		// this function calculates chisquared using a fit array
		int length=data.length;
		double tempc2=0.0;
		int nonfit=0;
		for(int i=0;i<length;i++){
			if(!fitmask[i]){
				tempc2+=(fit[i]-(double)data[i])*(fit[i]-(double)data[i])*weights[i];
			}else{
				nonfit++;
			}
		}
		return tempc2/((double)(length-nonfit-numfit));
	}
}
