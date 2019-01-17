/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class gridfit{
	NLLSfitinterface_v2 fitclass;
	public boolean output;
	public boolean escapePressed;

	/*
	 * This is a generic grid search fitting class Any calling class must either
	 * implement the NLLSfitinterface or have access to a class that implements
	 * that interface the interface defined fitfunc method is used to let this
	 * class know what the fitting function is Copyright Jay Unruh Stowers
	 * Institute for Medical Research 4/25/08
	 */

	public gridfit(NLLSfitinterface_v2 fitinterface){
		fitclass=fitinterface;
		output=true;
	}

	public float[] fitdata(double[] params,int[] fixes1,double[][] constraints,float[] data,float[] weights1,double[] stats){
		// this function fits the array data to an arbitrary function denoted by
		// the NLLSfitinterface
		// the returned array is the fit
		// stats returns the chisquared
		// the weights are set to 1 if weights1 is null
		// if fixes[i] is zero, the ith parameter is fit, if not, the ith
		// parameter is fixed
		// constraints is a 3 by nparams array containing the upper and lower
		// bounds for each parameter
		// and the parameter shift for the grid search
		int nparams=params.length;
		int npts=data.length;
		int fitparams=nparams;
		escapePressed=false;
		int[] fixes=new int[nparams];
		if(fixes1!=null){
			for(int i=0;i<nparams;i++){
				fitparams-=fixes1[i];
				fixes[i]=fixes1[i];
			}
		}
		for(int i=0;i<nparams;i++){
			if(fixes[i]==0){
				if(Math.abs(constraints[1][i]-constraints[0][i])<constraints[2][i]){
					fixes[i]=1;
					params[i]=0.5*(constraints[0][i]+constraints[1][i]);
					fitparams--;
				}
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
		stats[1]=search_param(0,fitparams,params,fixes,constraints,data,weights);
		StringBuffer resstring=new StringBuffer();
		resstring.append("final minc2 ="+(float)stats[1]);
		for(int i=0;i<params.length;i++){
			resstring.append(" , "+(float)params[i]);
		}
		if(output){
			fitclass.showresults(resstring.toString());
		}
		double[] fit=fitclass.fitfunc(params);
		float[] ffit=new float[fit.length];
		for(int i=0;i<fit.length;i++){
			ffit[i]=(float)fit[i];
		}
		return ffit;
	}

	public double search_param(int paramid,int numfit,double[] params,int[] fixes,double[][] constraints,float[] data,double[] weights){
		// this method recursively searches over parameter space starting with paramid and returns the
		// minimum parameter and chi squared
		// should think about multithreaded implementation
		// should think about premature breaks
		if(paramid==(params.length-1)){
			if(fixes[paramid]==1){
				return calculate_c2_params(params,numfit,data,weights);
			}else{
				params[paramid]=constraints[1][paramid];
				double minc2=calculate_c2_params(params,numfit,data,weights);
				double minparam=params[paramid];
				for(double x=constraints[1][paramid]-constraints[2][paramid];x>=constraints[0][paramid];x-=constraints[2][paramid]){
					params[paramid]=x;
					double tempc2=calculate_c2_params(params,numfit,data,weights);
					if(tempc2>0.0){
						if(tempc2<minc2||minc2==0.0f){
							minparam=params[paramid];
							minc2=tempc2;
						}
					}
					if(escapePressed) break;
				}
				params[paramid]=minparam;
				return minc2;
			}
		}else{
			int counter=1;
			if(fixes[paramid]==1){
				return search_param(paramid+1,numfit,params,fixes,constraints,data,weights);
			}else{
				params[paramid]=constraints[1][paramid];
				double minc2=search_param(paramid+1,numfit,params,fixes,constraints,data,weights);
				double[] minparams=params.clone();
				if(paramid<2){
					StringBuffer resstring=new StringBuffer();
					resstring.append(""+paramid+" , "+counter+" , "+(float)minc2);
					for(int i=0;i<params.length;i++){
						resstring.append(" , "+(float)params[i]);
					}
					if(output){
						fitclass.showresults(resstring.toString());
					}
				}
				for(double x=constraints[1][paramid]-constraints[2][paramid];x>=constraints[0][paramid];x-=constraints[2][paramid]){
					params[paramid]=x;
					double tempc2=search_param(paramid+1,numfit,params,fixes,constraints,data,weights);
					if(tempc2>0.0){
						if(tempc2<minc2||minc2==0.0f){
							minc2=tempc2;
							minparams=params.clone();
						}
					}
					counter++;
					if(paramid<2){
						StringBuffer resstring=new StringBuffer();
						resstring.append(""+paramid+" , "+counter+" , "+(float)tempc2);
						for(int i=0;i<params.length;i++){
							resstring.append(" , "+(float)params[i]);
						}
						if(output){
							fitclass.showresults(resstring.toString());
						}
					}
					if(escapePressed) break;
				}
				System.arraycopy(minparams,0,params,0,params.length);
				return minc2;
			}
		}
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
			double tempval=(fit[i]-data[i])*(fit[i]-data[i])*weights[i];
			if(!Double.isNaN(tempval)){
				tempc2+=(fit[i]-data[i])*(fit[i]-data[i])*weights[i];
			}
		}
		return tempc2/(length-numfit);
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
