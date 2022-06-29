/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class simplexfit_v2{
	NLLSfitinterface_v2 fitclass;
	double toler,simplex_size;
	int maxiter,miniter;

	/*
	 * This is a simplex fitting class Any calling class must either implement
	 * the NLLSfitinterface or have access to a class that implements that
	 * interface the interface defined fitfunc method is used to let this class
	 * know what the fitting function is Copyright Jay Unruh Stowers Institute
	 * for Medical Research 2/26/09
	 */

	public simplexfit_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,int miniter1,double simplex_size1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		miniter=miniter1;
		simplex_size=simplex_size1;
	}

	public simplexfit_v2(NLLSfitinterface_v2 fitinterface,int maxiter1,int miniter1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		miniter=miniter1;
		toler=0.0001;
		simplex_size=0.1;
	}

	public simplexfit_v2(NLLSfitinterface_v2 fitinterface,double toler1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=100;
		miniter=1;
		simplex_size=0.1;
	}

	public simplexfit_v2(NLLSfitinterface_v2 fitinterface){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=100;
		miniter=1;
		simplex_size=0.1;
	}

	public float[] fitdata(double[] params,int[] fixes1,double[][] constraints,float[] data,float[] weights1,double[] stats,boolean output){
		int nparams=params.length;
		int nfit=0;
		boolean[] fixes=new boolean[nparams];
		boolean[] fits=new boolean[nparams];
		for(int i=0;i<nparams;i++){
			if(fixes1[i]==0){
				nfit++;
				fixes[i]=false;
				fits[i]=true;
			}else{
				fixes[i]=true;
				fits[i]=false;
			}
		}
		double[][] simplex=new double[nfit][nfit+1];
		double[] c2=new double[nfit+1];
		double[] weights=new double[data.length];
		if(weights1==null){
			for(int i=0;i<data.length;i++){
				weights[i]=1.0;
			}
		}else{
			for(int i=0;i<data.length;i++){
				weights[i]=weights1[i];
			}
		}
		if(maxiter==0){
			stats[0]=0.0;
			stats[1]=calculate_c2_params(params,nfit,data,weights,constraints);
			double[] fit=fitclass.fitfunc(params);
			float[] ffit=new float[fit.length];
			for(int i=0;i<fit.length;i++){
				ffit[i]=(float)fit[i];
			}
			return ffit;
		}
		int counter=0;
		for(int i=0;i<nparams;i++){
			if(fits[i]){
				simplex[counter][0]=params[i];
				counter++;
			}
		}

		// now generate the simplex
		for(int i=0;i<nfit;i++){
			for(int j=0;j<nfit;j++){
				if(i==j){
					simplex[j][i+1]=simplex[j][0]+simplex_size*simplex[j][0];
				}else{
					simplex[j][i+1]=simplex[j][0];
				}
			}
		}
		// calculate the chi squared at each of the vertices
		for(int i=0;i<=nfit;i++){
			double[] tempparams=new double[nparams];
			counter=0;
			for(int j=0;j<nparams;j++){
				if(fits[j]){
					tempparams[j]=simplex[counter][i];
					counter++;
				}else{
					tempparams[j]=params[j];
				}
			}
			c2[i]=calculate_c2_params(tempparams,nfit,data,weights,constraints);
		}
		double c2old=c2[0];
		// find out which index is the best
		int best_index=0;
		for(int i=1;i<=nfit;i++){
			if(c2[i]<c2[best_index]){
				best_index=i;
			}
		}
		int iteration=1;
		int worst_index;
		boolean first=true;
		do{
			if(!first){
				iteration++;
				c2old=c2[best_index];
			}
			// find out which simplex is the worst
			worst_index=0;
			for(int i=1;i<=nfit;i++){
				if(c2[i]>c2[worst_index]){
					worst_index=i;
				}
			}
			// calculate the average of all simplexes but the worst
			double[] avg_all_but_worst=new double[nfit];
			for(int i=0;i<=nfit;i++){
				if(i!=worst_index){
					for(int j=0;j<nfit;j++){
						avg_all_but_worst[j]=avg_all_but_worst[j]+simplex[j][i]/nfit;
					}
				}
			}
			// first try a reflection
			counter=0;
			for(int i=0;i<nparams;i++){
				if(fits[i]){
					params[i]=2.0*avg_all_but_worst[counter]-simplex[counter][worst_index];
					counter++;
				}
			}
			double testnew=calculate_c2_params(params,nfit,data,weights,constraints);
			if(testnew<c2[worst_index]){
				// if it worked, keep it
				counter=0;
				for(int i=0;i<nparams;i++){
					if(fits[i]){
						simplex[counter][worst_index]=params[i];
						counter++;
					}
				}
				c2[worst_index]=testnew;
				if(c2[worst_index]<c2[best_index]){
					// if the worst point is now the best point, try an
					// expansion along that direction
					// first save the information
					best_index=worst_index;
					counter=0;
					for(int i=0;i<nparams;i++){
						if(fits[i]){
							params[i]=2.0*simplex[counter][worst_index]-avg_all_but_worst[counter];
							counter++;
						}
					}
					testnew=calculate_c2_params(params,nfit,data,weights,constraints);
					if(testnew<c2[best_index]){
						counter=0;
						for(int i=0;i<nparams;i++){
							if(fits[i]){
								simplex[counter][best_index]=params[i];
								counter++;
							}
						}
						c2[best_index]=testnew;
					}
				}
			}else{
				// if the reflection didn't work, try a contraction towards the
				// center
				counter=0;
				for(int i=0;i<nparams;i++){
					if(fits[i]){
						params[i]=0.5*(simplex[counter][worst_index]+avg_all_but_worst[counter]);
						counter++;
					}
				}
				testnew=calculate_c2_params(params,nfit,data,weights,constraints);
				if(testnew<c2[worst_index]){
					// if it worked save it
					counter=0;
					for(int i=0;i<nparams;i++){
						if(fits[i]){
							simplex[counter][worst_index]=params[i];
							counter++;
						}
					}
					c2[worst_index]=testnew;
					if(testnew<c2[best_index]){
						best_index=worst_index;
					}
				}else{
					// if the contraction didn't work, try an overall
					// contraction towards the current best vertex
					int temp_index=best_index;
					for(int i=0;i<=nfit;i++){
						if(i!=temp_index){
							counter=0;
							for(int j=0;j<nparams;j++){
								if(fits[j]){
									simplex[counter][i]=0.5*(simplex[counter][i]+simplex[counter][temp_index]);
									params[j]=simplex[counter][i];
									counter++;
								}
							}
							testnew=calculate_c2_params(params,nfit,data,weights,constraints);
							if(testnew<c2[best_index]){
								best_index=i;
								c2[best_index]=testnew;
							}
						}
					}
				}
			}
			// save the best vertex and report the chi squared
			counter=0;
			for(int i=0;i<nparams;i++){
				if(fits[i]){
					params[i]=simplex[counter][best_index];
					counter++;
				}
			}
			if(output){
				fitclass.showresults("Iteration "+iteration+" c2 = "+(float)c2[best_index]);
			}
			first=false;
		}while(((iteration<maxiter)&&(Math.abs((c2[best_index]-c2old)/c2old)>toler))||(iteration<miniter));
		stats[0]=iteration;
		stats[1]=c2[best_index];
		double[] fit=fitclass.fitfunc(params);
		float[] ffit=new float[fit.length];
		for(int i=0;i<fit.length;i++){
			ffit[i]=(float)fit[i];
		}
		return ffit;
	}

	public double calculate_c2_params(double[] params,int numfit,float[] data,double[] weights,double[][] constraints){
		// this function calculates chisquared using the parameter array
		int length=data.length;
		double tempc2=0.0;
		double[] tempfit=fitclass.fitfunc(params);
		for(int i=0;i<length;i++){
			tempc2+=(tempfit[i]-data[i])*(tempfit[i]-data[i])*weights[i];
		}
		// implement penalties for going outside of constraints
		for(int i=0;i<params.length;i++){
			if(params[i]<constraints[0][i]){
				tempc2*=10;
			}else{
				if(params[i]>constraints[1][i]){
					tempc2*=10;
				}
			}
		}
		return tempc2/(length-numfit);
	}

}
