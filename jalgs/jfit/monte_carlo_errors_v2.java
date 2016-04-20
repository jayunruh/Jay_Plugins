/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.jsim.rngs;

public class monte_carlo_errors_v2 implements NLLSfitinterface_v2{
	NLLSfitinterface_v2 fitclass;
	Object fitfunctionclass;
	double toler;
	int maxiter;
	boolean global;
	rngs random;

	public monte_carlo_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,boolean global1,double lambda1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		global=global1;
		if(global){
			fitfunctionclass=new NLLSglobalfit_v2(this,toler,maxiter,lambda1);
		}else{
			fitfunctionclass=new NLLSfit_v2(this,toler,maxiter,lambda1);
		}
	}

	public monte_carlo_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,boolean global1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		global=global1;
		if(global){
			fitfunctionclass=new NLLSglobalfit_v2(this,toler,maxiter);
		}else{
			fitfunctionclass=new NLLSfit_v2(this,toler,maxiter);
		}
	}

	public monte_carlo_errors_v2(NLLSfitinterface_v2 fitinterface,int maxiter1,boolean global1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		toler=0.0001;
		global=global1;
		if(global){
			fitfunctionclass=new NLLSglobalfit_v2(this,toler,maxiter);
		}else{
			fitfunctionclass=new NLLSfit_v2(this,toler,maxiter);
		}
	}

	public monte_carlo_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,boolean global1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=50;
		global=global1;
		if(global){
			fitfunctionclass=new NLLSglobalfit_v2(this,toler,maxiter);
		}else{
			fitfunctionclass=new NLLSfit_v2(this,toler,maxiter);
		}
	}

	public monte_carlo_errors_v2(NLLSfitinterface_v2 fitinterface,boolean global1){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=50;
		global=global1;
		if(global){
			fitfunctionclass=new NLLSglobalfit_v2(this,toler,maxiter);
		}else{
			fitfunctionclass=new NLLSfit_v2(this,toler,maxiter);
		}
	}

	/********************************
	 * gets errors by adding gaussian random numbers to a data set (according to the residuals and weights)
	 * outputs a params+1 x ntrials array with fit values and last column with the chi squared
	 * @param params: the fit parameters
	 * @param fixes: with 1 for fixed and 0 for varied
	 * @param constraints: a 2 x nparams array with upper and lower constraints
	 * @param data: the data
	 * @param weights: the weights (can be null)
	 * @param ntrials: number of trials (usually 500)
	 * @return
	 */
	public double[][] geterrors(double[] params,int[] fixes,double[][] constraints,float[] data,float[] weights,int ntrials){
		random=new rngs();
		int oldnfit=0;
		for(int i=0;i<fixes.length;i++)
			if(fixes[i]==0)
				oldnfit++;
		double[] stats=new double[2];
		float[] fit=((NLLSfit_v2)fitfunctionclass).fitdata(params,fixes,constraints,data,weights,stats,false);
		double minc2=stats[1];
		double[][] sampling=new double[oldnfit+1][ntrials];
		for(int i=0;i<ntrials;i++){
			// add random noise with a variance equal to minc2 times weights to
			// the data and then fit again
			float[] sim=fit.clone();
			if(weights!=null){
				for(int j=0;j<fit.length;j++){
					if(weights[j]>0.0f){
						float stdev=(float)Math.sqrt(minc2/weights[j]);
						sim[j]+=random.gasdev(0.0,stdev);
					}
				}
			}else{
				for(int j=0;j<fit.length;j++){
					float stdev=(float)Math.sqrt(minc2);
					sim[j]+=random.gasdev(0.0,stdev);
				}
			}
			double[] simparams=params.clone();
			double[] simstats=new double[2];
			((NLLSfit_v2)fitfunctionclass).fitdata(simparams,fixes,constraints,sim,weights,simstats,false);
			int counter=0;
			StringBuffer sb=new StringBuffer();
			sb.append(""+(i+1)+"\t");
			for(int j=0;j<fixes.length;j++){
				if(fixes[j]==0){
					sampling[counter][i]=simparams[j];
					sb.append(""+simparams[j]+"\t");
					counter++;
				}
			}
			sampling[oldnfit][i]=simstats[1];
			sb.append(""+simstats[1]);
			showresults(sb.toString());
		}
		return sampling;
	}

	public double[][] geterrorsglobal(double[][] params,int[][] vflmatrix,String[][] formulas,String[] paramsnames,double[][][] constraints,float[][] data,float[][] weights,int ntrials){
		double[] stats=new double[2];
		float[][] fit=((NLLSglobalfit_v2)fitfunctionclass).fitdata(params,vflmatrix,formulas,paramsnames,constraints,data,weights,stats,new double[params.length],false);
		double minc2=stats[1];
		int oldnfit=0;
		for(int i=0;i<params.length;i++){
			for(int j=0;j<params[0].length;j++){
				if(vflmatrix[i][j]==0)
					oldnfit++;
			}
		}
		double[][] sampling=new double[oldnfit+1][ntrials];
		for(int i=0;i<ntrials;i++){
			// add random noise with a variance equal to minc2 times weights to
			// the data and then fit again
			float[][] sim=new float[fit.length][];
			for(int j=0;j<fit.length;j++)
				sim[j]=fit[j].clone();
			for(int j=0;j<fit.length;j++){
				for(int k=0;k<fit[0].length;k++){
					if(weights[j][j]>0.0f){
						float stdev=(float)Math.sqrt(minc2/weights[j][k]);
						sim[j][k]+=random.gasdev(0.0,stdev);
					}
				}
			}
			double[][] simparams=new double[params.length][];
			for(int j=0;j<params.length;j++)
				simparams[j]=params[j].clone();
			double[] simstats=new double[2];
			((NLLSglobalfit_v2)fitfunctionclass).fitdata(simparams,vflmatrix,formulas,paramsnames,constraints,sim,weights,simstats,new double[params.length],false);
			int counter=0;
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<params.length;j++){
				for(int k=0;k<params[0].length;k++){
					if(vflmatrix[j][k]==0){
						sampling[counter][i]=simparams[j][k];
						sb.append(""+simparams[j]+" , ");
						counter++;
					}
				}
			}
			sampling[i][oldnfit]=stats[1];
			sb.append(""+stats[1]);
			showresults(sb.toString());
		}
		return sampling;
	}

	public double[] fitfunc(double[] params){
		return fitclass.fitfunc(params);
	}

	public void showresults(String results){
		fitclass.showresults(results);
	}

}
