/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class NLLSglobalfit implements NLLSfitinterface{
	NLLSfitinterface fitclass;
	NLLSfit singlefit;
	public double toler;
	public int maxiter;
	int fitlength,nparams;
	int[][] linkmap;
	public Object[][] fps;
	public double lambda;

	public NLLSglobalfit(NLLSfitinterface fitinterface,double toler1,int maxiter1,double lambda1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		singlefit=new NLLSfit(this,toler,maxiter,lambda1);
		lambda=lambda1;
	}

	public NLLSglobalfit(NLLSfitinterface fitinterface,double toler1,int maxiter1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=maxiter1;
		singlefit=new NLLSfit(this,toler,maxiter);
		lambda=0.0;
	}

	public NLLSglobalfit(NLLSfitinterface fitinterface,int maxiter1){
		fitclass=fitinterface;
		maxiter=maxiter1;
		toler=0.0001;
		singlefit=new NLLSfit(this,toler,maxiter);
		lambda=0.0;
	}

	public NLLSglobalfit(NLLSfitinterface fitinterface,double toler1){
		fitclass=fitinterface;
		toler=toler1;
		maxiter=50;
		singlefit=new NLLSfit(this,toler,maxiter);
		lambda=0.0;
	}

	public NLLSglobalfit(NLLSfitinterface fitinterface){
		fitclass=fitinterface;
		toler=0.0001;
		maxiter=50;
		singlefit=new NLLSfit(this,toler,maxiter);
		lambda=0.0;
	}

	public void changemaxiter(int newmaxiter){
		maxiter=newmaxiter;
		singlefit.maxiter=newmaxiter;
	}

	public float[][] fitdata(double[][] params,int[][] vflmatrix,String[][] formulas,String[] paramsnames,double[][][] constraints,float[][] data,float[][] weights,double[] stats,double[] c2vals,
			boolean output){
		// here we have to map two dimensional arrays to one dimensional arrays
		// with an appropriate linking scheme
		// the resulting 1D arrays are then sent to NLLSfit for optimization
		// note that all fits need to be the same length
		// start by mapping the data and weights
		fitlength=data[0].length;
		int nfits=data.length;
		fps=new Object[nfits][fitlength];
		float[] newdata=new float[nfits*fitlength];
		float[] newweights=new float[nfits*fitlength];
		linkmap=new int[nfits][fitlength];
		for(int i=0;i<nfits;i++){
			for(int j=0;j<fitlength;j++){
				newdata[j+i*fitlength]=data[i][j];
				if(weights!=null){
					newweights[j+i*fitlength]=weights[i][j];
				}else{
					newweights[j+i*fitlength]=1.0f;
				}
			}
		}
		// now find out how many parameters aren't linked
		nparams=params[0].length;
		int nfitparams=0;
		for(int i=0;i<nfits;i++){
			for(int j=0;j<nparams;j++){
				if(i==0){
					// don't allow the leftmost column to be linked (can be
					// functionally linked)
					if(vflmatrix[i][j]==2){
						vflmatrix[i][j]=0;
					}
				}
				if(vflmatrix[i][j]<2){
					nfitparams++;
				}
			}
		}
		double[] newparams=new double[nfitparams];
		int[] newfixes=new int[nfitparams];
		// now map the old parameters to the new parameter list
		int paramindex=0;
		for(int i=0;i<nfits;i++){
			for(int j=0;j<nparams;j++){
				if(vflmatrix[i][j]<2){
					newparams[paramindex]=params[i][j];
					linkmap[i][j]=paramindex;
					if(vflmatrix[i][j]==1){
						newfixes[paramindex]=1;
					}
					paramindex++;
				}else{
					// this one is linked--find the closest leftward non-linked
					// parameter to link it to
					if(vflmatrix[i][j]==2){
						for(int k=i-1;k>=0;k--){
							if(vflmatrix[k][j]!=2){
								linkmap[i][j]=linkmap[k][j];
								if(linkmap[i][j]<0){
									fps[i][j]=fps[k][j];
								}
								break;
							}
						}
					}else{
						// this one is vertically and functionally
						// linked--create a formula parser for it
						linkmap[i][j]=-1;
						fps[i][j]=new formula_parser(formulas[i][j],paramsnames);
					}
				}
			}
		}
		// now map the constraints--constraints on non-linked parameters
		// override constraints on linked parameters
		double[][] newconstraints;
		if(constraints!=null){
			newconstraints=new double[2][nfitparams];
			int counter2=0;
			for(int i=0;i<nfits;i++){
				for(int j=0;j<nparams;j++){
					if(vflmatrix[i][j]<2){
						newconstraints[0][counter2]=constraints[0][i][j];
						newconstraints[1][counter2]=constraints[1][i][j];
						counter2++;
					}
				}
			}
		}else{
			newconstraints=null;
		}
		float[] newfit=singlefit.fitdata(newparams,newfixes,newconstraints,newdata,newweights,stats,output);
		float[][] fit=new float[nfits][fitlength];
		for(int i=0;i<nfits;i++){
			for(int j=0;j<fitlength;j++){
				fit[i][j]=newfit[i*fitlength+j];
			}
		}
		for(int i=0;i<nfits;i++){
			for(int j=0;j<nparams;j++){
				if(linkmap[i][j]>=0){
					params[i][j]=newparams[linkmap[i][j]];
				}else{
					params[i][j]=0.0;
				}
			}
		}
		for(int i=0;i<nfits;i++){
			int numfit=0;
			for(int j=0;j<nparams;j++){
				if(vflmatrix[i][j]!=1){
					numfit++;
				}
			}
			double[] tempweights=new double[fitlength];
			double[] tempfit=new double[fitlength];
			for(int j=0;j<fitlength;j++){
				tempfit[j]=fit[i][j];
				if(weights!=null){
					tempweights[j]=weights[i][j];
				}else{
					tempweights[j]=1.0f;
				}
			}
			c2vals[i]=singlefit.calculate_c2_fit(tempfit,numfit,data[i],tempweights);
		}
		return fit;
	}

	public double fitfunc(double[] params,int indvar){
		// here we map the parameters to the appropriate two dimensional
		// parameter array
		// then based on the value of indvar, we send one column of that array
		// to the calling routine
		int fitindex=(int)((float)indvar/(float)fitlength);
		int newindvar=indvar-fitlength*fitindex;
		double[] newparams=new double[nparams];
		for(int i=0;i<nparams;i++){
			if(linkmap[fitindex][i]>=0){
				newparams[i]=params[linkmap[fitindex][i]];
			}
		}
		for(int i=0;i<nparams;i++){
			if(linkmap[fitindex][i]<0){
				newparams[i]=((formula_parser)fps[fitindex][i]).getfunc(newparams);
			}
		}
		return fitclass.fitfunc(newparams,newindvar);
	}

	public void showresults(String results){
		fitclass.showresults(results);
	}

}
