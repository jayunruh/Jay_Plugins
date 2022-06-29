/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class support_plane_errors_v2 implements NLLSfitinterface_v2{
	NLLSfitinterface_v2 fitclass;
	Object fitfunctionclass;
	double toler;
	int maxiter;
	boolean global;

	public support_plane_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,boolean global1,double lambda1){
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

	public support_plane_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,int maxiter1,boolean global1){
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

	public support_plane_errors_v2(NLLSfitinterface_v2 fitinterface,int maxiter1,boolean global1){
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

	public support_plane_errors_v2(NLLSfitinterface_v2 fitinterface,double toler1,boolean global1){
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

	public support_plane_errors_v2(NLLSfitinterface_v2 fitinterface,boolean global1){
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

	public double[][] geterrors(double[] params,int[] fixes,double[][] constraints,float[] data,float[] weights,double flim,double spacing,int errindex){
		double[] minparams=new double[params.length];
		System.arraycopy(params,0,minparams,0,params.length);
		double thisparam=minparams[errindex];
		int oldnfit=0;
		for(int i=0;i<fixes.length;i++)
			if(fixes[i]==0)
				oldnfit++;
		int newnfit=oldnfit-1;
		double corr=(double)(data.length-newnfit)/(double)(data.length-oldnfit);
		int oldfix=fixes[errindex];
		fixes[errindex]=1;
		double[] stats=new double[2];
		((NLLSfit_v2)fitfunctionclass).fitdata(params,fixes,constraints,data,weights,stats,false);
		double minc2=stats[1]*corr;
		showresults(""+thisparam+" , "+minc2);
		double[][] tempc2plot=new double[2][10000];
		tempc2plot[0][0]=thisparam;
		tempc2plot[1][0]=minc2;
		// start by finding the upper parameter limit
		int nupper=0;
		boolean searching=true;
		double uppererrval=0.0;
		double fprev=1.0;
		double f=1.0;
		do{
			thisparam+=spacing;
			nupper++;
			params[errindex]=thisparam;
			((NLLSfit_v2)fitfunctionclass).fitdata(params,fixes,constraints,data,weights,stats,false);
			fprev=f;
			f=stats[1]*corr/minc2;
			tempc2plot[0][nupper]=thisparam;
			tempc2plot[1][nupper]=stats[1]*corr;
			showresults(""+thisparam+" , "+stats[1]*corr);
			if(f>=flim){
				uppererrval=((flim-fprev)/(f-fprev))*spacing+thisparam-spacing;
				searching=false;
			}else{
				if(constraints!=null&&thisparam>=constraints[1][errindex]){
					uppererrval=constraints[1][errindex];
					searching=false;
				}else{
					if(nupper>=4999||Double.isNaN(f)||Double.isInfinite(f)){
						uppererrval=thisparam;
						searching=false;
					}
				}
			}
		}while(searching);
		// reset everything
		System.arraycopy(minparams,0,params,0,params.length);
		// now find the lower parameter limit
		int nlower=-1;
		searching=true;
		double lowererrval=0.0;
		fprev=1.0;
		f=1.0;
		thisparam=minparams[errindex];
		do{
			thisparam-=spacing;
			nlower++;
			params[errindex]=thisparam;
			((NLLSfit_v2)fitfunctionclass).fitdata(params,fixes,constraints,data,weights,stats,false);
			fprev=f;
			f=stats[1]*corr/minc2;
			tempc2plot[0][nlower+nupper+1]=thisparam;
			tempc2plot[1][nlower+nupper+1]=stats[1]*corr;
			showresults(""+thisparam+" , "+stats[1]*corr);
			if(f>=flim){
				lowererrval=thisparam+spacing-((flim-fprev)/(f-fprev))*spacing;
				searching=false;
			}else{
				if(constraints!=null&&thisparam<=constraints[0][errindex]){
					lowererrval=constraints[0][errindex];
					searching=false;
				}else{
					if(nlower>=4999||Double.isNaN(f)||Double.isInfinite(f)){
						lowererrval=thisparam;
						searching=false;
					}
				}
			}
		}while(searching);
		double[][] tempc2plot2=new double[2][nupper+nlower+3];
		for(int i=0;i<=nlower;i++){
			tempc2plot2[0][i+1]=tempc2plot[0][nupper+nlower+1-i];
			tempc2plot2[1][i+1]=tempc2plot[1][nupper+nlower+1-i];
		}
		for(int i=0;i<=nupper;i++){
			tempc2plot2[0][i+nlower+2]=tempc2plot[0][i];
			tempc2plot2[1][i+nlower+2]=tempc2plot[1][i];
		}
		// reset everything
		System.arraycopy(minparams,0,params,0,params.length);
		fixes[errindex]=oldfix;
		tempc2plot2[0][0]=lowererrval;
		tempc2plot2[1][0]=uppererrval;
		return tempc2plot2;
	}

	public double[][] geterrorsglobal(double[][] params,int[][] vflmatrix,String[][] formulas,String[] paramsnames,double[][][] constraints,float[][] data,float[][] weights,double flim,
			double spacing,int[] errindex){
		double[][] minparams=new double[params.length][params[0].length];
		for(int i=0;i<params.length;i++){
			System.arraycopy(params[i],0,minparams[i],0,params[0].length);
		}
		double thisparam=params[errindex[1]][errindex[0]];
		int oldfix=vflmatrix[errindex[1]][errindex[0]];
		vflmatrix[errindex[1]][errindex[0]]=1;
		double[] stats=new double[2];
		((NLLSglobalfit_v2)fitfunctionclass).fitdata(params,vflmatrix,formulas,paramsnames,constraints,data,weights,stats,new double[params.length],false);
		double minc2=stats[1];
		double[][] tempc2plot=new double[2][10000];
		tempc2plot[0][0]=thisparam;
		tempc2plot[1][0]=minc2;
		// start by finding the upper parameter limit
		int nupper=0;
		boolean searching=true;
		double uppererrval=0.0;
		double fprev=1.0;
		double f=1.0;
		do{
			thisparam+=spacing;
			nupper++;
			params[errindex[1]][errindex[0]]=thisparam;
			((NLLSglobalfit_v2)fitfunctionclass).fitdata(params,vflmatrix,formulas,paramsnames,constraints,data,weights,stats,new double[params.length],false);
			fprev=f;
			f=stats[1]/minc2;
			tempc2plot[0][nupper]=thisparam;
			tempc2plot[1][nupper]=stats[1];
			showresults(""+thisparam+" , "+stats[1]);
			if(f>=flim){
				uppererrval=((flim-fprev)/(f-fprev))*spacing+thisparam-spacing;
				searching=false;
			}else{
				if(constraints!=null&&thisparam>=constraints[1][errindex[1]][errindex[0]]){
					uppererrval=constraints[1][errindex[1]][errindex[0]];
					searching=false;
				}else{
					if(nupper>=4999){
						uppererrval=thisparam;
						searching=false;
					}
				}
			}
		}while(searching);
		// reset everything
		for(int i=0;i<minparams.length;i++){
			System.arraycopy(minparams[i],0,params[i],0,minparams[i].length);
		}
		// now find the lower parameter limit
		int nlower=-1;
		searching=true;
		double lowererrval=0.0;
		fprev=1.0;
		f=1.0;
		thisparam=minparams[errindex[1]][errindex[0]];
		do{
			thisparam-=spacing;
			nlower++;
			params[errindex[1]][errindex[0]]=thisparam;
			((NLLSglobalfit_v2)fitfunctionclass).fitdata(params,vflmatrix,formulas,paramsnames,constraints,data,weights,stats,new double[params.length],false);
			fprev=f;
			f=stats[1]/minc2;
			tempc2plot[0][nlower+nupper+1]=thisparam;
			tempc2plot[1][nlower+nupper+1]=stats[1];
			showresults(""+thisparam+" , "+stats[1]);
			if(f>=flim){
				lowererrval=thisparam+spacing-((flim-fprev)/(f-fprev))*spacing;
				searching=false;
			}else{
				if(constraints!=null&&thisparam<=constraints[0][errindex[1]][errindex[0]]){
					lowererrval=constraints[0][errindex[1]][errindex[0]];
					searching=false;
				}else{
					if(nlower>=4999){
						lowererrval=thisparam;
						searching=false;
					}
				}
			}
		}while(searching);
		double[][] tempc2plot2=new double[2][nupper+nlower+3];
		for(int i=0;i<=nlower;i++){
			tempc2plot2[0][i+1]=tempc2plot[0][nupper+nlower+1-i];
			tempc2plot2[1][i+1]=tempc2plot[1][nupper+nlower+1-i];
		}
		for(int i=0;i<=nupper;i++){
			tempc2plot2[0][i+nlower+2]=tempc2plot[0][i];
			tempc2plot2[1][i+nlower+2]=tempc2plot[1][i];
		}
		// reset everything
		for(int i=0;i<minparams.length;i++){
			System.arraycopy(minparams[i],0,params[i],0,minparams[i].length);
		}
		vflmatrix[errindex[1]][errindex[0]]=oldfix;
		tempc2plot2[0][0]=lowererrval;
		tempc2plot2[1][0]=uppererrval;
		return tempc2plot2;
	}

	public double[] fitfunc(double[] params){
		return fitclass.fitfunc(params);
	}

	public void showresults(String results){
		fitclass.showresults(results);
	}

}
