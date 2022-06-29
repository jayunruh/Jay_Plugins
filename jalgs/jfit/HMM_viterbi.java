/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.large_double;

public class HMM_viterbi{

	public double viterbi(float[] data,double[] states,int[] retpath,double emit_err,double[][] trans_p){
		int nstates=states.length;
		int npts=data.length;
		double[][] V=new double[npts][nstates];
		int[][] path=new int[nstates][npts];
		for(int i=0;i<nstates;i++){
			V[0][i]=emit_p(data[0],emit_err,states[i])/nstates;
			path[i][0]=i;
		}
		for(int t=1;t<npts;t++){
			for(int y=0;y<nstates;y++){
				double maxval=V[t-1][0]*trans_p[0][y];
				int maxindex=0;
				for(int k=1;k<nstates;k++){
					double temp=V[t-1][k]*trans_p[k][y];
					if(temp>maxval){
						maxval=temp;
						maxindex=k;
					}
				}
				V[t][y]=maxval*emit_p(data[t],emit_err,states[y]);
				path[y][t]=maxindex;
			}
		}
		double maxval=V[npts-1][0];
		int maxindex=0;
		for(int i=1;i<nstates;i++){
			if(V[npts-1][i]>maxval){
				maxval=V[npts-1][i];
				maxindex=i;
			}
		}
		retpath[npts-1]=path[maxindex][npts-1];
		for(int t=(npts-2);t>=0;t--){
			retpath[t]=path[retpath[t+1]][t+1];
		}
		return maxval;
	}

	public large_double viterbi_poisson(float[] data,double[] states,double background,int[] retpath,double[][] trans_p,double[] start_p){
		float[] newbackground=new float[data.length];
		for(int i=0;i<data.length;i++){
			newbackground[i]=(float)background;
		}
		return viterbi_poisson(data,states,newbackground,retpath,trans_p,start_p);
	}

	public large_double viterbi_poisson(float[] data,double[] states,float[] background,int[] retpath,double[][] trans_p,double[] start_p){
		int nstates=states.length;
		int npts=data.length;
		large_double[] V=new large_double[nstates];
		int[][] path=new int[nstates][npts];
		if(start_p==null){
			for(int i=0;i<nstates;i++){
				V[i]=new large_double(emit_p_poisson(data[0],states[i],background[0])/nstates);
				path[i][0]=i;
			}
		}else{
			for(int i=0;i<nstates;i++){
				V[i]=new large_double(emit_p_poisson(data[0],states[i],background[0])*start_p[i]);
				path[i][0]=i;
			}
		}
		for(int t=1;t<npts;t++){
			large_double[] Vnew=new large_double[nstates];
			for(int y=0;y<nstates;y++){
				large_double maxval=V[0].mult_copy(trans_p[0][y]);
				int maxindex=0;
				for(int k=1;k<nstates;k++){
					large_double temp=V[k].mult_copy(trans_p[k][y]);
					if(temp.compare(maxval)>0){
						maxval=temp;
						maxindex=k;
					}
				}
				Vnew[y]=maxval.mult_copy(emit_p_poisson(data[t],states[y],background[t]));
				path[y][t]=maxindex;
			}
			V=Vnew;
		}
		large_double maxval=V[0];
		int maxindex=0;
		for(int i=1;i<nstates;i++){
			if(V[i].compare(maxval)>0){
				maxval=V[i];
				maxindex=i;
			}
		}
		retpath[npts-1]=path[maxindex][npts-1];
		for(int t=(npts-2);t>=0;t--){
			retpath[t]=path[retpath[t+1]][t+1];
		}
		return maxval;
	}

	public large_double viterbi2(float[] data,double[] states,float background,float stdev,int[] retpath,double[][] trans_p,double[] start_p){
		int nstates=states.length;
		int npts=data.length;
		large_double[] V=new large_double[nstates];
		int[][] path=new int[nstates][npts];
		if(start_p==null){
			for(int i=0;i<nstates;i++){
				V[i]=new large_double(emit_p((double)data[0]-background,stdev,states[i])/nstates);
				path[i][0]=i;
			}
		}else{
			for(int i=0;i<nstates;i++){
				V[i]=new large_double(emit_p((double)data[0]-background,stdev,states[i])*start_p[i]);
				path[i][0]=i;
			}
		}
		for(int t=1;t<npts;t++){
			large_double[] Vnew=new large_double[nstates];
			for(int y=0;y<nstates;y++){
				large_double maxval=V[0].mult_copy(trans_p[0][y]);
				int maxindex=0;
				for(int k=1;k<nstates;k++){
					large_double temp=V[k].mult_copy(trans_p[k][y]);
					if(temp.compare(maxval)>0){
						maxval=temp;
						maxindex=k;
					}
				}
				Vnew[y]=maxval.mult_copy(emit_p(data[t]-background,stdev,states[y]));
				path[y][t]=maxindex;
			}
			V=Vnew;
		}
		large_double maxval=V[0];
		int maxindex=0;
		for(int i=1;i<nstates;i++){
			if(V[i].compare(maxval)>0){
				maxval=V[i];
				maxindex=i;
			}
		}
		retpath[npts-1]=path[maxindex][npts-1];
		for(int t=(npts-2);t>=0;t--){
			retpath[t]=path[retpath[t+1]][t+1];
		}
		return maxval;
	}

	public double emit_p(double observation,double emit_err,double state){
		// this returns the gaussian emission probability of seeing observation
		// for a given state and emit_err noise stdev
		return Math.exp(-(observation-state)*(observation-state)/(2.0*emit_err*emit_err))/(emit_err*Math.sqrt(2.0*Math.PI));
	}

	public double emit_p_poisson(double observation,double state,double background){
		// this returns the poisson emission probability of seeing observation
		// for a given state and background
		// observation here must be an integer
		double intensity=state+background;
		if(observation>100){
			double err=Math.sqrt(intensity);
			return emit_p(observation,err,intensity);
		}else{
			return intpow(intensity,(int)observation)*Math.exp(-intensity)/factorial((int)observation);
		}
	}

	public double intpow(double val,int exponent){
		if(exponent==0){
			return 1.0;
		}
		double retval=val;
		for(int i=1;i<exponent;i++){
			retval*=val;
		}
		return retval;
	}

	public double factorial(int number){
		double factval;
		int i;
		if(number>0){
			factval=1.0;
			for(i=1;i<=number;i++){
				factval*=i;
			}
			return factval;
		}else{
			return 1.0;
		}
	}

	public double[][] bleach_prob_matrix(int nmolecules,double rate){
		double[][] mat=new double[nmolecules+1][nmolecules+1];
		for(int i=nmolecules;i>1;i--){
			double psingle=i*rate;
			double pdouble=rate*rate*0.5*(i*i-i);
			mat[nmolecules-i][nmolecules-i+1]=psingle;
			mat[nmolecules-i][nmolecules-i+2]=pdouble;
			mat[nmolecules-i][nmolecules-i]=1.0-psingle-pdouble;
		}
		mat[nmolecules][nmolecules]=1.0;
		mat[nmolecules-1][nmolecules-1]=1.0-rate;
		mat[nmolecules-1][nmolecules]=rate;
		return mat;
	}

}
