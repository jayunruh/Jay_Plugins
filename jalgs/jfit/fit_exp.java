/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class fit_exp{
	public float tmult=1.05f;
	public float mint,maxt;
	public float minratio=2.0f;
	public boolean fix1,fix2;

	public fit_exp(float mint,float maxt){
		this.mint=mint;
		this.maxt=maxt;
	}

	public double[] fitdata(float[] data){
		float[] xvals=new float[data.length];
		for(int i=0;i<data.length;i++){
			xvals[i]=i;
		}
		return fitdata(data,xvals);
	}

	public double[] fitdata(float[] data,float[] xvals){
		// here we fit the decay by tau gridsearch with linear least squares for
		// the tau, amp, and baseline values
		double c2min=0.0;
		double[] minparams=new double[4];
		boolean first=true;
		float maxt1=maxt;
		if(mint==maxt || fix1) maxt1=mint*tmult;
		for(double t=mint;t<maxt1;t*=tmult){
			double[] coef=fit_linear(t,data,xvals);
			double[] fitparams={coef[0],coef[1],t};
			double c2val=c2(fitparams,data,xvals);
			if(first){
				c2min=c2val;
				minparams[0]=fitparams[0];
				minparams[1]=fitparams[1];
				minparams[2]=fitparams[2];
				minparams[3]=c2min;
				first=false;
			}else{
				if(c2val<c2min){
					c2min=c2val;
					minparams[0]=fitparams[0];
					minparams[1]=fitparams[1];
					minparams[2]=fitparams[2];
					minparams[3]=c2min;
				}
			}
		}
		return minparams;
	}

	public double[] fitdata2exp(float[] ac){
		float[] xvals=new float[ac.length];
		for(int i=0;i<ac.length;i++){
			xvals[i]=i;
		}
		return fitdata2exp(ac,xvals);
	}

	public double[] fitdata2exp(float[] ac,float[] xvals){
		// here we fit the decay by taud gridsearch with linear least squares
		// for the amp and baseline values
		double c2min=0.0;
		double[] minparams=new double[6];
		boolean first=true;
		float maxt1=maxt/minratio;
		if(fix1) maxt1=mint*tmult;
		for(double t1=mint;t1<maxt1;t1*=tmult){
			double mint2=t1*minratio;
			double maxt2=maxt;
			if(fix2){
				mint2=maxt;
				maxt2=maxt*tmult;
			}
			for(double t2=mint2;t2<maxt2;t2*=tmult){
				double[] coef=fit_linear2exp(t1,t2,ac,xvals);
				double[] fitparams={0.0,coef[0],t1,coef[1],t2};
				fitparams[0]=coef[0];
				fitparams[1]=coef[1];
				fitparams[3]=coef[2];
				double c2val=c2_2comp(fitparams,ac,xvals);
				if(first){
					c2min=c2val;
					System.arraycopy(fitparams,0,minparams,0,5);
					minparams[5]=c2min;
					first=false;
				}else{
					if(c2val<c2min){
						c2min=c2val;
						System.arraycopy(fitparams,0,minparams,0,5);
						minparams[5]=c2min;
					}
				}
			}
		}
		return minparams;
	}

	public double[] fitdata3exp(float[] ac){
		float[] xvals=new float[ac.length];
		for(int i=0;i<ac.length;i++){
			xvals[i]=i;
		}
		return fitdata3exp(ac,xvals);
	}

	public double[] fitdata3exp(float[] ac,float[] xvals){
		// here we fit the decay by taud gridsearch with linear least squares
		// for the amp and baseline values
		double c2min=0.0;
		double[] minparams=new double[8];
		boolean first=true;
		for(double t1=mint;t1<((maxt/minratio)/minratio);t1*=tmult){
			for(double t2=t1*minratio;t2<(maxt/minratio);t2*=tmult){
				for(double t3=t2*minratio;t3<maxt;t3*=tmult){
					double[] coef=fit_linear3exp(t1,t2,t3,ac,xvals);
					double[] fitparams={coef[0],coef[1],t1,coef[2],t2,coef[3],t3};
					double c2val=c2_3comp(fitparams,ac,xvals);
					if(first){
						c2min=c2val;
						System.arraycopy(fitparams,0,minparams,0,7);
						minparams[7]=c2min;
						first=false;
					}else{
						if(c2val<c2min){
							c2min=c2val;
							System.arraycopy(fitparams,0,minparams,0,7);
							minparams[7]=c2min;
						}
					}
				}
			}
		}
		return minparams;
	}

	public double[] fit_linear(double t,float[] data,float[] xvals){
		// here we fit the decay with linear least squares for the amp and
		// baseline values
		int length=data.length;
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		double[] tempparams={0.0,1.0,t};
		for(int i=0;i<length;i++){
			double fval=corfunc(tempparams,i,xvals);
			sumx2+=fval*fval;
			sumx+=fval;
			sumy+=data[i];
			sumxy+=data[i]*fval;
		}
		double divider=length*sumx2-sumx*sumx;
		double baseline=(sumx2*sumy-sumx*sumxy)/divider;
		double amp=(length*sumxy-sumx*sumy)/divider;
		double[] fitparams={baseline,amp};
		return fitparams;
	}

	public double[] fit_linear2exp(double t1,double t2,float[] data,float[] xvals){
		// here we fit the decay with linear least squares for the amp and
		// baseline values
		int corrlength=data.length;
		double[] tempparams={0.0,1.0,t1};
		double[] tempparams2={0.0,1.0,t2};
		double[][] cor1=new double[2][corrlength];
		for(int i=0;i<corrlength;i++){
			cor1[0][i]=corfunc(tempparams,i,xvals);
			cor1[1][i]=corfunc(tempparams2,i,xvals);
		}
		linleastsquares lls=new linleastsquares(cor1,true,0,corrlength-1);
		return lls.fitdata(data,null);
	}

	public double[] fit_linear3exp(double t1,double t2,double t3,float[] data,float[] xvals){
		// here we fit the decay with linear least squares for the amp and
		// baseline values
		int corrlength=data.length;
		double[] tempparams={0.0,1.0,t1};
		double[] tempparams2={0.0,1.0,t2};
		double[] tempparams3={0.0,1.0,t3};
		double[][] cor1=new double[3][corrlength];
		for(int i=0;i<corrlength;i++){
			cor1[0][i]=corfunc(tempparams,i,xvals);
			cor1[1][i]=corfunc(tempparams2,i,xvals);
			cor1[2][i]=corfunc(tempparams3,i,xvals);
		}
		linleastsquares lls=new linleastsquares(cor1,true,0,corrlength-1);
		return lls.fitdata(data,null);
	}

	public double corfunc(double[] params,int indvar,float[] xvals){
		return params[0]+params[1]*Math.exp(-xvals[indvar]/params[2]);
	}

	public double[] corfunc_array(double[] params,float[] xvals){
		double[] fit=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			fit[i]=corfunc(params,i,xvals);
		}
		return fit;
	}

	public float[] corfunc_arrayf(double[] params,float[] xvals){
		float[] fit=new float[xvals.length];
		for(int i=0;i<xvals.length;i++){
			fit[i]=(float)corfunc(params,i,xvals);
		}
		return fit;
	}

	public double c2(double[] params,float[] data,float[] xvals){
		double tempc2=0.0;
		for(int i=0;i<data.length;i++){
			double resid=corfunc(params,i,xvals)-data[i];
			tempc2+=resid*resid;
		}
		return tempc2;
	}

	public double c2_2comp(double[] params,float[] data,float[] xvals){
		double tempc2=0.0;
		double[] params1={0.0,params[1],params[2]};
		double[] params2={0.0,params[3],params[4]};
		for(int i=0;i<data.length;i++){
			double resid=corfunc(params1,i,xvals)+corfunc(params2,i,xvals)+params[0]-data[i];
			tempc2+=resid*resid;
		}
		return tempc2;
	}

	public double c2_3comp(double[] params,float[] data,float[] xvals){
		double tempc2=0.0;
		double[] params1={0.0,params[1],params[2]};
		double[] params2={0.0,params[3],params[4]};
		double[] params3={0.0,params[5],params[6]};
		for(int i=0;i<data.length;i++){
			double resid=corfunc(params1,i,xvals)+corfunc(params2,i,xvals)+corfunc(params3,i,xvals)+params[0]-data[i];
			tempc2+=resid*resid;
		}
		return tempc2;
	}

}
