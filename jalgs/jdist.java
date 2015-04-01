/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class jdist{
	// here we have utility methods for calculating limiting statistical values
	// most of the routines are adapted from numerical recipes

	public double ibeta(double a,double b,double x){
		// Returns the incomplete beta function Ix(a, b). Adapted from numerical
		// recipes.
		double beta=0.0;
		if(x<0.0||x>1.0){
			return Double.NaN;
		}
		if(x!=0.0||x!=1.0){
			beta=Math.exp(sf.gammaln(a+b)-sf.gammaln(a)-sf.gammaln(b)+a*Math.log(x)+b*Math.log(1.0-x));
		}
		if(x<(a+1.0)/(a+b+2.0)){
			return beta*betacf(a,b,x)/a;
		}else{
			return 1.0-beta*betacf(b,a,1.0-x)/b;
		}
	}

	private double betacf(double a,double b,double x){
		// adapted from numerical recipes
		int MAXIT=100;
		double EPS=3.0e-7;
		double FPMIN=1.0e-30;
		double qab=a+b;
		double qap=a+1.0;
		double qam=a-1.0;
		double c=1.0;
		double d=1.0-qab*x/qap;
		if(Math.abs(d)<FPMIN){
			d=FPMIN;
		}
		d=1.0/d;
		double h=d;
		int m=0;
		while(m<MAXIT){
			int m2=2*(m+1);
			double aa=(m+1)*(b-m-1)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if(Math.abs(d)<FPMIN){
				d=FPMIN;
			}
			c=1.0+aa/c;
			if(Math.abs(c)<FPMIN){
				c=FPMIN;
			}
			d=1.0/d;
			h*=d*c;
			aa=-(a+m+1)*(qab+m+1)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if(Math.abs(d)<FPMIN)
				d=FPMIN;
			c=1.0+aa/c;
			if(Math.abs(c)<FPMIN)
				c=FPMIN;
			d=1.0/d;
			double del=d*c;
			h*=del;
			if(Math.abs(del-1.0)<EPS){
				return h;
			}
			m++;
		}
		return Double.NaN;
	}

	public double FCumProb(double numdof,double dendof,double F){
		double xval=numdof*F/(numdof*F+dendof);
		return ibeta(0.5*numdof,0.5*dendof,xval);
	}

	public double FLimit(double numdof,double dendof,double prob){
		double df=0.001;
		double F=1.0-df;
		double cumprob=0.0;
		double lastcumprob=0.0;
		do{
			F+=df;
			lastcumprob=cumprob;
			cumprob=FCumProb(numdof,dendof,F);
		}while(cumprob<prob);
		double fraction=(prob-lastcumprob)/(cumprob-lastcumprob);
		return fraction*(F-df)+(1.0-fraction)*F;
	}

	public double tCumProbTwoSampEqVar(double dof1,double dof2,double mean1,double mean2,double sterr1,double sterr2,double nulldiff,boolean twotailed){
		double n1=dof1+1.0;
		double n2=dof2+1.0;
		// assume that variances are equal and that dof=n-1;
		double S2pooled=(dof1*n1*sterr1*sterr1+dof2*n2*sterr2*sterr2)/(dof1+dof2);
		double sterrdiff=Math.sqrt(S2pooled*(1.0/n1+1.0/n2));
		double t=(Math.abs(mean1-mean2)-nulldiff)/sterrdiff;
		return tCumProb(dof1+dof2,t,twotailed);
	}

	public double tCumProbTwoSampNVar(double dof1,double dof2,double mean1,double mean2,double sterr1,double sterr2,double nulldiff,boolean twotailed){
		double n1=dof1+1.0;
		double n2=dof2+1.0;
		// assume that variances are not equal and that dof=n-1;
		double sterrdiff=Math.sqrt(sterr1*sterr1+sterr2*sterr2);
		double C=sterr1*sterr1/(sterr1*sterr1+sterr2*sterr2);
		double dof=dof1*dof2/(dof2*C*C+dof1*(1.0-C)*(1.0-C));
		double t=(Math.abs(mean1-mean2)-nulldiff)/sterrdiff;
		return tCumProb(dof,t,twotailed);
	}

	public double tCumProbOneSamp(double dof,double mean,double value,double sterr,boolean twotailed){
		// t is the number of standard errors the value is away from the mean
		// the standard error is given by stdev/sqrt(n);
		double t=(value-mean)/sterr;
		return tCumProb(dof,t,twotailed);
	}
	
	public double tLim(double dof,double prob,boolean twotailed){
		//this gets the critical t value for a specific probability
		double dt=0.001;
		double t=0.0-dt;
		double cumprob=0.0;
		double lastcumprob=0.0;
		do{
			t+=dt;
			lastcumprob=cumprob;
			cumprob=tCumProb(dof,t,twotailed);
		}while(cumprob>prob);
		double fraction=(prob-cumprob)/(lastcumprob-cumprob);
		return t-fraction*dt;
	}

	public double tCumProb(double dof,double t,boolean twotailed){
		//this gets the probabity for a specific t value
		if(!twotailed){
			return 0.5*ibeta(0.5*dof,0.5,dof/(dof+t*t));
		}else{
			return ibeta(0.5*dof,0.5,dof/(dof+t*t));
		}
	}

}
