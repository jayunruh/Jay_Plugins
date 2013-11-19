/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class fit_corr{
	// here we fit a correlation function using hybrid grid search/linear least
	// squares algorithm
	public double mintd,maxtd,tdmult,r;
	public boolean third;
	public int psftype;

	public fit_corr(double mintd1,double maxtd1,double tdmult1,double r1){
		mintd=mintd1;
		maxtd=maxtd1;
		tdmult=tdmult1;
		r=r1;
		third=false;
		psftype=0;
	}

	public fit_corr(){
		mintd=0.001;
		maxtd=5.0;
		tdmult=1.05;
		r=5.0;
		third=false;
		psftype=0;
	}

	public fit_corr(double mintd1,double maxtd1,double tdmult1,double r1,boolean third1){
		mintd=mintd1;
		maxtd=maxtd1;
		tdmult=tdmult1;
		r=r1;
		third=third1;
		psftype=0;
	}

	public double[] fitac(float[] ac,float[] xvals,boolean skip,boolean base){
		double[] newxvals=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			newxvals[i]=(double)xvals[i];
		}
		return fitac(ac,newxvals,skip,base);
	}

	public double[] fitac(float[] ac,double[] xvals,boolean skip,boolean base){
		// here we fit the autocorrelation by taud gridsearch with linear least
		// squares for the g0 and baseline values
		double c2min=0.0;
		double[] minparams=new double[4];
		boolean first=true;
		for(double td=mintd;td<maxtd;td*=tdmult){
			double[] coef=fit_linear_ac(td,ac,xvals,skip,base);
			double[] fitparams={0.0,coef[0],td};
			if(base){
				fitparams[0]=coef[0];
				fitparams[1]=coef[1];
			}
			double c2val=c2(fitparams,ac,xvals,skip);
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

	public double[] fitac(float[] ac,double[] xvals,double td1,boolean fixtd,boolean skip,boolean base){
		if(!fixtd){
			return fitac(ac,xvals,skip,base);
		}else{
			double[] temp=fit_linear_ac(td1,ac,xvals,skip,base);
			if(base){
				double[] temp2={temp[0],temp[1],td1};
				double c2val=c2(temp2,ac,xvals,skip);
				double[] temp3={temp[0],temp[1],td1,c2val};
				return temp3;
			}else{
				double[] temp2={0.0,temp[0],td1};
				double c2val=c2(temp2,ac,xvals,skip);
				double[] temp3={0.0,temp[0],td1,c2val};
				return temp3;
			}
		}
	}

	public double[] fitac_2comp(float[] ac,float[] xvals,boolean skip,boolean base){
		double[] newxvals=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			newxvals[i]=(double)xvals[i];
		}
		return fitac_2comp(ac,newxvals,skip,base);
	}

	public double[] fitac_2comp(float[] ac,double[] xvals,boolean skip,boolean base){
		// here we fit the autocorrelation by taud gridsearch with linear least
		// squares for the g0 and baseline values
		double c2min=0.0;
		double[] minparams=new double[6];
		boolean first=true;
		for(double td1=mintd;td1<(maxtd/5.0);td1*=tdmult){
			for(double td2=mintd*5.0;td2<maxtd;td2*=tdmult){
				double[] coef=fit_linear_ac_2comp(td1,td2,ac,xvals,skip,base);
				double[] fitparams={0.0,coef[0],td1,coef[1],td2};
				if(base){
					fitparams[0]=coef[0];
					fitparams[1]=coef[1];
					fitparams[3]=coef[2];
				}
				double c2val=c2_2comp(fitparams,ac,xvals,skip);
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

	public double[] fitac_2comp(float[] ac,double[] xvals,double td11,boolean fixtd1,double td21,boolean fixtd2,boolean skip,boolean base){
		// here we fit the autocorrelation by taud gridsearch with linear least
		// squares for the g0 and baseline values
		if(!fixtd1&&!fixtd2){
			return fitac_2comp(ac,xvals,skip,base);
		}else{
			if(fixtd1&&fixtd2){
				double[] temp=fit_linear_ac_2comp(td11,td21,ac,xvals,skip,base);
				if(base){
					double[] temp2={temp[0],temp[1],td11,temp[2],td21};
					double c2val=c2_2comp(temp2,ac,xvals,skip);
					double[] temp3={temp[0],temp[1],td11,temp[2],td21,c2val};
					return temp3;
				}else{
					double[] temp2={0.0,temp[0],td11,temp[1],td21};
					double c2val=c2_2comp(temp2,ac,xvals,skip);
					double[] temp3={0.0,temp[0],td11,temp[1],td21,c2val};
					return temp3;
				}
			}else{
				double c2min=0.0;
				double[] minparams=new double[6];
				boolean first=true;
				double td1=td11;
				double tempmintd=td1*5.0;
				double tempmaxtd=maxtd;
				if(fixtd2){
					td1=td21;
					tempmintd=mintd;
					tempmaxtd=td1/5.0;
				}
				for(double td2=tempmintd;td2<tempmaxtd;td2*=tdmult){
					double[] coef=fit_linear_ac_2comp(td1,td2,ac,xvals,skip,base);
					double[] fitparams={0.0,coef[0],td1,coef[1],td2};
					if(base){
						fitparams[0]=coef[0];
						fitparams[1]=coef[1];
						fitparams[3]=coef[2];
					}
					double c2val=c2_2comp(fitparams,ac,xvals,skip);
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
				if(fixtd2){
					double temp=minparams[1];
					minparams[1]=minparams[3];
					minparams[3]=temp;
					minparams[2]=minparams[4];
					minparams[4]=td21;
				}
				return minparams;
			}
		}
	}

	public double[] fit_linear_ac(double td,float[] ac,float[] xvals,boolean skip,boolean base){
		// here we fit the autocorrelation with linear least squares for the g0
		// and baseline values
		double[] temp=new double[xvals.length];
		for(int i=0;i<xvals.length;i++)
			temp[i]=(double)xvals[i];
		return fit_linear_ac(td,ac,temp,skip,base);
	}

	public double[] fit_linear_ac(double td,float[] ac,double[] xvals,boolean skip,boolean base){
		// here we fit the autocorrelation with linear least squares for the g0
		// and baseline values
		int corrlength=ac.length;
		double sumx2=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double sumxy=0.0;
		double[] tempparams={0.0,1.0,td};
		int starti=0;
		if(skip){
			starti=1;
		}
		for(int i=starti;i<corrlength;i++){
			double corval=corfunc(tempparams,i,xvals);
			sumx2+=corval*corval;
			sumx+=corval;
			sumy+=(double)ac[i];
			sumxy+=(double)ac[i]*corval;
		}
		if(base){
			double divider=(double)(corrlength-starti)*sumx2-sumx*sumx;
			double baseline=(sumx2*sumy-sumx*sumxy)/divider;
			double g0=((double)(corrlength-starti)*sumxy-sumx*sumy)/divider;
			double[] fitparams={baseline,g0};
			return fitparams;
		}else{
			double g0=sumxy/sumx2;
			double[] fitparams={g0};
			return fitparams;
		}
	}

	public double[] fit_linear_ac_2comp(double td1,double td2,float[] ac,double[] xvals,boolean skip,boolean base){
		// here we fit the autocorrelation with linear least squares for the
		// g01, g02, and baseline values
		int corrlength=ac.length;
		double[] tempparams={0.0,1.0,td1};
		double[] tempparams2={0.0,1.0,td2};
		int starti=0;
		if(skip){
			starti=1;
		}
		double[][] cor1=new double[2][corrlength-starti];
		float[] data=new float[corrlength-starti];
		for(int i=starti;i<corrlength;i++){
			cor1[0][i-starti]=corfunc(tempparams,i,xvals);
			cor1[1][i-starti]=corfunc(tempparams2,i,xvals);
			data[i-starti]=ac[i];
		}
		linleastsquares lls=new linleastsquares(cor1,base,0,corrlength-starti-1);
		return lls.fitdata(data,null);
	}

	public double corfunc(double[] params,int indvar,double[] xvals){
		if(!third){
			if(psftype==0){
				return params[0]+params[1]/((1.0+xvals[indvar]/params[2])*Math.sqrt(1.0+xvals[indvar]/(params[2]*r*r)));
			}else{
				if(psftype==1){
					return params[0]+params[1]/(1.0+xvals[indvar]/params[2]);
				}else{
					return params[0]+params[1]/(Math.sqrt(1.0+xvals[indvar]/params[2])*Math.sqrt(1.0+xvals[indvar]/(params[2]*r*r)));
				}
			}
		}else{
			return params[0]+params[1]/((1.0+(4.0/3.0)*xvals[indvar]/params[2])*Math.sqrt(1.0+(4.0/3.0)*xvals[indvar]/(params[2]*r*r)));
		}
	}

	public double[] corfunc_array(double[] params,double[] xvals){
		double[] fit=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			fit[i]=corfunc(params,i,xvals);
		}
		return fit;
	}

	public float[] corfunc_arrayf(double[] params,double[] xvals){
		float[] fit=new float[xvals.length];
		for(int i=0;i<xvals.length;i++){
			fit[i]=(float)corfunc(params,i,xvals);
		}
		return fit;
	}

	public double c2(double[] params,float[] ac,float[] xvals,boolean skip){
		// here we fit the autocorrelation with linear least squares for the g0
		// and baseline values
		double[] temp=new double[xvals.length];
		for(int i=0;i<xvals.length;i++)
			temp[i]=(double)xvals[i];
		return c2(params,ac,temp,skip);
	}

	public double c2(double[] params,float[] data,double[] xvals,boolean skip){
		double tempc2=0.0;
		int startx=0;
		if(skip){
			startx=1;
		}
		for(int i=startx;i<data.length;i++){
			double resid=corfunc(params,i,xvals)-(double)data[i];
			tempc2+=resid*resid;
		}
		return tempc2;
	}

	public double c2_2comp(double[] params,float[] data,double[] xvals,boolean skip){
		double tempc2=0.0;
		int startx=0;
		if(skip){
			startx=1;
		}
		double[] params1={0.0,params[1],params[2]};
		double[] params2={0.0,params[3],params[4]};
		for(int i=startx;i<data.length;i++){
			double resid=corfunc(params1,i,xvals)+corfunc(params2,i,xvals)+params[0]-(double)data[i];
			tempc2+=resid*resid;
		}
		return tempc2;
	}

}
