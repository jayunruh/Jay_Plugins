/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.jstatistics;
import jalgs.matrixsolve;

import java.util.Arrays;

public class fit_gaussian{
	/*
	 * this class fits gaussian functions according to the methodof thompson et
	 * al (2002) Biophys J, vol 82, p 2775stdev is the PSF stdev (w0/2), fitsize
	 * is the size of thesquare image used for fitting, centertolerance is the
	 * centerresolution tolerance used as convergence criteriause this only if
	 * the baseline is well known
	 */
	public float stdev;
	public float toler;
	public int fitsize,fitsized2;
	//private float[] model_gaus=new float[10000];
	gausfunc gf;

	public fit_gaussian(float stdev1,int fitsize1,float toler1){
		stdev=stdev1;
		fitsize=fitsize1;
		fitsized2=(int)(fitsize/2.0f);
		toler=toler1;
		gf=new gausfunc();
		/*for(int i=0;i<10000;i++){
			double r=0.01*(double)i;
			model_gaus[i]=(float)Math.exp(-(r*r)/(2.0*stdev*stdev));
		}*/
	}
	
	/*************
	 * this is a constructor for those methods that don't require other public variables
	 */
	public fit_gaussian(){
		gf=new gausfunc();
		toler=0.0001f;
	}

	/*************
	 * This is the thompson et al version for 2D gaussians
	 * @param params
	 * @param fitstats
	 * @param data
	 * @param fitmask
	 * @param S
	 * @param baseline
	 * @return
	 */
	public float[] do_fit_gaussian(float[] params,float[] fitstats,float[] data,boolean[] fitmask,float S,float baseline){
		// params contains xc, yc, amplitude
		// fitstats returns iterations and chisquared
		// fit mask is true for masked pixels
		// S is the slope of var vs. intensity for the noise
		float newx=params[0];
		float oldx;
		float newy=params[1];
		float oldy;
		for(int i=0;i<200;i++){
			oldx=newx;
			oldy=newy;
			newx=0.0f;
			newy=0.0f;
			float temp1=0.0f;
			for(int j=0;j<fitsize;j++){
				float multiplier=(float)getinterpgaus(Math.abs(params[1]-j));
				for(int k=0;k<fitsize;k++){
					if(!fitmask[k+j*fitsize]){
						float temp=(data[k+j*fitsize]-baseline)*multiplier*(float)getinterpgaus(Math.abs(params[0]-k));
						temp1+=(temp*temp);
						newx+=temp*k;
						newy+=temp*j;
					}
				}
			}
			newx/=temp1;
			newy/=temp1;
			if(Math.abs(newx-oldx)<toler&&Math.abs(newy-oldy)<toler){
				fitstats[0]=i+1;
				break;
			}
		}
		float temp1=0.0f;
		float amplitude=0.0f;
		for(int j=0;j<fitsize;j++){
			float multiplier=(float)getinterpgaus(Math.abs(params[1]-j));
			for(int k=0;k<fitsize;k++){
				float temp=multiplier*(float)getinterpgaus(Math.abs(params[0]-k));
				temp1+=temp*temp;
				amplitude+=temp*(data[k+j*fitsize]-baseline);
			}
		}
		amplitude/=temp1;
		// amplitude*=2.0f*stdev*stdev*(float)Math.PI;
		params[0]=newx;
		params[1]=newy;
		params[2]=amplitude;
		float[] fit=new float[fitsize*fitsize];
		for(int j=0;j<fitsize;j++){
			float multiplier=(float)getinterpgaus(Math.abs(params[1]-j));
			for(int k=0;k<fitsize;k++){
				fit[k+j*fitsize]=amplitude*multiplier*(float)getinterpgaus(Math.abs(params[0]-k))+baseline;
			}
		}
		fitstats[1]=(float)calculate_c2_fit(fit,3,data,S,fitmask);
		return fit;
	}

	/*******************
	 * a 2D gaussian fit with a fixed stdev
	 * @param params
	 * @param fitstats
	 * @param data
	 * @param fitmask
	 * @param S
	 * @param maxiter
	 * @return
	 */
	public float[] do_fit_gaussian2(float[] params,float[] fitstats,float[] data,boolean[] fitmask,float S,int maxiter){
		// here we have a bit more complex version based on non-linear least
		// squares
		// here params is xc,yc,baseline,amp
		float[] fit=null;
		double oldc2,c2;
		c2=0.0f;
		int iter=0;
		if(maxiter==0){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)getinterpgaus(Math.abs(params[1]-i));
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)getinterpgaus(Math.abs(params[0]-j))+params[2];
					}
				}
			}
			c2=calculate_c2_fit(fit,4,data,S,fitmask);
		}
		double[][] jacobian=new double[4][4];
		double[] jvector=new double[4];
		iterloop: for(iter=0;iter<maxiter;iter++){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)getinterpgaus(Math.abs(params[1]-i));
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)getinterpgaus(Math.abs(params[0]-j))+params[2];
					}
				}
			}
			oldc2=c2;
			c2=calculate_c2_fit(fit,4,data,S,fitmask);
			if(iter>0){
				if(Math.abs((oldc2-c2)/c2)<toler){
					break iterloop;
				}
			}
			jacobian=new double[4][4];
			jvector=new double[4];
			for(int i=0;i<fitsize;i++){
				double tempy=((double)i-(double)params[1])/((double)stdev*(double)stdev);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						double temp=fit[j+i*fitsize];
						double temp2=data[j+i*fitsize];
						double tempx=((double)j-(double)params[0])/((double)stdev*(double)stdev);
						double[] d={(temp-params[2])*tempx,(temp-params[2])*tempy,1.0,(temp-params[2])/params[3]};
						for(int k=0;k<4;k++){
							for(int l=0;l<=k;l++){
								// jacobian[k][l]+=(double)((d[k]*d[l])/(temp2*S));
								jacobian[k][l]+=d[k]*d[l];
							}
							// jvector[k]+=(double)((d[k]*(temp2-temp))/(temp2*S));
							jvector[k]+=d[k]*(temp2-temp);
						}
					}
				}
			}
			for(int k=0;k<3;k++){
				for(int l=(k+1);l<4;l++){
					jacobian[k][l]=jacobian[l][k];
				}
			}
			double[] dparams=new double[4];
			(new matrixsolve()).gjsolve(jacobian,jvector,dparams,4);
			params[0]+=dparams[0];
			params[1]+=dparams[1];
			params[2]+=dparams[2];
			params[3]+=dparams[3];
			if(params[0]<(0.25f*fitsize)){
				params[0]=0.25f*fitsize;
			}
			if(params[0]>(0.75f*fitsize)){
				params[0]=0.75f*fitsize;
			}
			if(params[1]<(0.25f*fitsize)){
				params[1]=0.25f*fitsize;
			}
			if(params[1]>(0.75f*fitsize)){
				params[1]=0.75f*fitsize;
			}
			if(params[2]<-100.0f){
				params[2]=-100.0f;
			}
			if(params[3]<1.0f){
				params[3]=1.0f;
			}
		}
		fitstats[0]=iter;
		fitstats[1]=(float)c2;
		/*
		 * fit=new float[4*4]; for(int i=0;i<4;i++){ for(int j=0;j<4;j++){
		 * fit[j+i*4]=(float)jacobian[i][j]; } }
		 */
		return fit;
	}
	
	/*********
	 * a 2d gaussian fit where the stdev is fit as well
	 * @param params
	 * @param fitstats
	 * @param data
	 * @param fitmask
	 * @param S
	 * @param maxiter
	 * @return
	 */
	public float[] do_fit_gaussian3(float[] params,float[] fitstats,float[] data,boolean[] fitmask,float S,int maxiter){
		//this nlls version fits the stdev
		// here params is xc,yc,baseline,amp,stdev
		float[] fit=null;
		double oldc2,c2;
		c2=0.0f;
		int iter=0;
		if(maxiter==0){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]-i),params[4]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)gf.getinterpgaus(Math.abs(params[0]-j),params[4])+params[2];
					}
				}
			}
			c2=calculate_c2_fit(fit,5,data,S,fitmask);
		}
		double[][] jacobian=new double[5][5];
		double[] jvector=new double[5];
		iterloop: for(iter=0;iter<maxiter;iter++){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]-i),params[4]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)gf.getinterpgaus(Math.abs(params[0]-j),params[4])+params[2];
					}
				}
			}
			oldc2=c2;
			c2=calculate_c2_fit(fit,5,data,S,fitmask);
			if(iter>0){
				if(Math.abs((oldc2-c2)/c2)<toler){
					break iterloop;
				}
			}
			jacobian=new double[5][5];
			jvector=new double[5];
			for(int i=0;i<fitsize;i++){
				double tempy=((double)i-(double)params[1])/((double)params[4]*(double)params[4]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						double temp=fit[j+i*fitsize];
						double temp2=data[j+i*fitsize];
						double tempx=((double)j-(double)params[0])/((double)params[4]*(double)params[4]);
						double[] d={(temp-params[2])*tempx,(temp-params[2])*tempy,1.0,(temp-params[2])/params[3],(temp-params[2])*(tempx*tempx+tempy*tempy)*params[4]};
						for(int k=0;k<5;k++){
							for(int l=0;l<=k;l++){
								// jacobian[k][l]+=(double)((d[k]*d[l])/(temp2*S));
								jacobian[k][l]+=d[k]*d[l];
							}
							// jvector[k]+=(double)((d[k]*(temp2-temp))/(temp2*S));
							jvector[k]+=d[k]*(temp2-temp);
						}
					}
				}
			}
			for(int k=0;k<4;k++){
				for(int l=(k+1);l<5;l++){
					jacobian[k][l]=jacobian[l][k];
				}
			}
			double[] dparams=new double[5];
			(new matrixsolve()).gjsolve(jacobian,jvector,dparams,5);
			//now we implement the updates and constraints
			params[0]+=dparams[0];
			params[1]+=dparams[1];
			params[2]+=dparams[2];
			params[3]+=dparams[3];
			params[4]+=dparams[4];
			if(params[0]<(0.25f*fitsize)) params[0]=0.25f*fitsize;
			if(params[0]>(0.75f*fitsize)) params[0]=0.75f*fitsize;
			if(params[1]<(0.25f*fitsize)) params[1]=0.25f*fitsize;
			if(params[1]>(0.75f*fitsize)) params[1]=0.75f*fitsize;
			if(params[2]<-100.0f) params[2]=-100.0f;
			if(params[3]<1.0f) params[3]=1.0f;
			if(params[4]<0.9f*stdev) params[4]=0.9f*stdev;
			if(params[4]>4.0f*stdev) params[4]=4.0f*stdev;
		}
		fitstats[0]=iter;
		fitstats[1]=(float)c2;
		/*
		 * fit=new float[4*4]; for(int i=0;i<4;i++){ for(int j=0;j<4;j++){
		 * fit[j+i*4]=(float)jacobian[i][j]; } }
		 */
		return fit;
	}
	
	/**********
	 * this does a simple linear least squares with provided fits and a baseline
	 * @param data
	 * @param fit1
	 * @param fit2
	 * @param fitmask
	 * @return
	 */
	public float[] get_2fit_amps(float[] data,float[] fit1,float[] fit2,boolean[] fitmask){
		double[][] jacobian=new double[3][3];
		double[] jvector=new double[3];
		float[] ones=new float[fit1.length];
		Arrays.fill(ones,1.0f);
		float[][] indvars={ones,fit1,fit2};
		for(int i=0;i<3;i++){
			for(int j=0;j<=i;j++){
				for(int k=0;k<fit1.length;k++){
					if(!fitmask[k]) jacobian[i][j]+=(double)indvars[i][k]*(double)indvars[j][k];
				}
				if(j!=i){
					jacobian[j][i]=jacobian[i][j];
				}
			}
			for(int k=0;k<fit1.length;k++){
				if(!fitmask[k]) jvector[i]+=(double)data[k]*(double)indvars[i][k];
			}
		}
		double[] outcoef=new double[3];
		(new matrixsolve()).gjsolve(jacobian,jvector,outcoef,3);
		float[] outcoef2={(float)outcoef[0],(float)outcoef[1],(float)outcoef[2]};
		return outcoef2;
	}

	public float[] do_fit_gaussian_assym(float[] params,float[] fitstats,float[] data,boolean[] fitmask,float S,int maxiter){
		// here we have a bit more complex version based on non-linear least
		// squares
		// this version includes assymetry for 3D PALM analysis
		// here params is xc,yc,baseline,amp,xstdev,ystdev
		float[] fit=null;
		double oldc2,c2;
		c2=0.0f;
		int iter=0;
		if(params[4]<stdev){
			params[4]=stdev;
		}
		if(params[5]<stdev){
			params[5]=stdev;
		}
		if(params[4]<params[5]){
			params[4]=stdev;
		}else{
			params[5]=stdev;
		}
		double xstdevmult=params[4]/(double)stdev;
		double ystdevmult=params[5]/(double)stdev;
		if(maxiter==0){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)getinterpgaus(ystdevmult*Math.abs(params[1]-i));
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)getinterpgaus(xstdevmult*Math.abs(params[0]-j))+params[2];
					}
				}
			}
			c2=calculate_c2_fit(fit,6,data,S,fitmask);
		}
		double[][] jacobian=null;
		double[] jvector=null;
		iterloop: for(iter=0;iter<maxiter;iter++){
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)getinterpgaus(ystdevmult*Math.abs(params[1]-i));
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit[j+i*fitsize]=params[3]*multiplier*(float)getinterpgaus(xstdevmult*Math.abs(params[0]-j))+params[2];
					}
				}
			}
			oldc2=c2;
			c2=calculate_c2_fit(fit,6,data,S,fitmask);
			if(iter>0){
				if(Math.abs((oldc2-c2)/c2)<toler){
					break iterloop;
				}
			}
			jacobian=new double[6][6];
			jvector=new double[6];
			for(int i=0;i<fitsize;i++){
				double tempy=((double)i-(double)params[1])/((double)params[5]*(double)params[5]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						double temp=fit[j+i*fitsize];
						double temp2=data[j+i*fitsize];
						double tempx=((double)j-(double)params[0])/((double)params[4]*(double)params[4]);
						double[] d={(temp-params[2])*tempx,(temp-params[2])*tempy,1.0,(temp-params[2])/params[3],(temp-params[2])*tempx*tempx*params[4],(temp-params[2])*tempy*tempy*params[5]};
						for(int k=0;k<6;k++){
							for(int l=0;l<=k;l++){
								// jacobian[k][l]+=(double)((d[k]*d[l])/(temp2*S));
								jacobian[k][l]+=d[k]*d[l];
							}
							// jvector[k]+=(double)((d[k]*(temp2-temp))/(temp2*S));
							jvector[k]+=d[k]*(temp2-temp);
						}
					}
				}
			}
			for(int k=0;k<5;k++){
				for(int l=(k+1);l<6;l++){
					jacobian[k][l]=jacobian[l][k];
				}
			}
			double[] dparams=new double[6];
			(new matrixsolve()).gjsolve(jacobian,jvector,dparams,6);
			params[0]+=dparams[0];
			params[1]+=dparams[1];
			params[2]+=dparams[2];
			params[3]+=dparams[3];
			params[4]+=dparams[4];
			params[5]+=dparams[5];
			if(params[0]<(0.25f*fitsize)){
				params[0]=0.25f*fitsize;
			}
			if(params[0]>(0.75f*fitsize)){
				params[0]=0.75f*fitsize;
			}
			if(params[1]<(0.25f*fitsize)){
				params[1]=0.25f*fitsize;
			}
			if(params[1]>(0.75f*fitsize)){
				params[1]=0.75f*fitsize;
			}
			if(params[2]<-100.0f){
				params[2]=-100.0f;
			}
			if(params[3]<1.0f){
				params[3]=1.0f;
			}
			if(params[4]<stdev){
				params[4]=stdev;
			}
			if(params[4]>(0.5f*fitsize)){
				params[4]=(0.5f*fitsize);
			}
			if(params[5]<stdev){
				params[5]=stdev;
			}
			if(params[5]>(0.5f*fitsize)){
				params[5]=(0.5f*fitsize);
			}
			if(params[5]<params[4]){
				params[5]=stdev;
			}else{
				params[4]=stdev;
			}
		}
		fitstats[0]=iter;
		fitstats[1]=(float)c2;
		/*
		 * fit=new float[4*4]; for(int i=0;i<4;i++){ for(int j=0;j<4;j++){
		 * fit[j+i*4]=(float)jacobian[i][j]; } }
		 */
		return fit;
	}

	private double getinterpgaus(double r){
		/*int rp=(int)(r*100.0);
		double rem=r*100.0-(double)rp;
		if(rp<9999){
			return rem*(model_gaus[rp+1]-model_gaus[rp])+model_gaus[rp];
		}else{
			return 0.0;
		}*/
		return gf.getinterpgaus(r,1.0);
	}

	public double calculate_c2_fit(float[] fit,int numfit,float[] data,float S,boolean[] fitmask){
		// this function calculates chisquared using a fit array
		// var=S*Intensity
		int length=data.length;
		double tempc2=0.0;
		int nonfit=0;
		for(int i=0;i<length;i++){
			if(!fitmask[i]){
				if(data[i]>0.0f){
					tempc2+=(double)(fit[i]-data[i])*(fit[i]-data[i])*(1.0f/(data[i]*S));
				}else{
					tempc2+=(double)(fit[i]-data[i])*(fit[i]-data[i]);
				}
			}else{
				nonfit++;
			}
		}
		return tempc2/(length-nonfit-numfit);
	}
	
	/************
	 * this guesses parameters for a 1D data set, xvals1 can be null, returns base, xc,stdev,amp
	 * @param xvals
	 * @param yvals
	 * @param length
	 * @return
	 */
	public static double[] guess1DParams(float[] xvals1,float[] yvals,int length){
		float[] xvals=new float[length];
		if(xvals1!=null) xvals=xvals1;
		else for(int i=0;i<length;i++) xvals[i]=(float)i;
		float min=yvals[0];
		int maxloc=0;
		float max=yvals[0];
		for(int i=1;i<length;i++){
			if(yvals[i]<min){min=yvals[i];}
			if(yvals[i]>max){max=yvals[i]; maxloc=i;}
		}
		float halfmax=(max+min)/2.0f;
		float fwhm=0.0f;
		for(int i=(maxloc+1);i<length;i++){
			if(yvals[i]<halfmax){
				fwhm+=0.5f*(xvals[i]-xvals[maxloc]);
				break;
			}
			if(i==(length-1)){
				fwhm+=0.5f*(xvals[i]-xvals[maxloc]);
			}
		}
		for(int i=(maxloc-1);i>=0;i--){
			if(yvals[i]<halfmax){
				fwhm+=0.5f*(xvals[maxloc]-xvals[i]);
				break;
			}
			if(i==0){
				fwhm+=0.5f*(xvals[maxloc]-xvals[i]);
			}
		}

		//parameters are baseline,xc1,stdev1,amp1
		double[] params={(double)min,(double)xvals[maxloc],(double)fwhm/2.35,(double)(max-min)};
		return params;
	}
	
	/*****************
	 * generates constraints for the 1D data set, xvals1 can be null
	 * @param xvals1
	 * @param params: base,xc,stdev,amp
	 * @return
	 */
	public static double[][] get1DConstraints(float[] xvals1,double[] params,int length){
		double[][] constraints=new double[2][4];
		float[] xvals=new float[length];
		if(xvals1!=null) xvals=xvals1;
		else for(int i=0;i<length;i++) xvals[i]=(float)i;
		constraints[0][0]=params[0]-0.5*params[3]; constraints[1][0]=params[3]+params[0];
		constraints[0][1]=(double)xvals[0]; constraints[1][1]=(double)xvals[xvals.length-1];
		constraints[0][2]=0.2*(double)(xvals[1]-xvals[0]); constraints[1][2]=(double)(xvals[xvals.length-1]-xvals[0]);
		constraints[0][3]=0.0; constraints[1][3]=10.0*params[3];
		return constraints;
	}
	
	/****************
	 * a 1D gaussian fitting function with no use intervention
	 * @param xvals1: can be null
	 * @param yvals
	 * @param params: base, xc, stdev,amp
	 * @param stats
	 * @param constraints
	 * @param fixes
	 * @return
	 */
	public float[] run1DFit(float[] xvals1,float[] yvals,double[] params1,double[] stats,double[][] constraints1,int[] fixes){
		float[] tempx=new float[yvals.length];
		if(xvals1!=null) tempx=xvals1;
		else for(int i=0;i<tempx.length;i++) tempx[i]=(float)i;
		double[] params=params1;
		if(params1==null) params=guess1DParams(tempx,yvals,yvals.length); //note that if this happens, we will not return the params variables
		double sumparams=0.0; for(int i=0;i<params.length;i++) sumparams+=params[i]*params[i];
		if(sumparams==0.0) params=guess1DParams(tempx,yvals,yvals.length);
		double[][] constraints=constraints1;
		if(constraints1==null) constraints=get1DConstraints(tempx,params,yvals.length);
		int maxiter=10;
		float[] fit=new float[yvals.length];
		for(int i=0;i<yvals.length;i++) fit[i]=(float)gf.getinterpgaus(Math.abs(tempx[i]-params[1]),params[2])*(float)params[3]+(float)params[0];
		double[][] jacobian=new double[4][4];
		double[] jvector=new double[4];
		int iter=0;
		double c2=calculate_c2_1D_fit(fit,4,yvals,0.0f);
		double oldc2;
		iterloop: for(iter=0;iter<maxiter;iter++){
			fit=new float[yvals.length];
			for(int i=0;i<yvals.length;i++) fit[i]=(float)gf.getinterpgaus(Math.abs(tempx[i]-params[1]),params[2])*(float)params[3]+(float)params[0];
			oldc2=c2;
			c2=calculate_c2_1D_fit(fit,4,yvals,0.0f);
			if(iter>0){
				if(Math.abs((oldc2-c2)/c2)<toler){
					break iterloop;
				}
			}
			jacobian=new double[4][4];
			jvector=new double[4];
			for(int i=0;i<yvals.length;i++){
				double temp=fit[i];
				double temp2=yvals[i];
				double tempx2=((double)tempx[i]-(double)params[1])/(params[2]*params[2]);
				double[] d={1.0,(temp-params[0])*tempx2,(temp-params[0])*(tempx2*tempx2)*params[2],(temp-params[0])/params[3]};
				for(int k=0;k<4;k++){
					for(int l=0;l<=k;l++){
						jacobian[k][l]+=d[k]*d[l];
					}
					jvector[k]+=d[k]*(temp2-temp);
				}
			}
			for(int k=0;k<3;k++){
				for(int l=(k+1);l<4;l++){
					jacobian[k][l]=jacobian[l][k];
				}
			}
			double[] dparams=new double[4];
			(new matrixsolve()).gjsolve(jacobian,jvector,dparams,4);
			if(fixes[0]==0) params[0]+=dparams[0];
			if(fixes[1]==0) params[1]+=dparams[1];
			if(fixes[2]==0) params[2]+=dparams[2];
			if(fixes[3]==0) params[3]+=dparams[3];
			for(int i=0;i<4;i++){
				if(fixes[i]==0){
					if(Double.isNaN(params[i])) params[i]=constraints[0][i];
					if(Double.isInfinite(params[i])) params[i]=constraints[1][i];
					if(params[i]<constraints[0][i]) params[i]=constraints[0][i];
					if(params[i]>constraints[1][i]) params[i]=constraints[1][i];
				}
			}
		}
		stats[0]=(double)iter;
		stats[1]=c2;
		for(int i=0;i<params.length;i++) params1[i]=params[i];
		return fit;
	}
	
	public static double calculate_c2_1D_fit(float[] fit,int numfit,float[] data,float S){
		// this function calculates chisquared using a fit array
		// var=S*Intensity
		int length=data.length;
		double tempc2=0.0;
		for(int i=0;i<length;i++){
			if(data[i]>0.0f && S>0.0f){
				tempc2+=(double)(fit[i]-data[i])*(fit[i]-data[i])*(1.0f/(data[i]*S));
			}else{
				tempc2+=(double)(fit[i]-data[i])*(fit[i]-data[i]);
			}
		}
		return tempc2/(length-numfit);
	}

}
