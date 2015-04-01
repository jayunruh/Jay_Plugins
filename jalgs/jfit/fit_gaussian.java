/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

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
	
	/*******************************
	 * this version fits a pair of gaussians
	 * @param params
	 * @param fitstats
	 * @param data
	 * @param fitmask
	 * @param S
	 * @param maxiter
	 * @return
	 */
	public float[] do_fit_2gaussian(float[] params,float[] fitstats,float[] data,boolean[] fitmask,float S,int maxiter){
		//this nlls version fits the stdev
		//stdev was 4, amp was 3
		// here params is xc,yc,stdev,angle,dist
		// amp1, amp2, and baseline are determined from linear least squares
		float[] fit=null;
		double oldc2,c2;
		c2=0.0f;
		int iter=0;
		if(maxiter==0){
			float xd=0.5f*params[4]*(float)Math.cos(params[3]);
			float yd=0.5f*params[4]*(float)Math.sin(params[3]);
			float[] fit1=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]+yd-i),params[2]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit1[j+i*fitsize]=multiplier*(float)gf.getinterpgaus(Math.abs(params[0]+xd-j),params[2]);
					}
				}
			}
			float[] fit2=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]-yd-i),params[2]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit2[j+i*fitsize]=multiplier*(float)gf.getinterpgaus(Math.abs(params[0]-xd-j),params[2]);
					}
				}
			}
			float[] amps=get_2fit_amps(data,fit1,fit2,fitmask);
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize*fitsize;i++) fit[i]=amps[0]+amps[1]*fit1[i]+amps[2]*fit2[i];
			c2=calculate_c2_fit(fit,8,data,S,fitmask);
		}
		double[][] jacobian=new double[5][5];
		double[] jvector=new double[5];
		iterloop: for(iter=0;iter<maxiter;iter++){
			float xd=0.5f*params[4]*(float)Math.cos(params[3]);
			float yd=0.5f*params[4]*(float)Math.sin(params[3]);
			float[] fit1=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]+yd-i),params[2]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit1[j+i*fitsize]=multiplier*(float)gf.getinterpgaus(Math.abs(params[0]+xd-j),params[2]);
					}
				}
			}
			float[] fit2=new float[fitsize*fitsize];
			for(int i=0;i<fitsize;i++){
				float multiplier=(float)gf.getinterpgaus(Math.abs(params[1]-yd-i),params[2]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						fit2[j+i*fitsize]=multiplier*(float)gf.getinterpgaus(Math.abs(params[0]-xd-j),params[2]);
					}
				}
			}
			float[] amps=get_2fit_amps(data,fit1,fit2,fitmask);
			fit=new float[fitsize*fitsize];
			for(int i=0;i<fitsize*fitsize;i++) fit[i]=amps[0]+amps[1]*fit1[i]+amps[2]*fit2[i];

			oldc2=c2;
			c2=calculate_c2_fit(fit,8,data,S,fitmask);
			if(iter>0){
				if(Math.abs((oldc2-c2)/c2)<toler){
					break iterloop;
				}
			}
			jacobian=new double[5][5];
			jvector=new double[5];
			for(int i=0;i<fitsize;i++){
				double tempy1=((double)i-(double)params[1]-yd)/((double)params[2]*(double)params[2]);
				double tempy2=((double)i-(double)params[1]+yd)/((double)params[2]*(double)params[2]);
				for(int j=0;j<fitsize;j++){
					if(!fitmask[j+i*fitsize]){
						// here params was xc,yc,baseline,amp,stdev
						// here params is xc,yc,stdev,angle,dist
						int index=j+i*fitsize;
						double temp=fit[index];
						double temp2=data[index];
						double tempx1=((double)j-(double)params[0]-xd)/((double)params[2]*(double)params[2]);
						double tempx2=((double)j-(double)params[0]+xd)/((double)params[2]*(double)params[2]);
						double[] d={amps[1]*fit1[index]*tempx1+amps[2]*fit2[index]*tempx2,
								amps[1]*fit1[index]*tempy1+amps[2]*fit2[index]*tempy2,
								params[2]*(amps[1]*fit1[index]*(tempx1*tempx1+tempy1*tempy1)+amps[2]*fit2[index]*(tempx2*tempx2+tempy2*tempy2)),
								amps[1]*fit1[index]*(xd*tempy1-yd*tempx1)+amps[2]*fit2[index]*(yd*tempx2-xd*tempy2),
								(amps[1]*fit1[index]*(xd*tempx1-yd*tempy1)+amps[2]*fit2[index]*(xd*tempx2-yd*tempy2))/params[4]
								};
						//double[] d={(temp-params[2])*tempx,(temp-params[2])*tempy,1.0,(temp-params[2])/params[3],(temp-params[2])*(tempx+tempy)/(double)params[4]};
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

}
