/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.awt.Rectangle;

import jalgs.jsim.rngs;

public class sf{

	public static double gammaln(float x){
		return gammaln((double)x);
	}

	public static double gammaln(double x){
		// adapted from numerical recipes
		double[] coefficients={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
		double y=x;
		double a=x+5.5;
		a-=(x+0.5)*Math.log(a);
		double z=1.000000000190015;
		for(int j=0;j<=5;j++)
			z+=coefficients[j]/++y;
		return -a+Math.log(2.5066282746310005*z/x);
	}

	public static double gamma(float xx){
		return Math.exp(gammaln(xx));
	}

	public static double lnfact(int k){
		return gammaln((double)(k+1));
	}

	public static double lnpoisson(double xx,int k){
		return k*Math.log(xx)-xx-lnfact(k);
	}

	public static double poisson(double xx,int k){
		if(xx>0.0){
			return Math.exp(lnpoisson(xx,k));
		}else{
			if(k==0){
				return 1.0;
			}else{
				return 0.0;
			}
		}
	}

	public static double conf_hypgeo(double a10,double b10,double x10){
		// this is adapted from "computation of special functions" by Zhang and
		// Jin
		double hg=0.0;
		double a=a10;
		double b=b10;
		double x=x10;
		double a0=a;
		double a1=a;
		double x0=x;
		double r,rg,ta,tb,tba,xg,sum1,sum2,r1,r2,hg1,hg2,y0,y1;
		int m,nl,la,n,i,j,k;
		if(b==0.0||b==-Math.abs((int)b)){
			hg=1.0e300;
		}else{
			if(a==0.0||x==0.0){
				hg=1.0;
			}else{
				if(a==-1.0){
					hg=1-x/b;
				}else{
					if(a==b){
						hg=Math.exp(x);
					}else{
						if((a-b)==1.0){
							hg=(1.0+x/b)*Math.exp(x);
						}else{
							if(a==1.0&&b==2.0){
								hg=(Math.exp(x)-1.0)/x;
							}else{
								if(a==(int)a&&a<0.0){
									m=(int)(-a);
									r=1.0;
									hg=1.0;
									for(k=1;k<=m;k++){
										r=r*(a+k-1.0)/(k)/(b+k-1.0)*x;
										hg+=r;
									}
								}
							}
						}
					}
				}
			}
		}
		if(hg!=0.0){
			return hg;
		}
		if(x<0.0){
			a=b-a;
			a0=a;
			x=Math.abs(x);
		}
		nl=0;
		la=0;
		y0=0.0;
		y1=0.0;
		if(a>=2.0){
			nl=1;
			la=(int)a;
			a=a-la-1.0;
		}
		for(n=0;n<=nl;n++){
			if(a0>=2.0)
				a+=1.0;
			if(x<=(30.0+Math.abs(b))||a<0.0){
				hg=1.0;
				rg=1.0;
				for(j=1;j<=500;j++){
					rg=rg*(a+j-1.0)/(j*(b+j-1.0))*x;
					hg+=rg;
					if(Math.abs(rg/hg)<1.0E-15){
						break;
					}
				}
			}else{
				ta=gammafor((float)a);
				tb=gammafor((float)b);
				xg=b-a;
				tba=gammafor((float)xg);
				sum1=1.0;
				sum2=1.0;
				r1=1.0;
				r2=1.0;
				for(i=1;i<=8;i++){
					r1=-r1*(a-i-1.0)*(a-b+i)/(x*i);
					r2=-r2*(b-a+i-1.0)*(a-i)/(x*i);
					sum1+=r1;
					sum2+=r2;
				}
				hg1=tb/tba*Math.pow(x,-a)*Math.cos(Math.PI*a)*sum1;
				hg2=tb/ta*Math.exp(x)*Math.pow(x,a-b)*sum2;
				hg=hg1+hg2;
			}
			if(n==0)
				y0=hg;
			if(n==1)
				y1=hg;
		}
		if(a0>=2.0){
			for(i=1;i<=(la-1);i++){
				hg=((2.0*a-b+x)*y1+(b-a)*y0)/a;
				y0=y1;
				y1=hg;
				a+=1.0;
			}
		}
		if(x0>=2.0)
			hg*=Math.exp(x0);
		a=a1;
		x=x0;
		return hg;
	}

	public static double gammafor(double x){
		// this is adapted from "computation of special functions" by Zhang and
		// Jin
		double[] g={1.0,0.5772156649015329,-0.6558780715202538,-0.420026350340952e-1,0.1665386113822915,-.421977345555443e-1,-0.96219715278770e-2,0.72189432466630e-2,-0.11651675918591e-2,
				-0.2152416741149e-3,0.1280502823882e-3,-0.201348547807e-4,-0.12504934821e-5,0.11330272320e-5,-0.2056338417e-6,0.61160950e-8,0.50020075e-8,-0.11812746e-8,0.1043427e-9,0.77823e-11,
				-0.36968e-11,0.51e-12,-0.206e-13,-0.54e-14,0.14e-14,0.1e-15};
		double ga=1.0;
		double z,r,gr;
		int m1,k,m;
		r=0.0;
		if(x==(int)x){
			if(x>0.0){
				ga=1.0;
				m1=(int)x-1;
				for(k=2;k<=m1;k++){
					ga*=k;
				}
			}else{
				ga=1.0E300;
			}
		}else{
			if(Math.abs(x)>1.0){
				z=Math.abs(x);
				m=(int)z;
				r=1.0;
				for(k=1;k<=m;k++){
					r*=(z-k);
				}
				z-=m;
			}else{
				z=x;
			}
			gr=g[25];
			for(k=25;k>=1;k--){
				gr=gr*z+g[k-1];
			}
			ga=1.0/(gr*z);
			if(Math.abs(x)>1.0){
				ga*=r;
				if(x<0.0)
					ga=-Math.PI/(x*ga*Math.sin(Math.PI*x));
			}
		}
		return ga;
	}

	public static double conf_hypgeo_integral(double a,double b,double z){
		// here b must be greater than a
		return Math.exp(lnconf_hypgeo_integral(a,b,z));
	}

	public static double lnconf_hypgeo_integral(double a,double b,double z){
		// here b must be greater than a
		double lnmult=gammaln(b)-gammaln(a)-gammaln(b-a);
		double integral=conf_hypgeo_integral_only(a,b,z);
		return lnmult+Math.log(integral);
	}

	public static double conf_hypgeo_integral_only(double a,double b,double z){
		double integral=0.0;
		for(double u=0.01;u<1.0;u+=0.01){
			integral+=Math.exp(z*u)*Math.pow(u,a-1.0)*Math.pow(1.0-u,b-a-1.0)*0.01;
		}
		return integral;
	}

	public static double erfc(double x){
		// complimentary error function adapted from numerical recipes
		double a,b,retval;
		b=Math.abs(x);
		a=1.0/(1.0+0.5*b);
		retval=Math.exp(-b*b-1.26551223+a*(1.00002368+a*(0.37409196+a*(0.09678418+a*(-0.18628806+a*(0.27886807+a*(-1.13520398+a*(1.48851587+a*(-0.82215223+a*0.17087277)))))))));
		retval*=a;
		if(x>=0.0)
			return retval;
		else
			return 2.0-retval;
	}

	public static double erf(double x){
		return 1.0-erfc(x);
	}
	
	public static double dixonqlim(int dof){
		double[] lims={0.97,0.829,0.71,0.625,0.568,0.526,0.493,0.466,0.444,0.426,0.41,0.396,0.384,0.374,0.365,0.356,0.349};
		if(dof<3) return 0.0;
		if(dof>19) return 0.9592*Math.pow(dof,-0.344);
		return lims[dof-3];
	}
	
	public static int dixonqtest(float[] data,boolean max){
		float[] temp=find2minmax(data);
		float min1=temp[0]; float min2=temp[1]; float max1=temp[3]; float max2=temp[4];
		int minind=(int)temp[2]; int maxind=(int)temp[5];
		if(max){
			float q=(max1-max2)/(max1-min1);
			double qlim=dixonqlim(data.length);
			if(q>qlim) return maxind;
			else return -1;
		} else {
			float q=(min2-min1)/(max1-min1);
			double qlim=dixonqlim(data.length);
			if(q>qlim) return minind;
			else return -1;
		}
	}
	
	public static float[] find2minmax(float[] data){
		float max1=data[0]; float max2=data[1]; int maxind=0;
		if(max2>max1){
			max1=data[1]; max2=data[0]; maxind=1;
		}
		float min1=data[0]; float min2=data[1]; int minind=0;
		if(min2<min1){
			min1=data[1]; min2=data[0]; minind=1;
		}
		for(int i=0;i<data.length;i++){
			if(data[i]>max1){
				max2=max1; max1=data[i]; maxind=i;
			}
			if(data[i]<min1){
				min2=min1; min1=data[i]; minind=i;
			}
		}
		return new float[]{min1,min2,minind,max1,max2,maxind};
	}
	
	public static float[] findminmax(float[] data){
		float min=data[0]; int minind=0;
		float max=data[0]; int maxind=0;
		for(int i=1;i<data.length;i++){
			if(data[i]<min){min=data[i]; minind=i;}
			if(data[i]>max){max=data[i]; maxind=i;}
		}
		return new float[]{min,minind,max,maxind};
	}
	
	public static float grubbsglim(int N,float alpha,boolean twotailed){
		int dof=N-2;
		float p=alpha/(2*N);
		float tlim=(float)(new jdist()).tLim(dof, p, twotailed);
		return (N-1)*(float)Math.sqrt(tlim*tlim/(N*(N-2)+tlim*tlim));
	}
	
	public static int[] grubbstest(float[] data,float prob){
		float avg=jstatistics.getstatistic("Avg",data,null);
		float stdev=jstatistics.getstatistic("stdev",data,null);
		float[] temp=findminmax(data);
		float g=(temp[2]-avg)/stdev;
		float glim=grubbsglim(data.length,1.0f-prob,true);
		int[] indices={(int)temp[1],(int)temp[3]};
		if(g<=glim) indices[1]=-1;
		g=(avg-temp[0])/stdev;
		if(g<=glim) indices[0]=-1;
		return indices;
	}
	
	public static int[] mcoutliers(float[] data,int ntrials,float prob){
		//here we do a monte carlo test for outliers
		float avg=jstatistics.getstatistic("Avg",data,null);
		float stdev=jstatistics.getstatistic("stdev",data,null);
		float[] temp=findminmax(data);
		int countabove=0;
		int countbelow=0;
		rngs random=new rngs();
		for(int i=0;i<ntrials;i++){
			boolean above=false;
			boolean below=false;
			for(int j=0;j<data.length;j++){
				float temp2=(float)random.gasdev(avg,stdev);
				if(!above && temp2>temp[2]){
					countabove++;
					above=true;
				}
				if(!below && temp2<temp[0]){
					countbelow++;
					below=true;
				}
				if(above && below) break;
			}
		}
		float fabove=(float)countabove/(float)ntrials;
		float fbelow=(float)countbelow/(float)ntrials;
		int[] indices={(int)temp[1],(int)temp[3],countbelow,countabove};
		if(fabove>(1.0f-prob)) indices[1]=-1;
		if(fbelow>(1.0f-prob)) indices[0]=-1;
		return indices;
	}
	
	public static Object[] match_x_avg_diff(float[] x1,float[] data1,float[] x2,float[] data2){
		//this matches the average x value from two data sets by reducing the higher one and interpolating
		float avgx1=jstatistics.getstatistic("Avg",x1,null);
		float avgx2=jstatistics.getstatistic("Avg",x2,null);
		float avg1=jstatistics.getstatistic("Avg",data1,null);
		float avg2=jstatistics.getstatistic("Avg",data2,null);
		float sem1=jstatistics.getstatistic("sterr",data1,null);
		float sem2=jstatistics.getstatistic("sterr",data2,null);
		if(avg2==avg1) return new Object[]{new float[]{avg1,avg2,data1.length,data2.length,sem1,sem2,0.5f},x1,data1,x2,data2};
		//first sort the data sets in increasing order
		float[] s1=x1.clone();
		int[] ord1=jsort.javasort_order(s1);
		float[] s2=x2.clone();
		int[] ord2=jsort.javasort_order(s2);
		float[] datas1=new float[data1.length];
		for(int i=0;i<data1.length;i++) datas1[i]=data1[ord1[i]];
		float[] datas2=new float[data2.length];
		for(int i=0;i<data2.length;i++) datas2[i]=data2[ord2[i]];
		if(avgx2>avgx1){
			float prevavg2=avg2;
			float prevavgx2=avgx2;
			int n=data2.length;
			while(avgx2>avgx1 && n>1){
				prevavg2=avg2;
				prevavgx2=avgx2;
				avgx2=(avgx2*(float)n-s2[n-1])/(float)(n-1);
				avg2=(avg2*(float)n-datas2[n-1])/(float)(n-1);
				n--;
			}
			if(n==1) return null;
			float[] news2=new float[n+1]; System.arraycopy(s2,0,news2,0,n+1);
			float[] newdatas2=new float[n+1]; System.arraycopy(datas2,0,newdatas2,0,n+1);
			//interpolate between the two avgs
			float interpfrac=(avgx1-avgx2)/(prevavgx2-avgx2);
			float interp=avg2+(prevavg2-avg2)*interpfrac;
			float sterr1=jstatistics.getstatistic("sterr",data1,null);
			float sterr2=jstatistics.getstatistic("sterr",newdatas2,datas2.length,1,new Rectangle(0,0,n,1),null);
			float prevsterr2=jstatistics.getstatistic("sterr",datas2,datas2.length,1,new Rectangle(0,0,n+1,1),null);
			float interpsterr=sterr2+(prevsterr2-sterr2)*interpfrac;
			double prevp=(new jdist()).tCumProbTwoSampEqVar(data1.length,n+1,avg1,prevavg2,sterr1,prevsterr2,0.0,true);
			double p=(new jdist()).tCumProbTwoSampEqVar(data1.length,n,avg1,avg2,sterr1,sterr2,0.0,true);
			double interpp=p+(prevp-p)*interpfrac;
			float[] retvals={avg1,interp,data1.length,n,sterr1,interpsterr,(float)interpp};
			return new Object[]{retvals,s1,datas1,news2,newdatas2};
		} else {
			Object[] temp=match_x_avg_diff(x2,data2,x1,data1);
			float[] retvals=(float[])temp[0];
			float[] newretvals={retvals[1],retvals[0],retvals[3],retvals[2],retvals[5],retvals[4],retvals[6]};
			return new Object[]{newretvals,temp[3],temp[4],temp[1],temp[2]};
		}
	}

}
