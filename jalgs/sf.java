/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

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
		return (double)k*Math.log(xx)-xx-lnfact(k);
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
										r=r*(a+(double)k-1.0)/((double)k)/(b+(double)k-1.0)*x;
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

}
