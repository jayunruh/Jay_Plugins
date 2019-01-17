/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class besselI{
	// this is a lookup table based class for calculating bessel functions
	// of 0,first, and second order for variables up to 999.8
	public double[] bessel0vals,bessel1vals,bessel2vals;

	public besselI(){
		bessel0vals=new double[10000];
		bessel1vals=new double[10000];
		bessel2vals=new double[10000];
		for(int i=0;i<10000;i++){
			bessel0vals[i]=bessi0(0.1*i);
			bessel1vals[i]=bessi1(0.1*i);
			// bessel2vals[i]=bessel2(0.1*(double)i);
		}
	}

	public double besselval(double x,int order){
		// here we look up the bessel functions from the table calculated at the
		// beginning
		int index;
		double prevval,nextval,fraction,dindex;
		dindex=x*10.0;
		index=(int)dindex;
		if(index>9998){
			return 0.0;
		}
		if(order==0){
			prevval=bessel0vals[index];
			nextval=bessel0vals[index+1];
		}else{
			if(order==1){
				prevval=bessel1vals[index];
				nextval=bessel1vals[index+1];
			}else{
				prevval=bessel2vals[index];
				nextval=bessel2vals[index+1];
			}
		}
		fraction=dindex-index;
		return(prevval+fraction*(nextval-prevval));
	}

	double bessi0(double x){
		// Returns the modified Bessel function I0(x) for any real x. Adapted
		// from Computer Approximations by Hart 1978.
		double y,ax,ans;
		if((ax=Math.abs(x))<3.75){
			y=x/3.75;
			y*=y;
			ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
		}else{
			y=3.75/ax;
			ans=(Math.exp(ax)/Math.sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
		}
		return ans;
	}

	double bessi1(double x){
		// Returns the Bessel function I1(x) for any real x. Adapted from
		// Computer Approximations by Hart 1978.
		double ax,ans;
		double y;

		if((ax=Math.abs(x))<3.75){
			y=x/3.75;
			y*=y;
			ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
		}else{
			y=3.75/ax;
			ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
			ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
			ans*=(Math.exp(ax)/Math.sqrt(ax));
		}
		return x<0.0?-ans:ans;
	}

}
