/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class besselK{
	// this is a lookup table based class for calculating bessel functions
	// of 0,first, and second order for variables up to 999.8
	public double[] bessel0vals,bessel1vals,bessel2vals;
	public besselI besseli;

	public besselK(){
		besseli=new besselI();
		bessel0vals=new double[10000];
		bessel1vals=new double[10000];
		bessel2vals=new double[10000];
		for(int i=0;i<10000;i++){
			bessel0vals[i]=bessk0(0.1*i);
			bessel1vals[i]=bessk1(0.1*i);
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

	double bessk0(double x){
		// Returns the Bessel function K0(x) for any real x. Adapted from
		// Computer Approximations by Hart 1978.
		double y,ans;

		if(x<=2.0){
			y=x*x/4.0;
			ans=(-Math.log(x/2.0)*besseli.besselval(x,0))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5))))));
		}else{
			y=2.0/x;
			ans=(Math.exp(-x)/Math.sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
		}
		return ans;
	}

	double bessk1(double x){
		// Returns the Bessel function Y1(x) for any real x. Adapted from
		// Computer Approximations by Hart 1978.
		double y,ans;

		if(x<=2.0){
			y=x*x/4.0;
			ans=(Math.log(x/2.0)*besseli.besselval(x,1))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))));
		}else{
			y=2.0/x;
			ans=(Math.exp(-x)/Math.sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
		}
		return ans;
	}

}
