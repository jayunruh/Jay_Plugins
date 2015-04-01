/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class bessel{
	// this is a lookup table based class for calculating bessel functions
	// of 0,first, and second order for variables up to 999.8
	public double[] bessel0vals,bessel1vals,bessel2vals;

	public bessel(){
		bessel0vals=new double[10000];
		bessel1vals=new double[10000];
		bessel2vals=new double[10000];
		for(int i=0;i<10000;i++){
			bessel0vals[i]=bessel0(0.1*i);
			bessel1vals[i]=bessel1(0.1*i);
			bessel2vals[i]=bessel2(0.1*i);
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

	public double bessel0(double x)
	// Returns the Bessel function J0(x) for any real x. Adapted from Computer
	// Approximations by Hart 1978.
	{
		double ax,z;
		double xx,y,ans,ans1,ans2; // Accumulate polynomials in double
		// precision.
		if((ax=Math.abs(x))<8.0){ // Direct rational function fit.
			y=x*x;
			ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
			ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
			ans=ans1/ans2;
		}else{ // Fitting function.
			z=8.0/ax;
			y=z*z;
			xx=ax-0.785398164;
			ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
			ans2=-0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
			ans=Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
		}
		return ans;
	}

	public double bessel1(double x)
	// Returns the Bessel function J1(x) for any real x. Adapted from Computer
	// Approximations by Hart 1978.
	{
		double ax,z;
		double xx,y,ans,ans1,ans2; // Accumulate polynomials in double
		// precision.
		if((ax=Math.abs(x))<8.0){ // Direct rational approximation.
			y=x*x;
			ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
			ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
			ans=ans1/ans2;
		}else{ // Fitting function
			z=8.0/ax;
			y=z*z;
			xx=ax-2.356194491;
			ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
			ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
			ans=Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
			if(x<0.0)
				ans=-ans;
		}
		return ans;
	}

	public double bessel2(double x)
	// note that this is not as accurate as the functions above
	{
		int k;
		double funcval;
		funcval=0.0;
		if(x<40.0){
			for(k=0;k<91;k++){
				funcval+=intpow(((-0.25)*x*x),k)/(factorial(k)*factorial(k)*(k+1.0)*(k+2.0));
			}
		}else{
			funcval=Math.sqrt(2.0/(Math.PI*x))*Math.cos(x-1.25*Math.PI);
		}
		return funcval*((x*x)/4.0);
	}

	public double factorial(int x){
		int i;
		double funcval;
		if(x>0){
			funcval=1;
			for(i=1;i<=x;i++){
				funcval*=i;
			}
		}else{
			funcval=1.0;
		}
		return funcval;
	}

	public double intpow(double x,int power){
		// calculates a double value to an integer power (faster than pow(x,y))
		int i;
		double funcval;
		if(power>0){
			funcval=x;
			for(i=2;i<=power;i++){
				funcval*=x;
			}
			return funcval;
		}else{
			return 1.0;
		}
	}
}
