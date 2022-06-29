/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class bessely{
	// this is a lookup table based class for calculating bessel functions
	// of 0,first, and second order for variables up to 999.8
	public double[] bessel0vals,bessel1vals,bessel2vals;
	public bessel besselj;

	public bessely(){
		besselj=new bessel();
		bessel0vals=new double[10000];
		bessel1vals=new double[10000];
		bessel2vals=new double[10000];
		for(int i=0;i<10000;i++){
			bessel0vals[i]=bessy0(0.1*i);
			bessel1vals[i]=bessy1(0.1*i);
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

	double bessy0(double x){
		// Returns the Bessel function Y0(x) for any real x. Adapted from
		// Computer Approximations by Hart 1978.
		double z,xx,y,ans,ans1,ans2;

		if(x<8.0){
			y=x*x;
			ans1=-2957821389.0+y*(7062834065.0+y*(-512359803.6+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
			ans2=40076544269.0+y*(745249964.8+y*(7189466.438+y*(47447.26470+y*(226.1030244+y*1.0))));
			ans=(ans1/ans2)+0.636619772*besselj.besselval(x,0)*Math.log(x);
		}else{
			z=8.0/x;
			y=z*z;
			xx=x-0.785398164;
			ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
			ans2=-0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6+y*(-0.934945152e-7))));
			ans=Math.sqrt(0.636619772/x)*(Math.sin(xx)*ans1+z*Math.cos(xx)*ans2);
		}
		return ans;
	}

	double bessy1(double x){
		// Returns the Bessel function Y1(x) for any real x. Adapted from
		// Computer Approximations by Hart 1978.
		double z,xx,y,ans,ans1,ans2;
		if(x<8.0){
			y=x*x;
			ans1=x*(-0.4900604943e13+y*(0.1275274390e13+y*(-0.5153438139e11+y*(0.7349264551e9+y*(-0.4237922726e7+y*0.8511937935e4)))));
			ans2=0.2499580570e14+y*(0.4244419664e12+y*(0.3733650367e10+y*(0.2245904002e8+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
			ans=(ans1/ans2)+0.636619772*(besselj.besselval(x,1)*Math.log(x)-1.0/x);
		}else{
			z=8.0/x;
			y=z*z;
			xx=x-2.356194491;
			ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
			ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
			ans=Math.sqrt(0.636619772/x)*(Math.sin(xx)*ans1+z*Math.cos(xx)*ans2);
		}
		return ans;
	}

}
