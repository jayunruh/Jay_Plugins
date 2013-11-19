/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.sf;

public class erffunc{
	public float[] model_erf;

	public erffunc(){
		model_erf=new float[1000];
		for(int i=0;i<1000;i++){
			double r=0.01*(double)i;
			model_erf[i]=(float)sf.erf(r);
		}
	}

	public double getinterperf2(double r,double stdev){
		// here is an "unscaled" erf
		return 2.0*getinterperf(r,stdev)-1.0;
	}

	public double getinterperf(double r,double stdev){
		double rrel=Math.abs(r)/stdev;
		int rp=(int)(rrel*100.0);
		double rem=rrel*100.0-(double)rp;
		double interp=1.0;
		if(rp<999){
			interp=rem*(model_erf[rp+1]-model_erf[rp])+model_erf[rp];
		}
		if(r>=0.0){
			return 0.5*interp+0.5;
		}else{
			return -0.5*interp+0.5;
		}
	}

	public float[] get_func(double startr,int nr,double dr,double stdev){
		float[] temp=new float[nr];
		for(int i=0;i<nr;i++){
			temp[i]=(float)getinterperf(startr+i*dr,stdev);
		}
		return temp;
	}

	public float[] get_func(float[] rvals,double stdev){
		float[] temp=new float[rvals.length];
		for(int i=0;i<rvals.length;i++){
			temp[i]=(float)getinterperf((double)rvals[i],stdev);
		}
		return temp;
	}

}
