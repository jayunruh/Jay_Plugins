/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class pch{
	double dx=0.01;
	public double[] pch1;
	public int klength;
	public int psfflag; // psf flag is 0 for 3D gaussian, 1 for gaussian
	// lorentzian squared, 2 for 2D gaussian
	private double[] cbright,cnums;
	private double cback;
	private int cnlength;
	private boolean ready;

	public pch(int klength1,int psfflag1){
		klength=klength1;
		psfflag=psfflag1;
		ready=false;
	}

	public double getpch(int k,double[] brightnesses,double[] numbers,double background,int nlength){
		if(ready){
			if(cbright==null){
				cbright=new double[brightnesses.length];
				cnums=new double[brightnesses.length];
				ready=false;
			}else{
				for(int i=0;i<brightnesses.length;i++){
					if((brightnesses[i]!=cbright[i]||numbers[i]!=cnums[i])||(nlength!=cnlength||background!=cback)){
						ready=false;
						break;
					}
				}
			}
		}
		if(!ready){
			if(cbright==null){
				cbright=new double[brightnesses.length];
				cnums=new double[brightnesses.length];
			}
			cnlength=nlength;
			cback=background;
			for(int i=0;i<brightnesses.length;i++){
				cbright[i]=brightnesses[i];
				cnums[i]=numbers[i];
			}
			update_pch(brightnesses,numbers,background,nlength);
			ready=true;
		}
		return pch1[k];
	}

	public void update_pch(double[] brightnesses,double[] numbers,double background,int nlength){
		pch1=singlespecies(brightnesses[0],numbers[0],nlength);
		for(int i=1;i<brightnesses.length;i++){
			if(numbers[i]!=0.0){
				double[] temppch2=singlespecies(brightnesses[i],numbers[i],nlength);
				convolute_with(pch1,temppch2);
			}
		}
		if(background!=0.0){
			double[] poissonvals=new double[klength+1];
			for(int i=0;i<=klength;i++){
				poissonvals[i]=poisson(background,i);
			}
			convolute_with(pch1,poissonvals);
		}
	}

	public double[] singlespecies(double brightness,double number,int nlength){
		double[] temppch=new double[klength+1];
		double[] probsingles1=new double[klength+1];
		double[] probsingles2=new double[klength+1];
		double[] probsingles3=new double[klength+1];
		if(psfflag==0){
			for(int i=1;i<=klength;i++){
				probsingles1[i]=p3DG(i,brightness);
			}
		}else{
			if(psfflag==1){
				for(int i=1;i<=klength;i++){
					probsingles1[i]=p2GL(i,brightness);
				}
			}else{
				for(int i=1;i<=klength;i++){
					probsingles1[i]=p2DG(i,brightness);
				}
			}
		}
		// now calculate the 1 particle probability for k = 0 using the
		// normalization condition
		double sum=0.0;
		for(int i=1;i<=klength;i++){
			sum+=probsingles1[i];
		}
		probsingles1[0]=1.0-sum;

		// now calculate the average of the n particle distributions for each k
		// value weighted by
		// the poissonian distribution of n in the focal volume given an average
		// n value
		// the zero particle distribution is 1 for k = 0 and 0 for all other k
		// values
		for(int i=0;i<nlength;i++){
			if(i>0){
				for(int j=0;j<=klength;j++){
					probsingles2[j]=probsingles3[j];
				}
				convolute(probsingles1,probsingles2,probsingles3,klength+1);
			}else{
				for(int j=0;j<=klength;j++){
					probsingles3[j]=0.0;
				}
				probsingles3[0]=1.0;
			}
			for(int j=0;j<=klength;j++){
				temppch[j]+=probsingles3[j]*poisson(number,i);
			}
		}
		return temppch;
	}

	public double p2GL(int k,double e){
		// this function uses equation 15 from Chen et al 1999 to calculate the
		// probability of receiving k counts from a single particle
		// with brightness, i, given a gaussian lorentzian squared point spread
		// function
		double value,integral,x;
		value=gamma(k,(4.0*e)/(Math.PI*Math.PI));
		integral=(value*dx)/2.0;
		x=0.0;
		do{
			x=x+dx;
			value=(1+x*x)*gamma(k,(4.0*e)/(Math.PI*Math.PI*(1+x*x)*(1+x*x)));
			integral=integral+value*dx;
		}while(value>(0.0001*integral));
		return integral*(Math.PI/(2.0*factorial(k)));
	}

	public double p3DG(int k,double e){
		// this function uses equation 16 from Chen et al 1999 to calculate the
		// probability of receiving k counts from a single particle
		// with brightness, i, given a 3D gaussian point spread function
		double value,integral,x;
		value=gamma(k,e);
		integral=(value*dx)/2.0;
		x=0.0;
		do{
			x=x+dx;
			value=gamma(k,(e*Math.exp((-2.0)*x*x)));
			integral=integral+value*dx;
		}while(value>(0.0001*integral));
		return integral*(Math.pow(2.0,1.5)/(Math.sqrt(Math.PI)*factorial(k)));
	}

	public double p2DG(int k,double e){
		// since the volume of the psf doesn't matter, we will choose w0=1
		// for 2D integration, dx*dy=r*dr*dtheta
		double value,integral,x;
		value=dx*poisson(e,k);
		integral=(value*dx)/2.0;
		x=dx;
		do{
			x+=dx;
			value=x*poisson(e*Math.exp(-2.0*x*x),k);
			integral+=value*dx;
		}while(value>(0.0001*integral));
		return 4.0*integral;
	}

	public double poisson(double average,int counts){
		// this function calculates the poisson probability of receiving counts
		// given average counts
		// first calculate the factorial of counts
		if(average>0.0){
			return (Math.pow(average,counts)*(Math.exp(-average)))/factorial(counts);
		}else{
			if(counts==0){
				return 1.0;
			}else{
				return 0.0;
			}
		}
	}

	public double gamma(int counts,double x){
		// this function calculates incomplete gamma function
		int i;
		double gval;
		gval=0.0;
		for(i=0;i<counts;i++){
			gval+=Math.pow(x,i)/factorial(i);
		}
		return factorial(counts-1)*(1-Math.exp(-x)*gval);
	}

	public double factorial(int number){
		double factval;
		int i;
		if(number>0){
			factval=1.0;
			for(i=1;i<=number;i++){
				factval*=(double)i;
			}
			return factval;
		}else{
			return 1.0;
		}
	}

	public void convolute(double[] array1,double[] array2,double[] retarray,int numpts){
		int i,j;
		for(i=0;i<numpts;i++){
			retarray[i]=0.0;
			for(j=0;j<=i;j++){
				retarray[i]+=array1[i-j]*array2[j];
			}
		}
	}

	public void convolute_with(double[] array1,double[] array2){
		// here we convolve array2 into array1
		int i,j;
		double[] temparray1=new double[array1.length];
		System.arraycopy(array1,0,temparray1,0,array1.length);
		for(i=0;i<array1.length;i++){
			array1[i]=0.0;
			for(j=0;j<=i;j++){
				array1[i]+=temparray1[i-j]*array2[j];
			}
		}
	}

	public double[] largefact(int number){
		double[] factval=new double[2];
		int i;
		if(number>0){
			factval[0]=1.0;
			for(i=1;i<=number;i++){
				factval[0]*=(double)i;
				if(factval[0]>1000.0){
					factval[1]+=3;
					factval[0]/=1000.0;
				}
			}
			return factval;
		}else{
			factval[0]=1.0;
			factval[1]=0.0;
			return factval;
		}
	}
}
