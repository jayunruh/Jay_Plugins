/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

import jalgs.jsort;

import java.util.Random;

public class rngs{
	public Random jrnd;
	/*
	 * here are some random number generators if you want it to construct from a
	 * previous seed use the first constructor Copyright Jay Unruh Stowers
	 * Institute for Medical Research 4/25/08
	 */

	// static variables for gasdev
	int prev;
	double oldg;

	// static variables for poidev
	double oldavg,g;

	public rngs(int seed1){
		// seed=seed1;
		jrnd=new Random(seed1);
		prev=0;
		oldavg=-1.0;
	}

	public rngs(){
		// randomize_seed();
		jrnd=new Random();
		prev=0;
		oldavg=-1.0;
	}

	public double gasdev(double avg,double stdev){
		// this is an adaptation of the gaussian random number generator from
		// numerical recipes
		double a,rr,v1,v2;
		// if we don't already have a number from the previous run:
		if(prev==0){
			do{
				v1=2.0*jrnd.nextDouble()-1.0;
				v2=2.0*jrnd.nextDouble()-1.0;
				rr=v1*v1+v2*v2;
			}while(rr>=1.0||rr==0.0);
			a=Math.sqrt(-2.0*Math.log(rr)/rr);
			oldg=v1*a;
			prev=1;
			return v2*a*stdev+avg;
		}else{
			prev=0;
			return oldg*stdev+avg;
		}
	}

	public int poidev(double avg){
		// generates poisson random numbers as integers given mean value avg
		double em,t;
		// generate exponentially distributed event
		// times and the number of these required to exceed avg
		// gives us the poisson random number
		if(avg>1000.0){
			return (int)gasdev(avg,Math.sqrt(avg));
		}
		if(avg!=oldavg){
			oldavg=avg;
			g=Math.exp(-avg);
		}
		em=-1.0;
		t=1.0;
		// can multiply uniform deviates rather than adding exponential deviates
		do{
			em+=1.0;
			t*=jrnd.nextDouble();
		}while(t>g);
		return (int)em;
	}

	public double expdev(double avg){
		return -avg*Math.log(jrnd.nextDouble());
	}

	public double unidev(double upper,double lower){
		return jrnd.nextDouble()*(upper-lower)+lower;
	}

	public int random_integerize(float val){
		if(val%1.0f==0.0f) return (int)val;
		int add=jrnd.nextBoolean()?1:0;
		return (int)val+add;
	}

	public double analogdev(double avg,double gain,double offset,double readstdev){
		int photons=poidev(avg);
		double signal=0.0;
		for(int i=0;i<photons;i++){
			signal+=expdev(gain);
		}
		return signal+gasdev(offset,readstdev);
	}

	public double analogdev(int photons,double gain,double offset,double readstdev){
		double signal=0.0;
		for(int i=0;i<photons;i++){
			signal+=expdev(gain);
		}
		return signal+gasdev(offset,readstdev);
	}

	public double multianalogdev(double avg,float[] gain,float[] amp,double offset,double readstdev){
		int photons=poidev(avg);
		return multianalogdev(photons,gain,amp,offset,readstdev);
	}

	public double multianalogdev(int photons,float[] gain,float[] amp,double offset,double readstdev){
		double signal=0.0;
		for(int i=0;i<photons;i++){
			if(amp[0]==1.0f){
				signal+=expdev(gain[0]);
			}else{
				double temp=unidev(1.0,0.0);
				int counter=-1;
				double temp2=0.0;
				do{
					counter++;
					temp2+=amp[counter];
				}while(temp>temp2);
				signal+=expdev(gain[counter]);
			}
		}
		return signal+gasdev(offset,readstdev);
	}

	public int arbdev(double[] dist){
		// integrate the distribution until the integral is greater than a
		// uniform random number between 0 and 1
		double integral=0.0;
		int counter=-1;
		double dumdouble=unidev(1.0,0.0);
		do{
			counter++;
			integral+=dist[counter];
		}while(counter<dist.length&&integral<dumdouble);
		if(counter>=dist.length){
			counter=dist.length-1;
		}
		return counter;
	}

	public double arbdev(double[] dist,double[] xvals){
		// integrate the distribution until the integral is greater than a
		// uniform random number between 0 and 1
		double integral=0.0;
		int counter=-1;
		double dumdouble=unidev(1.0,0.0);
		do{
			counter++;
			integral+=dist[counter];
		}while(counter<dist.length&&integral<dumdouble);
		if(counter>=dist.length){
			counter=dist.length-1;
		}
		return xvals[counter];
	}

	public double[] random_circle(double radius){
		double[] temp={unidev(radius,-radius),unidev(radius,-radius)};
		double length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]);
		while(length>radius){
			double[] temp2={unidev(radius,-radius),unidev(radius,-radius)};
			temp=temp2;
			length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]);
		}
		return temp;
	}
	
	public double[] random_sphere(double radius) {
		double[] vec = random_3D_vector();
    	vec[0]*=radius;
    	vec[1]*=radius;
    	vec[2]*=radius;
    	return vec;
	}

	public double[] random_3D_vector(){
		double[] temp={unidev(1.0,-1.0),unidev(1.0,-1.0),unidev(1.0,-1.0)};
		double length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
		while(length>1.0){
			double[] temp2={unidev(1.0,-1.0),unidev(1.0,-1.0),unidev(1.0,-1.0)};
			temp=temp2;
			length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
		}
		temp[0]/=length;
		temp[1]/=length;
		temp[2]/=length;
		return temp;
	}

	public double[] random_2D_vector(){
		double[] temp={unidev(1.0,-1.0),unidev(1.0,-1.0)};
		double length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]);
		while(length>1.0){
			double[] temp2={unidev(1.0,-1.0),unidev(1.0,-1.0)};
			temp=temp2;
			length=Math.sqrt(temp[0]*temp[0]+temp[1]*temp[1]);
		}
		temp[0]/=length;
		temp[1]/=length;
		return temp;
	}

	public float[] random_arrayf(int length,float min,float max){
		float[] temp=new float[length];
		for(int i=0;i<length;i++){
			temp[i]=(float)unidev(max,min);
		}
		return temp;
	}

	public double[] random_arrayd(int length,float min,float max){
		double[] temp=new double[length];
		for(int i=0;i<length;i++){
			temp[i]=unidev(max,min);
		}
		return temp;
	}

	public byte[] random_arrayb(int length){
		byte[] temp=new byte[length];
		jrnd.nextBytes(temp);
		return temp;
	}

	public int[] random_order(int length){
		float[] temp=random_arrayf(length,0.0f,1.0f);
		return jsort.get_javasort_order(temp);
	}
	
	public float[] random_fractions(int ncat){
		float[] fractions=new float[ncat];
		float sum=0.0f;
		for(int i=0;i<ncat;i++){
			fractions[i]=(float)unidev(1.0,0.0);
			sum+=fractions[i];
		}
		for(int i=0;i<ncat;i++) fractions[i]/=sum;
		return fractions;
	}

}
