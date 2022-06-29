/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class binning_function{
	public double[] b3func=new double[89];
	double three_1_5;
	public double rtemp;
	public double r;

	public binning_function(double r1){
		three_1_5=Math.pow(3.0,1.5);
		r=r1;
		rtemp=1.0/(r*r);
		initialize_function();
	}

	private void initialize_function(){
		double step=1.05;
		int counter=0;
		for(double T_td=0.001;T_td<10000.0;T_td*=1.2){
			double integral=0.0;
			for(double i=0.00001;i<T_td;i*=step){
				double istep=(i*step-i/step)/2.0;
				if(i==0.00001){
					istep=0.00001;
				}
				for(double j=0.00001;j<T_td-i;j*=step){
					double jstep=(j*step-j/step)/2.0;
					if(j==0.00001){
						jstep=0.00001;
					}
					integral+=(T_td-i-j)*g3(i,j)*istep*jstep;
				}
			}
			b3func[counter]=6.0*integral;
			counter++;
		}
	}

	public double b3(double taud,double T){
		double T_taud=T/taud;
		double x=Math.log(T_taud*1000.0)/Math.log(1.2);
		if(x<0.0){
			return T*T*T;
		}
		if(x>=88.0){
			return 0.0;
		}
		int xprev=(int)x;
		double rem=x-xprev;
		double interp=rem*(b3func[xprev+1]-b3func[xprev])+b3func[xprev];
		return taud*taud*taud*interp;
	}

	public double b2(double taud,double T){
		double T_td=T/taud;
		double b=Math.sqrt(r*r-1);
		double c=Math.sqrt(T_td+r*r);
		double temp=((4.0*r)/b)*(2.0*r*b-2.0*b*c+(T_td+1.0)*Math.log(((r+b)*(b-c))/((b-r)*(b+c))));
		return taud*taud*temp;
	}

	public double g3(double x2,double x3){
		return three_1_5/((4.0*x2*x3+4.0*(x2+x3)+3.0)*Math.sqrt(4.0*x2*x3*rtemp*rtemp+4.0*(x2+x3)*rtemp+3.0));
	}
}
