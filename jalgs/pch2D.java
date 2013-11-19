/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class pch2D{
	int klengtha,klengthb,psfflag; // psf flag is 0 for one photon and 1 for two
								   // photon and 2 for 2D diffusion
	boolean ready;
	pch pchfunc;
	double[][] pch1;
	double[] cbrighta,cbrightb,cnums;
	double cbacka,cbackb;
	int cnlength;

	public pch2D(int klength1,int klength2,int psfflag1){
		klengtha=klength1;
		klengthb=klength2;
		psfflag=psfflag1;
		ready=false;
		pchfunc=new pch(klengtha+klengthb,psfflag);
	}

	public double getpch(int ka,int kb,double[] brightnessesa,double[] brightnessesb,double[] numbers,double backgrounda,double backgroundb,int nlength){
		if(ready){
			if(cbrighta==null){
				cbrighta=new double[brightnessesa.length];
				cbrightb=new double[brightnessesa.length];
				cnums=new double[brightnessesa.length];
				ready=false;
			}else{
				for(int i=0;i<brightnessesa.length;i++){
					if(((brightnessesa[i]!=cbrighta[i]||brightnessesb[i]!=cbrightb[i])||numbers[i]!=cnums[i])||(nlength!=cnlength||(backgrounda!=cbacka||backgroundb!=cbackb))){
						ready=false;
						break;
					}
				}
			}
		}
		if(!ready){
			if(cbrighta==null){
				cbrighta=new double[brightnessesa.length];
				cbrightb=new double[brightnessesa.length];
				cnums=new double[brightnessesa.length];
			}
			cnlength=nlength;
			cbacka=backgrounda;
			cbackb=backgroundb;
			for(int i=0;i<brightnessesa.length;i++){
				cbrighta[i]=brightnessesa[i];
				cbrightb[i]=brightnessesb[i];
				cnums[i]=numbers[i];
			}
			update_pch(brightnessesa,brightnessesb,numbers,backgrounda,backgroundb,nlength);
			ready=true;
		}
		return pch1[ka][kb];
	}

	public void update_pch(double[] brightnessesa,double[] brightnessesb,double[] numbers,double backgrounda,double backgroundb,int nlength){
		pch1=singlespecies_2D(brightnessesa[0],brightnessesb[0],numbers[0],nlength);
		for(int i=1;i<brightnessesa.length;i++){
			if(numbers[i]!=0.0){
				double[][] temppch2=singlespecies_2D(brightnessesa[i],brightnessesb[i],numbers[i],nlength);
				convolute_2D_with(pch1,temppch2);
			}
		}
		if(backgrounda!=0.0){
			double[][] poissonvals=new double[klengtha+1][klengthb+1];
			for(int i=0;i<=klengtha;i++){
				poissonvals[i][0]=pchfunc.poisson(backgrounda,i);
			}
			convolute_2D_with(pch1,poissonvals);
		}
		if(backgroundb!=0.0){
			double[][] poissonvals=new double[klengtha+1][klengthb+1];
			for(int i=0;i<=klengthb;i++){
				poissonvals[0][i]=pchfunc.poisson(backgroundb,i);
			}
			convolute_2D_with(pch1,poissonvals);
		}
	}

	public double[][] singlespecies_2D(double brightnessa,double brightnessb,double number,int nlength){
		double[][] temppch=new double[klengtha+1][klengthb+1];
		double[] pchsingle=pchfunc.singlespecies(brightnessa+brightnessb,number,nlength);
		for(int i=0;i<=klengtha;i++){
			for(int j=0;j<=klengthb;j++){
				temppch[i][j]=binomial(i,i+j,brightnessa/(brightnessa+brightnessb))*pchsingle[i+j];
			}
		}
		return temppch;
	}

	public double binomial(int number_true,int trials,double true_prob){
		return (pchfunc.factorial(trials)/(pchfunc.factorial(number_true)*pchfunc.factorial(trials-number_true)))*Math.pow(true_prob,number_true)*Math.pow(1.0-true_prob,trials-number_true);
	}

	public void convolute_2D(double[][] matrix1,double[][] matrix2,double[][] retmatrix){
		int length1=matrix1.length;
		int length2=matrix1[0].length;
		for(int i=0;i<length1;i++){
			for(int j=0;j<length2;j++){
				retmatrix[i][j]=0.0;
				for(int k=0;k<=i;k++){
					for(int l=0;l<=j;l++){
						retmatrix[i][j]=retmatrix[i][j]+matrix1[k][l]*matrix2[i-k][j-l];
					}
				}
			}
		}
	}

	public void convolute_2D_with(double[][] matrix1,double[][] matrix2){
		int length1=matrix1.length;
		int length2=matrix1[0].length;
		double[][] tempmatrix=new double[length1][length2];
		for(int i=0;i<length1;i++){
			System.arraycopy(matrix1[i],0,tempmatrix[i],0,length2);
		}
		for(int i=0;i<length1;i++){
			for(int j=0;j<length2;j++){
				matrix1[i][j]=0.0;
				for(int k=0;k<=i;k++){
					for(int l=0;l<=j;l++){
						matrix1[i][j]+=tempmatrix[k][l]*matrix2[i-k][j-l];
					}
				}
			}
		}

	}
}
