/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class formula_parser{
	public String formula;
	public double[][] coefficients;
	public boolean[][] isterm;
	public int[][] nummult;
	public int numadd,numsub;

	// this is a simple formula parser that only allows for expanded polynomials
	// everything to the right of a divider is considered as denominator
	// therefore there must always be only one divider per polynomial term
	public formula_parser(String formula1,String[] paramsnames){
		// start by removing the equals sign if necessary
		if(formula1.charAt(0)=='='){
			formula=formula1.substring(1,formula1.length());
		}else{
			formula=formula1.toString();
		}
		// if the term starts with a -, add a 1.0 at the beginning
		if(formula.charAt(0)=='-'){
			String tempformula="1.0";
			tempformula.concat(formula);
			formula=tempformula;
		}
		String[][] terms=new String[10][10];
		numadd=0;
		numsub=0;
		// now parse the add components
		String[] temp=parseadd(formula);
		numadd=temp.length;
		if(numadd>5){
			numadd=5;
		}
		// now the subtract components
		int counter=0;
		for(int i=0;i<numadd;i++){
			terms[i][0]=temp[i];
			String[] temp2=parsesubtract(terms[i][0]);
			terms[i][0]=temp2[0];
			if(temp2.length>1){
				for(int j=1;j<temp2.length;j++){
					if(counter>=5){
						break;
					}
					terms[counter+5][0]=temp2[j];
					counter++;
				}
			}
		}
		numsub=counter;
		// now parse the addition and subtraction components into numerators and
		// denominators
		boolean[] isdiv=new boolean[10];
		for(int i=0;i<numadd;i++){
			String[] temp2=parsedivide(terms[i][0]);
			terms[i][0]=temp2[0];
			if(temp2.length>1){
				terms[i][5]=temp2[1];
				isdiv[i]=true;
			}
		}
		for(int i=0;i<numsub;i++){
			String[] temp2=parsedivide(terms[i+5][0]);
			terms[i+5][0]=temp2[0];
			if(temp2.length>1){
				terms[i+5][5]=temp2[1];
				isdiv[i+5]=true;
			}
		}
		// finally parse the numerators and denominators into individual
		// coefficients
		nummult=new int[2][10];
		for(int i=0;i<numadd;i++){
			String[] temp2=parsemultiply(terms[i][0]);
			nummult[0][i]=temp2.length;
			if(nummult[0][i]>5){
				nummult[0][i]=5;
			}
			for(int j=0;j<nummult[0][i];j++){
				terms[i][j]=temp2[j];
			}
			if(isdiv[i]){
				temp2=parsemultiply(terms[i][5]);
				nummult[1][i]=temp2.length;
				if(nummult[1][i]>5){
					nummult[1][i]=5;
				}
				for(int j=0;j<nummult[1][i];j++){
					terms[i][j+5]=temp2[j];
				}
			}
		}
		for(int i=0;i<numsub;i++){
			String[] temp2=parsemultiply(terms[i+5][0]);
			nummult[0][i+5]=temp2.length;
			if(nummult[0][i+5]>5){
				nummult[0][i+5]=5;
			}
			for(int j=0;j<nummult[0][i+5];j++){
				terms[i+5][j]=temp2[j];
			}
			if(isdiv[i+5]){
				temp2=parsemultiply(terms[i+5][5]);
				nummult[1][i+5]=temp2.length;
				if(nummult[1][i+5]>5){
					nummult[1][i+5]=5;
				}
				for(int j=0;j<nummult[1][i+5];j++){
					terms[i+5][j+5]=temp2[j];
				}
			}
		}
		// take the numbers and turn them into coefficients
		// the isterm matrix tells the function whether or not to look for a
		// parameter
		// parameter names are converted into numbers corresponding to the
		// parameter list
		coefficients=new double[10][10];
		isterm=new boolean[10][10];
		for(int i=0;i<numadd;i++){
			for(int j=0;j<nummult[0][i];j++){
				for(int k=0;k<paramsnames.length;k++){
					if(terms[i][j].equals(paramsnames[k])){
						coefficients[i][j]=k;
						isterm[i][j]=true;
						break;
					}
				}
				if(!isterm[i][j]){
					coefficients[i][j]=Double.parseDouble(terms[i][j]);
				}
			}
			if(isdiv[i]){
				for(int j=5;j<5+nummult[1][i];j++){
					for(int k=0;k<paramsnames.length;k++){
						if(terms[i][j].equals(paramsnames[k])){
							coefficients[i][j]=k;
							isterm[i][j]=true;
							break;
						}
					}
					if(!isterm[i][j]){
						coefficients[i][j]=Double.parseDouble(terms[i][j]);
					}
				}
			}
		}
		for(int i=5;i<5+numsub;i++){
			for(int j=0;j<nummult[0][i];j++){
				for(int k=0;k<paramsnames.length;k++){
					if(terms[i][j].equals(paramsnames[k])){
						coefficients[i][j]=k;
						isterm[i][j]=true;
						break;
					}
				}
				if(!isterm[i][j]){
					coefficients[i][j]=Double.parseDouble(terms[i][j]);
				}
			}
			if(isdiv[i]){
				for(int j=5;j<5+nummult[1][i];j++){
					for(int k=0;k<paramsnames.length;k++){
						if(terms[i][j].equals(paramsnames[k])){
							coefficients[i][j]=k;
							isterm[i][j]=true;
							break;
						}
					}
					if(!isterm[i][j]){
						coefficients[i][j]=Double.parseDouble(terms[i][j]);
					}
				}
			}
		}
	}

	public String[] parseadd(String formula1){
		String[] temp=new String[5];
		String tempformula=formula1.toString();
		int dumindex=0;
		int counter=0;
		do{
			dumindex=tempformula.indexOf('+');
			if(dumindex>=0){
				temp[counter]=tempformula.substring(0,dumindex);
				counter++;
				tempformula=tempformula.substring(dumindex+1,tempformula.length());
			}else{
				temp[counter]=tempformula.toString();
				counter++;
			}
		}while(dumindex>=0&&counter<5);
		String[] temp2=new String[counter];
		for(int i=0;i<counter;i++){
			temp2[i]=temp[i].toString();
		}
		return temp2;
	}

	public String[] parsesubtract(String formula1){
		String[] temp=new String[5];
		String tempformula=formula1.toString();
		int dumindex=0;
		int counter=0;
		do{
			dumindex=tempformula.indexOf('-');
			if(dumindex>=0){
				temp[counter]=tempformula.substring(0,dumindex);
				counter++;
				tempformula=tempformula.substring(dumindex+1,tempformula.length());
			}else{
				temp[counter]=tempformula.toString();
				counter++;
			}
		}while(dumindex>=0&&counter<5);
		String[] temp2=new String[counter];
		for(int i=0;i<counter;i++){
			temp2[i]=temp[i].toString();
		}
		return temp2;
	}

	public String[] parsemultiply(String formula1){
		String[] temp=new String[5];
		String tempformula=formula1.toString();
		int dumindex=0;
		int counter=0;
		do{
			dumindex=tempformula.indexOf('*');
			if(dumindex>=0){
				temp[counter]=tempformula.substring(0,dumindex);
				counter++;
				tempformula=tempformula.substring(dumindex+1,tempformula.length());
			}else{
				temp[counter]=tempformula.toString();
				counter++;
			}
		}while(dumindex>=0&&counter<5);
		String[] temp2=new String[counter];
		for(int i=0;i<counter;i++){
			temp2[i]=temp[i].toString();
		}
		return temp2;
	}

	public String[] parsedivide(String formula1){
		String[] temp=new String[5];
		String tempformula=formula1.toString();
		int dumindex=0;
		int counter=0;
		do{
			dumindex=tempformula.indexOf('/');
			if(dumindex>=0){
				temp[counter]=tempformula.substring(0,dumindex);
				counter++;
				tempformula=tempformula.substring(dumindex+1,tempformula.length());
			}else{
				temp[counter]=tempformula.toString();
				counter++;
			}
		}while(dumindex>=0&&counter<5);
		String[] temp2=new String[counter];
		for(int i=0;i<counter;i++){
			temp2[i]=temp[i].toString();
		}
		return temp2;
	}

	public double getfunc(double[] params){
		double value=0.0;
		for(int i=0;i<numadd;i++){
			double multvalue=1.0;
			for(int j=0;j<nummult[0][i];j++){
				if(isterm[i][j]){
					multvalue*=params[(int)coefficients[i][j]];
				}else{
					multvalue*=coefficients[i][j];
				}
			}
			double multvalue2=1.0;
			for(int j=5;j<5+nummult[1][i];j++){
				if(isterm[i][j]){
					multvalue2*=params[(int)coefficients[i][j]];
				}else{
					multvalue2*=coefficients[i][j];
				}
			}
			value+=multvalue/multvalue2;
		}
		for(int i=5;i<5+numsub;i++){
			double multvalue=1.0;
			for(int j=0;j<nummult[0][i];j++){
				if(isterm[i][j]){
					multvalue*=params[(int)coefficients[i][j]];
				}else{
					multvalue*=coefficients[i][j];
				}
			}
			double multvalue2=1.0;
			for(int j=5;j<5+nummult[1][i];j++){
				if(isterm[i][j]){
					multvalue2*=params[(int)coefficients[i][j]];
				}else{
					multvalue2*=coefficients[i][j];
				}
			}
			value-=multvalue/multvalue2;
		}
		return value;
	}

}
