/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public interface NLLSfitinterface_v3{
	// this interface allows the NLLS fitting routine to gain access to a
	// generic fitting function
	//this version 
	/*
	 * Copyright Jay Unruh Stowers Institute for Medical Research 4/25/08
	 */
	double[] fitfunc(double[] params);
	
	void applyconstraints(double[] params,int[] fixes);
	
	//this function should return the derivatives for each parameter (except the fixed ones) at all positions
	//fit is provided to expediate calculation if possible
	//the first index is the fit coordinate and the second is the parameter index
	double[][] derivfunc(double[] params,int[] fixes,double[] fit); 

	void showresults(String results);

}
