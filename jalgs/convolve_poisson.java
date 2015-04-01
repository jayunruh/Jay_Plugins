/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class convolve_poisson{

	public static double[] cp(double[] data){
		double[] newdata=new double[data.length];
		for(int i=0;i<data.length;i++){
			if(data[i]>0.0){
				for(int j=0;j<data.length;j++){
					newdata[j]+=data[i]*sf.poisson(i,j);
				}
			}
		}
		return newdata;
	}

	public static float[] cp(float[] data){
		double[] newdata=new double[data.length];
		float[] fdata=new float[data.length];
		for(int i=0;i<data.length;i++){
			if(data[i]>0.0){
				for(int j=0;j<data.length;j++){
					newdata[j]+=data[i]*sf.poisson(i,j);
				}
			}
		}
		for(int i=0;i<data.length;i++){
			fdata[i]=(float)newdata[i];
		}
		return fdata;
	}

}
