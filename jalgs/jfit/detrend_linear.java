/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

import jalgs.jstatistics;

public class detrend_linear{
	public linleastsquares lls;
	public int length,segments;
	public int seglength,lengthleft;

	public detrend_linear(int length1,int segments1){
		length=length1;
		segments=segments1;
		seglength=(int)length/segments;
		lengthleft=length-segments*seglength;
		float[][] indvars=new float[2][seglength];
		for(int i=0;i<seglength;i++){
			indvars[0][i]=1.0f;
			indvars[1][i]=(float)i;
		}
		lls=new linleastsquares(indvars);
	}

	public float[] detrend_array(float[] array){
		float avg=jstatistics.getstatistic("Avg",array,null);
		float[] outdata=new float[length];
		float[] coef=new float[2];
		for(int i=0;i<segments;i++){
			float[] subarray=new float[seglength];
			System.arraycopy(array,i*seglength,subarray,0,seglength);
			float[][] temp=detrend_subarray(subarray);
			System.arraycopy(temp[0],0,outdata,i*seglength,seglength);
			coef=temp[1];
		}
		float subval=coef[0]+seglength*coef[1];
		for(int i=0;i<lengthleft;i++){
			outdata[i+seglength*segments]=array[i+seglength*segments]-subval;
		}
		for(int i=0;i<length;i++){
			outdata[i]+=avg;
		}
		return outdata;
	}

	public float[] get_trend(float[] array){
		float[] outdata=new float[length];
		for(int i=0;i<segments;i++){
			float[] subarray=new float[seglength];
			System.arraycopy(array,i*seglength,subarray,0,seglength);
			double[] coef=lls.fitdata(subarray,null);
			for(int j=0;j<seglength;j++){
				outdata[j+i*seglength]=(float)coef[0]+(float)j*(float)coef[1];
			}
		}
		float subval=outdata[(segments-1)*seglength+seglength-1];
		for(int i=0;i<lengthleft;i++){
			outdata[i+seglength*segments]=subval;
		}
		return outdata;
	}

	private float[][] detrend_subarray(float[] subarray){
		double[] coef=lls.fitdata(subarray,null);
		float[] newdata=new float[seglength];
		for(int i=0;i<seglength;i++){
			newdata[i]=subarray[i]-(float)coef[0]-(float)i*(float)coef[1];
		}
		float[][] newdata2=new float[2][];
		newdata2[0]=newdata;
		float[] stats={(float)coef[0],(float)coef[1]};
		newdata2[1]=stats;
		return newdata2;
	}
}
