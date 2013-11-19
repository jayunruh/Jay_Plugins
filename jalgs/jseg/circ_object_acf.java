/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class circ_object_acf{

	public float[] calc_acf(float[] obj){
		int size=2*obj.length;
		float[] acf=new float[size];
		float max=0.0f;
		for(int d=0;d<size;d++){
			for(int i=0;i<2*size;i++){
				for(int j=0;j<size;j++){
					float r1=(float)Math.sqrt((j-size/2)*(j-size/2)+(i-size/2)*(i-size/2));
					float r2=(float)Math.sqrt((j-size/2)*(j-size/2)+(i-size/2-d)*(i-size/2-d));
					acf[d]+=interpolate(obj,r1)*interpolate(obj,r2);
				}
			}
			if(acf[d]>max){
				max=acf[d];
			}
		}
		if(max>0.0f){
			for(int d=0;d<size;d++){
				acf[d]/=max;
			}
		}
		return acf;
	}

	public float calc_disk_acf(float radius,float distance){
		return (float)(2.0*radius*radius*Math.acos(0.5*distance/radius)-0.5*Math.sqrt(4.0*radius*radius-distance*distance));
	}

	public float[] calc_acf_image(float[] obj){
		float[] racf=calc_acf(obj);
		return radial2image(racf);
	}

	public float[] radial2image(float[] radial){
		int size=radial.length*2;
		float[] image=new float[size*size];
		for(int i=0;i<size;i++){
			for(int j=0;j<size;j++){
				float r=(float)Math.sqrt((j-size/2)*(j-size/2)+(i-size/2)*(i-size/2));
				image[j+i*size]=interpolate(radial,r);
			}
		}
		return image;
	}

	public float interpolate(float[] function,float position){
		if(position==0.0f){
			return function[0];
		}
		if(position==(float)(function.length-1)){
			return function[function.length-1];
		}
		if(position>(float)(function.length-1)){
			return 0.0f;
		}
		float temp=position;
		if(position<0.0f){
			temp=-position;
		}
		int prev=(int)temp;
		float rem=temp-(float)prev;
		return function[prev]+rem*(function[prev+1]-function[prev]);
	}

}
