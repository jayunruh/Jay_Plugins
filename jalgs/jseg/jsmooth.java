/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.interpolation;
import jalgs.jfit.gausfunc;

public class jsmooth{
	// here we have simple multidimensional boxcar smoothing and binning

	public static float[] smooth1D(float[] data){
		float[] retdata=new float[data.length];
		for(int i=1;i<(data.length-1);i++){
			retdata[i]=(data[i-1]+data[i]+data[i+1])/3.0f;
		}
		retdata[0]=(data[0]+data[1])/2.0f;
		retdata[data.length-1]=(data[data.length-2]+data[data.length-1])/2.0f;
		return retdata;
	}

	public static float[] smooth2D(float[] data,int width,int height){
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				retdata[j+i*width]=algutils.getNeighbors2Avg(data,j,i,width,height);
			}
		}
		return retdata;
	}
	
	public static float[] smooth2Dnot0(float[] data,int width,int height){
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(data[j+i*width]!=0.0f) retdata[j+i*width]=algutils.getNeighbors2Stat(data,j,i,width,height,"Not0Avg");
			}
		}
		return retdata;
	}
	
	public static float[] smooth2Dnot0(float[] data,int width,int height,int size){
		//here size is the size of the smoothing kernel
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(data[j+i*width]!=0.0f) retdata[j+i*width]=algutils.get_region_stat("Not0Avg",data,j,i,size,size,width,height);
			}
		}
		return retdata;
	}
	
	public static float[] smooth2D(float[] data,int width,int height,int size,String stat){
		//here size is the size of the smoothing kernel
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				retdata[j+i*width]=algutils.get_region_stat(stat,data,j,i,size,size,width,height);
			}
		}
		return retdata;
	}
	
	public static float[] smooth2DCircle(Object data,int width,int height,int size,String stat){
		//here size is the size of the smoothing kernel
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				retdata[j+i*width]=algutils.get_circle_stat(stat,data,j,i,size,width,height);
			}
		}
		return retdata;
	}
	
	public static float[] smooth2Dobjects(float[] data,int width,int height,int size,float[] objects){
		//here size is the size of the smoothing kernel
		float[] retdata=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(data[j+i*width]!=0.0f){
					float id=objects[j+i*width];
					float[] datareg=algutils.get_region(data,j,i,size,size,width,height);
					float[] objreg=algutils.get_region(objects,j,i,size,size,width,height);
					double avg=0.0; int cnt=0;
					for(int k=0;k<objreg.length;k++){
						if(objreg[k]==id){avg+=datareg[k]; cnt++;}
					}
					retdata[j+i*width]=(float)(avg/cnt);
				}
			}
		}
		return retdata;
	}

	public static Object[] smooth3D(Object[] data,int width,int height){
		Object[] retdata=new Object[data.length];
		for(int z=0;z<data.length;z++){
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					retdata[j+i*width]=algutils.getNeighbors2Avg(data,j,i,z,width,height);
				}
			}
		}
		return retdata;
	}

	public static Object bin2D(Object data,int width,int height,int binby,boolean avg){
		float[] temp=algutils.convert_arr_float2(data);
		float[] temp2=bin2D(temp,width,height,binby,avg);
		if(data instanceof float[])
			return temp2;
		else if(data instanceof short[])
			return algutils.convert_arr_short(temp2);
		else if(data instanceof int[])
			return algutils.convert_arr_int(temp2);
		else
			return algutils.convert_arr_byte(temp2);
	}

	public static float[] bin2D(float[] data,int width,int height,int binby,boolean avg){
		int newwidth=(int)((float)width/(float)binby);
		int newheight=(int)((float)height/(float)binby);
		float[] newdata=new float[newwidth*newheight];
		for(int i=0;i<newheight;i++){
			for(int j=0;j<binby;j++){
				for(int k=0;k<newwidth;k++){
					for(int l=0;l<binby;l++){
						newdata[k+i*newwidth]+=data[l+k*binby+(j+i*binby)*width];
					}
				}
			}
		}
		float temp=binby*binby;
		if(avg){
			for(int i=0;i<newheight*newwidth;i++){
				newdata[i]/=temp;
			}
		}
		return newdata;
	}
	
	public static Object enlarge2D(Object data,int width,int height,int enlargeby,boolean interp){
		float[] temp=algutils.convert_arr_float2(data);
		float[] temp2=enlarge2D(temp,width,height,enlargeby,interp);
		if(data instanceof float[])
			return temp2;
		else if(data instanceof short[])
			return algutils.convert_arr_short(temp2);
		else if(data instanceof int[])
			return algutils.convert_arr_int(temp2);
		else
			return algutils.convert_arr_byte(temp2);
	}
	
	public static float[] enlarge2D(float[] data,int width,int height,int enlargeby,boolean interp){
		int newwidth=width*enlargeby;
		int newheight=height*enlargeby;
		float[] newdata=new float[newwidth*newheight];
		if(!interp){
    		for(int i=0;i<newheight;i++){
    			int oldy=(int)((float)i/(float)enlargeby);
    			for(int j=0;j<newwidth;j++){
    				int oldx=(int)((float)j/(float)enlargeby);
    				newdata[j+i*newwidth]=data[oldx+oldy*width];
    			}
    		}
		} else {
    		for(int i=0;i<newheight;i++){
    			float oldy=(float)i/(float)enlargeby;
    			for(int j=0;j<newwidth;j++){
    				float oldx=(float)j/(float)enlargeby;
    				newdata[j+i*newwidth]=interpolation.interp2D(data,width,height,oldx,oldy);
    			}
    		}
		}
		return newdata;
	}

	public static void blur1D(float[] data,float stdev){
		gausfunc gf=new gausfunc();
		blur1D(data,stdev,gf);
	}

	public static void blur1D(float[] data,float stdev,gausfunc gf){
		float[] blurfunc=generate_symblurfunc(stdev,gf);
		convsym1D(data,blurfunc);
	}

	public static float[] generate_symblurfunc(float stdev,gausfunc gf){
		// this generates a symmetric gaussian blurring function
		int funclength=(int)(3.0f*stdev);
		if(funclength<2)
			funclength=2;
		float[] blurfunc=gf.get_func(0,funclength,1.0,stdev);
		// normalize as though the mirror image is present
		float sum=blurfunc[0];
		for(int i=1;i<blurfunc.length;i++){
			sum+=2.0f*blurfunc[i];
		}
		for(int i=0;i<blurfunc.length;i++)
			blurfunc[i]/=sum;
		return blurfunc;
	}

	public static void convsym1D(float[] data,float[] symblurfunc){
		// here we convolve a function with a 1D symetrical function
		// assume that the data set is flat outside the viewing region
		float[] newdata=new float[data.length];
		for(int i=0;i<data.length;i++){
			newdata[i]+=symblurfunc[0]*data[i];
			for(int j=1;j<symblurfunc.length;j++){
				int pos=i+j;
				if(pos>=data.length)
					pos=data.length-1;
				newdata[i]+=symblurfunc[j]*data[pos];
				pos=i-j;
				if(pos<0)
					pos=0;
				newdata[i]+=symblurfunc[j]*data[pos];
			}
		}
		System.arraycopy(newdata,0,data,0,data.length);
	}

	public static void blur2D(float[] data,float stdev,int width,int height){
		gausfunc gf=new gausfunc();
		blur2D(data,stdev,width,height,gf);
	}

	public static void blur2D(float[] data,float stdev,int width,int height,gausfunc gf){
		float[] symblurfunc=generate_symblurfunc(stdev,gf);
		float[] newdata=new float[width*height];
		// first blur in the x direction while copying into newdata
		for(int k=0;k<height;k++){
			int offset=k*width;
			for(int i=0;i<width;i++){
				newdata[offset+i]+=symblurfunc[0]*data[offset+i];
				for(int j=1;j<symblurfunc.length;j++){
					int pos=i+j;
					if(pos>=width)
						pos=width-1;
					newdata[offset+i]+=symblurfunc[j]*data[offset+pos];
					pos=i-j;
					if(pos<0)
						pos=0;
					newdata[offset+i]+=symblurfunc[j]*data[offset+pos];
				}
			}
		}
		// now copy back into the original array, bluring in the y direction
		for(int i=0;i<width*height;i++)
			data[i]=0.0f;
		for(int k=0;k<width;k++){
			for(int i=0;i<height;i++){
				int offset=i*width;
				data[k+offset]+=symblurfunc[0]*newdata[k+offset];
				for(int j=1;j<symblurfunc.length;j++){
					int pos=i+j;
					if(pos>=height)
						pos=height-1;
					data[k+offset]+=symblurfunc[j]*newdata[k+pos*width];
					pos=i-j;
					if(pos<0)
						pos=0;
					data[k+offset]+=symblurfunc[j]*data[k+pos*width];
				}
			}
		}
	}

	public static void blur3D(Object[] data,float stdev,float stdevz,int width,int height,gausfunc gf){
		// here it is assumed that the slices are floating point
		for(int i=0;i<data.length;i++)
			blur2D((float[])data[i],stdev,width,height,gf);
		if(data.length>1){
			float[] symblurfunc=generate_symblurfunc(stdevz,gf);
			for(int i=0;i<width*height;i++){
				float[] zprofile=(float[])algutils.get_stack_col(data,width,height,i,data.length);
				convsym1D(zprofile,symblurfunc);
				algutils.set_stack_col(data,width,height,i,data.length,zprofile);
			}
		}
	}

	public static void mexhat(float[] data,float stdev1,float stdev2,int width,int height,gausfunc gf){
		float[] data1=data.clone();
		blur2D(data1,stdev1,width,height,gf);
		float[] data2=data.clone();
		blur2D(data2,stdev2,width,height,gf);
		// now subtract the larger stdev (the second) from the smaller (the
		// first)
		for(int i=0;i<data.length;i++)
			data[i]=data1[i]-data2[i];
	}

	public static void mexhat3D(Object[] data,float stdev,float stdevz,float stdevratio,int width,int height,gausfunc gf){
		mexhat3D(data,stdev,stdevz,stdevratio,width,height,gf,false);
	}

	public static void mexhat3D(Object[] data,float stdev,float stdevz,float stdevratio,int width,int height,gausfunc gf,boolean trunc){
		Object[] data1=algutils.clone_obj_array(data);
		blur3D(data1,stdev,stdevz,width,height,gf);
		Object[] data2=algutils.clone_obj_array(data);
		blur3D(data2,stdev*stdevratio,stdevz*stdevratio,width,height,gf);
		// now subtract the larger stdev (the second) from the smaller (the
		// first)
		for(int i=0;i<data.length;i++){
			float[] temp1=(float[])data1[i];
			float[] temp2=(float[])data2[i];
			float[] temp3=(float[])data[i];
			for(int j=0;j<width*height;j++){
				temp3[j]=temp1[j]-temp2[j];
				if(trunc&&temp3[j]<0.0f)
					temp3[j]=0.0f;
			}
		}
	}

}
