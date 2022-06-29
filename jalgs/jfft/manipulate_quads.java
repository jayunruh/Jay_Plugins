/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;


public class manipulate_quads{

	public float[] shiftxcenter(float[] data,int width,int height){
		float[] temp=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				int dumx=j+width/2;
				if(dumx>=width){
					dumx-=width;
				}
				temp[j+i*width]=data[dumx+width*i];
			}
		}
		return temp;
	}

	public float[] shiftxycenter(float[] data,int width,int height){
		float[] temp=new float[width*height];
		for(int i=0;i<height;i++){
			int dumy=i+height/2;
			if(dumy>=height){
				dumy-=height;
			}
			for(int j=0;j<width;j++){
				int dumx=j+width/2;
				if(dumx>=width){
					dumx-=width;
				}
				temp[j+i*width]=data[dumx+width*dumy];
			}
		}
		return temp;
	}

	public float[] shiftycenter(float[] data,int width,int height){
		float[] temp=new float[width*height];
		for(int i=0;i<height;i++){
			int dumy=i+height/2;
			if(dumy>=height){
				dumy-=height;
			}
			for(int j=0;j<width;j++){
				temp[j+i*width]=data[j+width*dumy];
			}
		}
		return temp;
	}
	
	public static float[] shiftxycenter(float[] data,int xc,int yc,int width,int height){
		int shiftx=xc-width/2;
		int shifty=yc-height/2;
		float[] newdata=new float[width*height];
		for(int i=0;i<height;i++){
			int ypos=i+shifty;
			if(ypos>height) ypos-=height;
			if(ypos<0) ypos+=height;
			for(int j=0;j<width;j++){
				int xpos=j+shiftx;
				if(xpos>width) xpos-=width;
				if(xpos<0) xpos+=width;
				newdata[j+i*width]=data[xpos+ypos*width];
			}
		}
		return newdata;
	}
	
	public static Object[] shiftxyzcenter(Object[] stack,int xc,int yc,int zc,int width,int height){
		int shiftx=xc-width/2;
		int shifty=yc-height/2;
		int shiftz=zc-stack.length/2;
		float[][] newdata=new float[stack.length][width*height];
		for(int i=0;i<stack.length;i++){
			int zpos=i+shiftz;
			if(zpos>=stack.length) zpos-=stack.length;
			if(zpos<0) zpos+=stack.length;
			for(int j=0;j<height;j++){
				int ypos=j+shifty;
				if(ypos>=height) ypos-=height;
				if(ypos<0) ypos+=height;
				int yindex=ypos*width;
				int jindex=j*width;
				for(int k=0;k<width;k++){
					int xpos=k+shiftx;
					if(xpos>=width) xpos-=width;
					if(xpos<0) xpos+=width;
					newdata[i][k+jindex]=((float[])stack[zpos])[xpos+yindex];
				}
			}
		}
		return newdata;
	}

	public float[] avgquadrants(float[] data,int width,int height){
		// note that this only works for nonshifted data
		float[] temp=new float[data.length];
		temp[0]=data[0];
		for(int i=1;i<height/2;i++){
			float avg=(data[i*width]+data[(height-i)*width])/2.0f;
			temp[i*width]=avg;
			temp[(height-i)*width]=avg;
		}
		for(int i=1;i<width/2;i++){
			float avg=(data[i]+data[width-i])/2.0f;
			temp[i]=avg;
			temp[width-i]=avg;
		}
		for(int i=1;i<height/2;i++){
			for(int j=1;j<width/2;j++){
				float avg=(data[i*width+j]+data[i*width+width-j]+data[(height-i)*width+j]+data[(height-i)*width+width-j])/4.0f;
				temp[i*width+j]=avg;
				temp[i*width+width-j]=avg;
				temp[(height-i)*width+j]=avg;
				temp[(height-i)*width+width-j]=avg;
			}
		}
		return temp;
	}

	public float[] duplicatequadrants(float[] data,int width,int height){
		// note that this only works for nonshifted data
		int newwidth=width*2;
		int newheight=height*2;
		float[] temp=new float[newwidth*newheight];
		temp[0]=data[0];
		for(int i=1;i<newheight/2;i++){
			temp[i*newwidth]=data[i*width];
			temp[(newheight-i)*newwidth]=data[i*width];
		}
		for(int i=1;i<newwidth/2;i++){
			temp[i]=data[i];
			temp[newwidth-i]=data[i];
		}
		for(int i=1;i<newheight/2;i++){
			for(int j=1;j<newwidth/2;j++){
				float temp2=data[i*width+j];
				temp[i*newwidth+j]=temp2;
				temp[i*newwidth+newwidth-j]=temp2;
				temp[(newheight-i)*newwidth+j]=temp2;
				temp[(newheight-i)*newwidth+newwidth-j]=temp2;
			}
		}
		return temp;
	}

}
