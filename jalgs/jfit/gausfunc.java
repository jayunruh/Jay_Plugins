/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class gausfunc{

	public float[] model_gaus;

	public gausfunc(){
		model_gaus=new float[1000];
		for(int i=0;i<1000;i++){
			double r=0.01*(double)i;
			model_gaus[i]=(float)Math.exp(-(r*r)/2.0);
		}
	}

	public double getinterpgaus(double r,double stdev){
		if(stdev==0.0){
			if(r==0.0)
				return 1.0f;
			else
				return 0.0f;
		}
		double rrel=r/stdev;
		int rp=(int)(rrel*100.0);
		double rem=rrel*100.0-(double)rp;
		if(rp<999){
			return rem*(model_gaus[rp+1]-model_gaus[rp])+model_gaus[rp];
		}else{
			return 0.0;
		}
	}

	public float[] get_func(double startr,int nr,double dr,double stdev){
		float[] temp=new float[nr];
		for(int i=0;i<nr;i++){
			temp[i]=(float)getinterpgaus(Math.abs(startr+i*dr),stdev);
		}
		return temp;
	}

	public float[] get_norm_func(double startr,int nr,double dr,double stdev){
		float[] temp=get_func(startr,nr,dr,stdev);
		double sum=0.0;
		for(int i=0;i<nr;i++){
			sum+=(double)temp[i];
		}
		for(int i=0;i<nr;i++){
			temp[i]/=sum;
		}
		return temp;
	}

	public float[] get_func(float[] rvals,double stdev){
		float[] temp=new float[rvals.length];
		for(int i=0;i<rvals.length;i++){
			temp[i]=(float)getinterpgaus((double)rvals[i],stdev);
		}
		return temp;
	}

	public float[][] get_2D_func(double startx,int nx,double dx,double starty,int ny,double dy,double stdev){
		float[][] temp=new float[ny][nx];
		for(int i=0;i<ny;i++){
			float yval=(float)getinterpgaus(Math.abs(starty+i*dy),stdev);
			for(int j=0;j<nx;j++){
				temp[i][j]=yval*(float)getinterpgaus(Math.abs(startx+j*dx),stdev);
			}
		}
		return temp;
	}

	public float[][] get_2D_func(double startr,int nr,double dr,double stdev){
		return get_2D_func(startr,nr,dr,startr,nr,dr,stdev);
	}

	public float[][] get_2D_norm_func(double startr,int nr,double dr,double stdev){
		return get_2D_norm_func(startr,nr,dr,startr,nr,dr,stdev);
	}

	public float[][] get_2D_norm_func(double startx,int nx,double dx,double starty,int ny,double dy,double stdev){
		float[][] temp=get_2D_func(startx,nx,dx,starty,ny,dy,stdev);
		double sum=0.0;
		for(int i=0;i<ny;i++){
			for(int j=0;j<nx;j++){
				sum+=(double)temp[i][j];
			}
		}
		for(int i=0;i<ny;i++){
			for(int j=0;j<nx;j++){
				temp[i][j]/=sum;
			}
		}
		return temp;
	}

	public float[] get_2D_func2(double startx,int nx,double dx,double starty,int ny,double dy,double stdev){
		float[] temp=new float[ny*nx];
		for(int i=0;i<ny;i++){
			float yval=(float)getinterpgaus(Math.abs(starty+i*dy),stdev);
			for(int j=0;j<nx;j++){
				temp[j+i*nx]=yval*(float)getinterpgaus(Math.abs(startx+j*dx),stdev);
			}
		}
		return temp;
	}

	public float[] get_2D_func2(double startr,int nr,double dr,double stdev){
		return get_2D_func2(startr,nr,dr,startr,nr,dr,stdev);
	}

	public float[] get_2D_norm_func2(double startr,int nr,double dr,double stdev){
		return get_2D_norm_func2(startr,nr,dr,startr,nr,dr,stdev);
	}

	public float[] get_2D_norm_func2(double startx,int nx,double dx,double starty,int ny,double dy,double stdev){
		float[] temp=get_2D_func2(startx,nx,dx,starty,ny,dy,stdev);
		double sum=0.0;
		for(int i=0;i<nx*ny;i++){
			sum+=(double)temp[i];
		}
		for(int i=0;i<nx*ny;i++){
			temp[i]/=sum;
		}
		return temp;
	}

	public void draw_2D_func(float[] image,double xc,double yc,int width,int height,double stdev,float amp){
		draw_2D_func(image,xc,yc,width,height,stdev,amp,4.0f);
	}

	public void draw_2D_func(float[] image,double xc,double yc,int width,int height,double stdev,float amp,float cutoff){
		int ystart=(int)(yc-cutoff*stdev);
		if(ystart<0)
			ystart=0;
		int yend=(int)(yc+cutoff*stdev);
		if(yend>=height)
			yend=(height-1);
		int xstart=(int)(xc-cutoff*stdev);
		if(xstart<0)
			xstart=0;
		int xend=(int)(xc+cutoff*stdev);
		if(xend>=width)
			xend=(width-1);
		for(int i=ystart;i<=yend;i++){
			float yval=(float)getinterpgaus(Math.abs((double)i-yc),stdev);
			for(int j=xstart;j<=xend;j++){
				float xval=(float)getinterpgaus(Math.abs((double)j-xc),stdev);
				image[j+i*width]+=amp*xval*yval;
			}
		}
	}

	public void draw_3D_func(Object[] stack,double xc,double yc,double zc,int width,int height,double stdevr,double stdevz,float amp){
		draw_3D_func(stack,xc,yc,zc,width,height,stdevr,stdevz,amp,4.0f);
	}

	public void draw_3D_func(Object[] stack,double xc,double yc,double zc,int width,int height,double stdevr,double stdevz,float amp,float cutoff){
		int zstart=(int)(zc-cutoff*stdevz);
		if(zstart<0)
			zstart=0;
		int zend=(int)(zc+cutoff*stdevz);
		if(zend>=stack.length)
			zend=(stack.length-1);
		for(int i=zstart;i<=zend;i++){
			float[] temp=(float[])stack[i];
			float amp2=amp*(float)getinterpgaus(Math.abs((double)i-zc),stdevz);
			draw_2D_func(temp,xc,yc,width,height,stdevr,amp2,cutoff);
		}
	}

}
