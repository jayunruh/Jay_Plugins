/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import j3D.*;

import jalgs.jdataio;

import java.awt.*;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import ij.process.*;

public class Plot3D implements EMFexport_interface{
	// private Label coordinates;
	private float[][] xValues;
	private float[][] yValues;
	private float[][][] zValues;
	private int[][] npts;

	protected float xMin,xMax,yMin,yMax,zMin,zMax,xScale,yScale,zScale;
	protected float logxmin,logymin,logxscale,logyscale,logxmax,logymax,logzmin,logzmax,logzscale;
	protected int maxxpts,maxypts,nseries;
	protected int selected;
	protected double rotx,roty,rotz;
	protected int[] shapes,colors;
	protected String xLabel,yLabel,zLabel;
	protected boolean logx,logy,logz;
	public static final int WIDTH=250;
	public static final int HEIGHT=150;
	public static final int TICK_LENGTH=3; // length of ticks
	public final Color gridColor=new Color(0xc0c0c0); // light gray
	public static final int LEFT_MARGIN=125;
	public static final int RIGHT_MARGIN=125;
	public static final int TOP_MARGIN=175;
	public static final int BOTTOM_MARGIN=175;
	public static final int shapesize=5;

	public Plot3D(){
	} // empty constructor for subclasses

	public Plot3D(String xLabel1,String yLabel1,String zLabel1,float[][] xValues1,float[][] yValues1,float[][][] zValues1,Object npts1){
		xValues=xValues1;
		yValues=yValues1;
		zValues=zValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		nseries=zValues.length;
		if(npts1 instanceof int[][]){
			npts=(int[][])npts1;
			maxxpts=npts[0][0];
			maxypts=npts[1][0];
			for(int i=1;i<nseries;i++){
				if(npts[0][i]>maxxpts){
					maxxpts=npts[0][i];
				}
				if(npts[1][i]>maxypts){
					maxypts=npts[1][i];
				}
			}
		}else{
			if(npts1!=null){
				npts=new int[2][nseries];
				npts[0][0]=((int[])npts1)[0];
				npts[1][0]=((int[])npts1)[1];
				for(int i=1;i<nseries;i++){
					npts[0][i]=npts[0][i];
					npts[1][i]=npts[1][i];
				}
				maxxpts=npts[0][0];
				maxypts=npts[1][0];
			}else{
				npts=new int[2][nseries];
				npts[0][0]=zValues[0].length;
				npts[1][0]=zValues[0][0].length;
				maxxpts=npts[0][0];
				maxypts=npts[1][0];
				for(int i=1;i<nseries;i++){
					npts[0][i]=npts[0][0];
					npts[1][i]=npts[1][0];
				}
			}
		}
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
		logx=false;
		logy=false;
		logz=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
		}
		selected=-1;
		rotx=-60.0;
		roty=0.0;
		rotz=-45.0;
	}

	public Plot3D(String xLabel1,String yLabel1,String zLabel1,float[] xValues1,float[] yValues1,float[][] zValues1){
		xValues=new float[1][zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			xValues[0][i]=xValues1[i];
		}
		yValues=new float[1][zValues1[0].length];
		for(int i=0;i<zValues1[0].length;i++){
			yValues[0][i]=yValues1[i];
		}
		zValues=new float[1][zValues1.length][zValues1[0].length];
		for(int i=0;i<zValues1.length;i++){
			for(int j=0;j<zValues1[0].length;j++){
				zValues[0][i][j]=zValues1[i][j];
			}
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		maxxpts=zValues1.length;
		maxypts=zValues1[0].length;
		npts=new int[2][1];
		npts[0][0]=maxxpts;
		npts[1][0]=maxypts;
		nseries=1;
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
		logx=false;
		logy=false;
		logz=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
		}
		selected=-1;
		rotx=-60.0;
		roty=0.0;
		rotz=-45.0;
	}

	public Plot3D(String xLabel1,String yLabel1,String zLabel1,float[][] zValues1,int startxy){
		xValues=new float[1][zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			xValues[0][i]=(float)(i+startxy);
		}
		yValues=new float[1][zValues1[0].length];
		for(int i=0;i<zValues1[0].length;i++){
			yValues[0][i]=(float)(i+startxy);
		}
		zValues=new float[1][zValues1.length][zValues1[0].length];
		for(int i=0;i<zValues1.length;i++){
			for(int j=0;j<zValues1[0].length;j++){
				zValues[0][i][j]=zValues1[i][j];
			}
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		maxxpts=zValues1.length;
		maxypts=zValues1[0].length;
		npts=new int[2][1];
		npts[0][0]=maxxpts;
		npts[1][0]=maxypts;
		nseries=1;
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
		logx=false;
		logy=false;
		logz=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
		}
		selected=-1;
		rotx=-60.0;
		roty=0.0;
		rotz=-45.0;
	}

	public Plot3D(String xLabel1,String yLabel1,String zLabel1,float[][][] zValues1,int startxy,Object npts1){
		zValues=zValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		nseries=zValues.length;
		if(npts1 instanceof int[][]){
			npts=(int[][])npts1;
			maxxpts=npts[0][0];
			maxypts=npts[1][0];
			for(int i=1;i<nseries;i++){
				if(npts[0][i]>maxxpts){
					maxxpts=npts[0][i];
				}
				if(npts[1][i]>maxypts){
					maxypts=npts[1][i];
				}
			}
		}else{
			if(npts1!=null){
				npts=new int[2][nseries];
				npts[0][0]=((int[])npts1)[0];
				npts[1][0]=((int[])npts1)[1];
				for(int i=1;i<nseries;i++){
					npts[0][i]=npts[0][0];
					npts[1][i]=npts[1][0];
				}
				maxxpts=npts[0][0];
				maxypts=npts[1][0];
			}else{
				npts=new int[2][nseries];
				npts[0][0]=zValues[0].length;
				npts[1][0]=zValues[0][0].length;
				maxxpts=npts[0][0];
				maxypts=npts[1][0];
				for(int i=1;i<nseries;i++){
					npts[0][i]=npts[0][0];
					npts[1][i]=npts[1][0];
				}
			}
		}
		xValues=new float[nseries][maxxpts];
		yValues=new float[nseries][maxypts];
		for(int i=0;i<nseries;i++){
			for(int j=0;j<npts[0][i];j++){
				xValues[i][j]=(float)(j+startxy);
			}
			for(int j=0;j<npts[1][i];j++){
				yValues[i][j]=(float)(j+startxy);
			}
		}
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
		logx=false;
		logy=false;
		logz=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
		}
		selected=-1;
		rotx=-60.0;
		roty=0.0;
		rotz=-45.0;
	}

	public Plot3D(InputStream is){
		init_from_is(is);
	}

	public Plot3D(String filename){
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(filename));
			init_from_is(is);
			is.close();
		}catch(IOException e){
			return;
		}
	}

	private void init_from_is(InputStream is){
		jdataio jdio=new jdataio();
		jdio.readstring(is); // read the label
		jdio.readintelint(is); // now the identifier
		xLabel=jdio.readstring(is);
		yLabel=jdio.readstring(is);
		zLabel=jdio.readstring(is);
		nseries=jdio.readintelint(is);
		maxxpts=jdio.readintelint(is);
		maxypts=jdio.readintelint(is);
		npts=new int[2][nseries];
		xValues=new float[nseries][maxxpts];
		yValues=new float[nseries][maxypts];
		zValues=new float[nseries][maxxpts][maxypts];
		shapes=new int[nseries];
		colors=new int[nseries];
		xMin=jdio.readintelfloat(is);
		xMax=jdio.readintelfloat(is);
		yMin=jdio.readintelfloat(is);
		yMax=jdio.readintelfloat(is);
		zMin=jdio.readintelfloat(is);
		zMax=jdio.readintelfloat(is);
		logx=jdio.readintelint(is)==1;
		logy=jdio.readintelint(is)==1;
		logz=jdio.readintelint(is)==1;
		for(int l=0;l<nseries;l++){
			npts[0][l]=jdio.readintelint(is);
			npts[1][l]=jdio.readintelint(is);
			shapes[l]=jdio.readintelint(is);
			colors[l]=jdio.readintelint(is);
			jdio.readintelfloatfile(is,npts[0][l],xValues[l]);
			jdio.readintelfloatfile(is,npts[1][l],yValues[l]);
			for(int i=0;i<npts[0][l];i++)
				jdio.readintelfloatfile(is,npts[1][l],zValues[l][i]); // y
			// values
		}
		selected=-1;
		rotx=-60.0;
		roty=0.0;
		rotz=-45.0;
	}

	public static boolean is_this(String filename){
		int temp=-1;
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(filename));
			jdataio jdio=new jdataio();
			jdio.readstring(is); // read the label
			temp=jdio.readintelint(is); // now the identifier
			is.close();
		}catch(IOException e){
			return false;
		}
		if(temp==1)
			return true;
		else
			return false;
	}

	protected float[] findminmax(float[][] arr,int[] npts1){
		float[] temp=new float[2];
		temp[0]=arr[0][0];
		temp[1]=arr[0][0];
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[i];j++){
				if(arr[i][j]<temp[0]){
					temp[0]=arr[i][j];
				}
				if(arr[i][j]>temp[1]){
					temp[1]=arr[i][j];
				}
			}
		}
		return temp;
	}

	protected float[] findminmax2(float[][][] arr,int[][] npts1){
		float[] temp=new float[2];
		temp[0]=arr[0][0][0];
		temp[1]=arr[0][0][0];
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[0][i];j++){
				for(int k=0;k<npts1[1][i];k++){
					if(arr[i][j][k]<temp[0]){
						temp[0]=arr[i][j][k];
					}
					if(arr[i][j][k]>temp[1]){
						temp[1]=arr[i][j][k];
					}
				}
			}
		}
		return temp;
	}

	protected float findmingt0(float[][] arr,int[] npts1,float max){
		float temp=max;
		if(max<=0.0f){
			return 0.0f;
		}
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[i];j++){
				if(arr[i][j]<temp&&arr[i][j]>0.0f){
					temp=arr[i][j];
				}
			}
		}
		return temp;
	}

	protected float findmingt02(float[][][] arr,int[][] npts1,float max){
		float temp=max;
		if(max<=0.0f){
			return 0.0f;
		}
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[0][i];j++){
				for(int k=0;k<npts[1][i];k++){
					if(arr[i][j][k]<temp&&arr[i][j][k]>0.0f){
						temp=arr[i][j][k];
					}
				}
			}
		}
		return temp;
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1,double zMin1,double zMax1){
		xMin=(float)xMin1;
		xMax=(float)xMax1;
		yMin=(float)yMin1;
		yMax=(float)yMax1;
		zMin=(float)zMin1;
		zMax=(float)zMax1;
	}

	public void setLimits(float[] limits){
		xMin=limits[0];
		xMax=limits[1];
		yMin=limits[2];
		yMax=limits[3];
		zMin=limits[4];
		zMax=limits[5];
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1,boolean logz1){
		logx=logx1;
		logy=logy1;
		logz=logz1;
	}

	public void autoscale(){
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
	}

	public void xautoscale(){
		float[] temp=findminmax(xValues,npts[0]);
		xMin=temp[0];
		xMax=temp[1];
	}

	public void yautoscale(){
		float[] temp=findminmax(yValues,npts[1]);
		yMin=temp[0];
		yMax=temp[1];
	}

	public void zautoscale(){
		float[] temp=findminmax2(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
	}

	public void setrotation(double xrot1,double yrot1,double zrot1){
		rotx=xrot1;
		roty=yrot1;
		rotz=zrot1;
	}

	public void setrotation(double[] rotation){
		rotx=rotation[0];
		roty=rotation[1];
		rotz=rotation[2];
	}

	public double[] getrotation(){
		double[] rotation={rotx,roty,rotz};
		return rotation;
	}

	public void updateSeries(float[] xValues1,float[] yValues1,float[][] zValues1,int series,boolean rescale){
		int xlength=zValues1.length;
		int ylength=zValues1[0].length;
		npts[0][series]=xlength;
		npts[1][series]=ylength;
		if(xlength>maxxpts||ylength>maxypts){
			if(xlength>maxxpts){
				maxxpts=xlength;
			}
			if(ylength>maxypts){
				maxypts=ylength;
			}
			float[][] newxValues=new float[nseries][maxxpts];
			float[][] newyValues=new float[nseries][maxypts];
			float[][][] newzValues=new float[nseries][maxxpts][maxypts];
			for(int i=0;i<series;i++){
				for(int j=0;j<npts[0][i];j++){
					newxValues[i][j]=xValues[i][j];
				}
				for(int j=0;j<npts[1][i];j++){
					newyValues[i][j]=yValues[i][j];
				}
				for(int j=0;j<npts[0][i];j++){
					for(int k=0;k<npts[1][i];k++){
						newzValues[i][j][k]=zValues[i][j][k];
					}
				}
			}
			for(int j=0;j<npts[0][series];j++){
				newxValues[series][j]=xValues1[j];
			}
			for(int j=0;j<npts[1][series];j++){
				newyValues[series][j]=yValues1[j];
			}
			for(int j=0;j<npts[0][series];j++){
				for(int k=0;k<npts[1][series];k++){
					newzValues[series][j][k]=zValues1[j][k];
				}
			}
			for(int i=series+1;i<nseries;i++){
				for(int j=0;j<npts[0][i];j++){
					newxValues[i][j]=xValues[i][j];
				}
				for(int j=0;j<npts[1][i];j++){
					newyValues[i][j]=yValues[i][j];
				}
				for(int j=0;j<npts[0][i];j++){
					for(int k=0;k<npts[1][i];k++){
						newzValues[i][j][k]=zValues[i][j][k];
					}
				}
			}
			xValues=newxValues;
			yValues=newyValues;
			zValues=newzValues;
			if(rescale){
				autoscale();
			}
		}else{
			for(int i=0;i<xlength;i++){
				xValues[series][i]=xValues1[i];
			}
			for(int i=0;i<ylength;i++){
				yValues[series][i]=yValues1[i];
			}
			for(int i=0;i<xlength;i++){
				for(int j=0;j<ylength;j++){
					zValues[series][i][j]=zValues1[i][j];
				}
			}
			if(rescale){
				autoscale();
			}
		}
	}

	public void updateSeries(float[][] zValues1,int series,boolean rescale){
		float[] xValues1=getXValues(series);
		float[] yValues1=getYValues(series);
		updateSeries(xValues1,yValues1,zValues1,series,rescale);
	}

	public void deleteSeries(int series,boolean rescale){
		nseries-=1;
		float[][] newxValues;
		float[][] newyValues;
		float[][][] newzValues;
		int[][] newnpts=new int[2][nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		int newmaxxpts=0;
		int newmaxypts=0;
		if(npts[0][series]==maxxpts){
			for(int i=0;i<=nseries;i++){
				if(i!=series){
					if(npts[0][i]>newmaxxpts){
						newmaxxpts=npts[0][i];
					}
				}
			}
		}else{
			newmaxxpts=maxxpts;
		}
		if(npts[1][series]==maxypts){
			for(int i=0;i<=nseries;i++){
				if(i!=series){
					if(npts[1][i]>newmaxypts){
						newmaxypts=npts[1][i];
					}
				}
			}
		}else{
			newmaxypts=maxypts;
		}
		newxValues=new float[nseries][newmaxxpts];
		newyValues=new float[nseries][newmaxypts];
		newzValues=new float[nseries][newmaxxpts][newmaxypts];
		for(int i=0;i<series;i++){
			newnpts[0][i]=npts[0][i];
			newnpts[1][i]=npts[1][i];
			newshapes[i]=shapes[i];
			newcolors[i]=colors[i];
			for(int j=0;j<newmaxxpts;j++){
				newxValues[i][j]=xValues[i][j];
			}
			for(int j=0;j<newmaxypts;j++){
				newyValues[i][j]=yValues[i][j];
			}
			for(int j=0;j<newmaxxpts;j++){
				for(int k=0;k<newmaxypts;k++){
					newzValues[i][j][k]=zValues[i][j][k];
				}
			}
		}
		for(int i=series+1;i<=nseries;i++){
			newnpts[0][i-1]=npts[0][i];
			newnpts[1][i-1]=npts[1][i];
			newshapes[i-1]=shapes[i];
			newcolors[i-1]=colors[i];
			for(int j=0;j<newmaxxpts;j++){
				newxValues[i-1][j]=xValues[i][j];
			}
			for(int j=0;j<newmaxypts;j++){
				newyValues[i-1][j]=yValues[i][j];
			}
			for(int j=0;j<newmaxxpts;j++){
				for(int k=0;k<newmaxypts;k++){
					newzValues[i-1][j][k]=zValues[i][j][k];
				}
			}
		}
		maxxpts=newmaxxpts;
		maxypts=newmaxypts;
		npts=newnpts;
		xValues=newxValues;
		yValues=newyValues;
		zValues=newzValues;
		shapes=newshapes;
		colors=newcolors;
		if(rescale){
			autoscale();
		}
		// IJ.showMessage("Plot Deleted");
		if(selected>=nseries){
			selected=-1;
		}
		// IJ.showMessage("Selected = "+selected);
	}

	public void addPoints(float[] xValues1,float[] yValues1,float[][] zValues1,boolean rescale){
		nseries++;
		float[][] newxValues;
		float[][] newyValues;
		float[][][] newzValues;
		int[][] newnpts=new int[2][nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		if(yValues1.length>maxypts||xValues1.length>maxxpts){
			newxValues=new float[nseries][xValues1.length];
			newyValues=new float[nseries][yValues1.length];
			newzValues=new float[nseries][xValues1.length][yValues1.length];
			for(int i=0;i<(nseries-1);i++){
				newnpts[0][i]=npts[0][i];
				newnpts[1][i]=npts[1][i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				for(int j=0;j<npts[0][i];j++){
					newxValues[i][j]=xValues[i][j];
				}
				for(int j=0;j<npts[1][i];j++){
					newyValues[i][j]=yValues[i][j];
				}
				for(int j=0;j<npts[0][i];j++){
					for(int k=0;k<npts[1][i];k++){
						newzValues[i][j][k]=zValues[i][j][k];
					}
				}
			}
			maxxpts=xValues1.length;
			maxypts=yValues1.length;
			newnpts[0][nseries-1]=maxxpts;
			newnpts[1][nseries-1]=maxypts;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			for(int j=0;j<maxxpts;j++){
				newxValues[nseries-1][j]=xValues1[j];
			}
			for(int j=0;j<maxypts;j++){
				newyValues[nseries-1][j]=yValues1[j];
			}
			for(int j=0;j<maxxpts;j++){
				for(int k=0;k<maxypts;k++){
					newzValues[nseries-1][j][k]=zValues1[j][k];
				}
			}
		}else{
			newxValues=new float[nseries][maxxpts];
			newyValues=new float[nseries][maxypts];
			newzValues=new float[nseries][maxxpts][maxypts];
			for(int i=0;i<(nseries-1);i++){
				newnpts[0][i]=npts[0][i];
				newnpts[1][i]=npts[1][i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				for(int j=0;j<maxxpts;j++){
					newxValues[i][j]=xValues[i][j];
				}
				for(int j=0;j<maxypts;j++){
					newyValues[i][j]=yValues[i][j];
				}
				for(int j=0;j<maxxpts;j++){
					for(int k=0;k<maxypts;k++){
						newzValues[i][j][k]=zValues[i][j][k];
					}
				}
			}
			newnpts[0][nseries-1]=xValues1.length;
			newnpts[1][nseries-1]=yValues1.length;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			for(int j=0;j<xValues1.length;j++){
				newxValues[nseries-1][j]=xValues1[j];
			}
			for(int j=0;j<yValues1.length;j++){
				newyValues[nseries-1][j]=yValues1[j];
			}
			for(int j=0;j<xValues1.length;j++){
				for(int k=0;k<yValues1.length;k++){
					newzValues[nseries-1][j][k]=zValues1[j][k];
				}
			}
		}
		npts=newnpts;
		shapes=newshapes;
		colors=newcolors;
		xValues=newxValues;
		yValues=newyValues;
		zValues=newzValues;
		if(selected>=nseries){
			selected=-1;
		}
		if(rescale){
			autoscale();
		}
	}

	public void addPoints(float[][] zValues1,boolean rescale,int startxy){
		float[] xValues1=new float[zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			xValues1[i]=(float)(i+startxy);
		}
		float[] yValues1=new float[zValues1[0].length];
		for(int i=0;i<zValues1[0].length;i++){
			yValues1[i]=(float)(i+startxy);
		}
		addPoints(xValues1,yValues1,zValues1,rescale);
	}

	public int getWidth(){
		return WIDTH+LEFT_MARGIN+RIGHT_MARGIN;
	}

	public int getHeight(){
		return HEIGHT+TOP_MARGIN+BOTTOM_MARGIN;
	}

	protected void drawPlot(renderer jr){

		logxmin=0;
		logymin=0;
		logxscale=0;
		logyscale=0;
		logxmax=0;
		logymax=0;
		logzmin=0;
		logzmax=0;
		logzscale=0;
		if(logx){
			if(xMin<=0.0f){
				logxmin=(float)Math.log((double)findmingt0(xValues,npts[0],xMax));
			}else{
				logxmin=(float)Math.log((double)xMin);
			}
			logxmax=(float)Math.log((double)xMax);
			logxscale=(float)WIDTH/(logxmax-logxmin);
		}
		if(logy){
			if(yMin<=0.0f){
				logymin=(float)Math.log((double)findmingt0(yValues,npts[1],yMax));
			}else{
				logymin=(float)Math.log((double)yMin);
			}
			logymax=(float)Math.log((double)yMax);
			logyscale=(float)WIDTH/(logymax-logymin);
		}
		if(logz){
			if(zMin<=0.0f){
				logzmin=(float)Math.log((double)findmingt02(zValues,npts,zMax));
			}else{
				logzmin=(float)Math.log((double)zMin);
			}
			logzmax=(float)Math.log((double)zMax);
			logzscale=(float)HEIGHT/(logzmax-logzmin);
		}
		// IJ.showMessage("testdraw1");

		xScale=(float)WIDTH/(xMax-xMin);
		yScale=(float)WIDTH/(yMax-yMin);
		zScale=(float)HEIGHT/(zMax-zMin);

		drawAxisLabels(jr);
		if(!logx&&!logy&&!logz){
			for(int j=0;j<nseries;j++){
				Color tempcolor=getColor(colors[j]);
				int xpoints[]=new int[npts[0][j]];
				int ypoints[]=new int[npts[1][j]];
				int zpoints[][]=new int[npts[0][j]][npts[1][j]];
				for(int i=0;i<npts[0][j];i++){
					xpoints[i]=LEFT_MARGIN+(int)((xValues[j][i]-xMin)*xScale);
					if(xpoints[i]<LEFT_MARGIN){
						xpoints[i]=LEFT_MARGIN;
					}
					if(xpoints[i]>LEFT_MARGIN+WIDTH){
						xpoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[1][j];i++){
					ypoints[i]=LEFT_MARGIN+(int)((yValues[j][i]-yMin)*yScale);
					if(ypoints[i]<LEFT_MARGIN){
						ypoints[i]=LEFT_MARGIN;
					}
					if(ypoints[i]>LEFT_MARGIN+WIDTH){
						ypoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[0][j];i++){
					for(int k=0;k<npts[1][j];k++){
						zpoints[i][k]=TOP_MARGIN+HEIGHT-(int)((zValues[j][i][k]-zMin)*zScale);
						if(zpoints[i][k]<TOP_MARGIN){
							zpoints[i][k]=TOP_MARGIN;
						}
						if(zpoints[i][k]>TOP_MARGIN+HEIGHT){
							zpoints[i][k]=TOP_MARGIN+HEIGHT;
						}
					}
				}
				if(j!=selected){
					if(shapes[j]==0){
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],tempcolor);
					}else{
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],shapes[j],tempcolor);
					}
				}else{
					if(shapes[j]==0){
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],1,tempcolor);
					}else{
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],tempcolor);
					}
				}
			}
		}else{
			for(int j=0;j<nseries;j++){
				Color tempcolor=getColor(colors[j]);
				int xpoints[]=new int[npts[0][j]];
				int ypoints[]=new int[npts[1][j]];
				int zpoints[][]=new int[npts[0][j]][npts[1][j]];
				for(int i=0;i<npts[0][j];i++){
					if(logx){
						float xtemp;
						if(xValues[j][i]>0.0f){
							xtemp=(float)Math.log((double)xValues[j][i]);
						}else{
							xtemp=logxmin;
						}
						xpoints[i]=LEFT_MARGIN+(int)((xtemp-logxmin)*logxscale);
					}else{
						xpoints[i]=LEFT_MARGIN+(int)((xValues[j][i]-xMin)*xScale);
					}
					if(xpoints[i]<LEFT_MARGIN){
						xpoints[i]=LEFT_MARGIN;
					}
					if(xpoints[i]>LEFT_MARGIN+WIDTH){
						xpoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[1][j];i++){
					if(logy){
						float ytemp;
						if(yValues[j][i]>0.0f){
							ytemp=(float)Math.log((double)yValues[j][i]);
						}else{
							ytemp=logymin;
						}
						ypoints[i]=LEFT_MARGIN+(int)((ytemp-logymin)*logyscale);
					}else{
						ypoints[i]=LEFT_MARGIN+(int)((yValues[j][i]-yMin)*yScale);
					}
					if(ypoints[i]<LEFT_MARGIN){
						ypoints[i]=LEFT_MARGIN;
					}
					if(ypoints[i]>LEFT_MARGIN+WIDTH){
						ypoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[0][j];i++){
					for(int k=0;k<npts[1][j];k++){
						if(logz){
							float ztemp;
							if(zValues[j][i][k]>0.0f){
								ztemp=(float)Math.log((double)zValues[j][i][k]);
							}else{
								ztemp=logzmin;
							}
							zpoints[i][k]=TOP_MARGIN+HEIGHT-(int)((ztemp-logzmin)*logzscale);
						}else{
							zpoints[i][k]=TOP_MARGIN+HEIGHT-(int)((zValues[j][i][k]-zMin)*zScale);
						}
						if(zpoints[i][k]<TOP_MARGIN){
							zpoints[i][k]=TOP_MARGIN;
						}
						if(zpoints[i][k]>TOP_MARGIN+HEIGHT){
							zpoints[i][k]=TOP_MARGIN+HEIGHT;
						}
					}
				}
				if(j!=selected){
					if(shapes[j]==0){
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],tempcolor);
					}else{
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],shapes[j],tempcolor);
					}
				}else{
					if(shapes[j]==0){
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],1,tempcolor);
					}else{
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[0][j],npts[1][j],tempcolor);
					}
				}
			}
		}
		jr.setrotation((int)rotx,(int)roty,(int)rotz);
	}

	private void drawAxisLabels(renderer jr){

		// calculate the appropriate label numbers
		float[] xticklabels=new float[4];
		float[] yticklabels=new float[4];
		float[] zticklabels=new float[4];
		for(int i=0;i<4;i++){
			if(logx){
				float tempx=logxmin+((float)i/3.0f)*(logxmax-logxmin);
				xticklabels[i]=(float)Math.exp((double)tempx);
			}else{
				xticklabels[i]=xMin+((float)i/3.0f)*(xMax-xMin);
			}
			if(logy){
				float tempy=logymin+((float)i/3.0f)*(logymax-logymin);
				yticklabels[i]=(float)Math.exp((double)tempy);
			}else{
				yticklabels[i]=yMin+((float)i/3.0f)*(yMax-yMin);
			}
			if(logz){
				float tempz=logzmin+((float)i/3.0f)*(logzmax-logzmin);
				zticklabels[i]=(float)Math.exp((double)tempz);
			}else{
				zticklabels[i]=zMin+((float)i/3.0f)*(zMax-zMin);
			}
		}

		// draw the z axis labels
		String s=jutils.formatted_string((double)zticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+HEIGHT-5,Color.BLACK);
		s=jutils.formatted_string((double)zticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+(int)((2*HEIGHT)/3)-5,Color.BLACK);
		s=jutils.formatted_string((double)zticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+(int)(HEIGHT/3)-5,Color.BLACK);
		s=jutils.formatted_string((double)zticklabels[3]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN-5,Color.BLACK);

		jr.addText3D(zLabel,LEFT_MARGIN+WIDTH+25,LEFT_MARGIN,TOP_MARGIN+HEIGHT/2-15,Color.BLACK);

		// now the x axis labels
		s=jutils.formatted_string((double)xticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)xticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+(int)(WIDTH/3)-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)xticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+(int)(2*WIDTH/3)-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)xticklabels[3]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addText3D(xLabel,LEFT_MARGIN+WIDTH/2,LEFT_MARGIN+WIDTH+60,TOP_MARGIN+HEIGHT,Color.BLACK);

		// and the y axis labels
		s=jutils.formatted_string((double)yticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)yticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+(int)(WIDTH/3)+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)yticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+(int)(2*WIDTH/3)+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string((double)yticklabels[3]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+WIDTH+10,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addText3D(yLabel,LEFT_MARGIN+WIDTH+60,LEFT_MARGIN+WIDTH/2,TOP_MARGIN+HEIGHT,Color.BLACK);

		// finally the grid lines
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+HEIGHT,Color.BLACK);
		jr.addLine3D(LEFT_MARGIN+WIDTH/3,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN+WIDTH/3,LEFT_MARGIN,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN+(2*WIDTH)/3,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN+(2*WIDTH)/3,LEFT_MARGIN,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN,Color.BLACK);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+HEIGHT/3,LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN+HEIGHT/3,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+(2*HEIGHT)/3,LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN+(2*HEIGHT)/3,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+WIDTH/3,TOP_MARGIN,LEFT_MARGIN,LEFT_MARGIN+WIDTH/3,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+(2*WIDTH)/3,TOP_MARGIN,LEFT_MARGIN,LEFT_MARGIN+(2*WIDTH)/3,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN,LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN,LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN,Color.BLACK);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+HEIGHT/3,LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT/3,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+(2*HEIGHT)/3,LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN+(2*HEIGHT)/3,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN,TOP_MARGIN+HEIGHT,LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addLine3D(LEFT_MARGIN+WIDTH/3,LEFT_MARGIN,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH/3,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN+(2*WIDTH)/3,LEFT_MARGIN,TOP_MARGIN+HEIGHT,LEFT_MARGIN+(2*WIDTH)/3,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN+WIDTH,LEFT_MARGIN,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+WIDTH/3,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH,LEFT_MARGIN+WIDTH/3,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+(2*WIDTH)/3,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH,LEFT_MARGIN+(2*WIDTH)/3,TOP_MARGIN+HEIGHT,gridColor);
		jr.addLine3D(LEFT_MARGIN,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,LEFT_MARGIN+WIDTH,LEFT_MARGIN+WIDTH,TOP_MARGIN+HEIGHT,Color.BLACK);
	}

	private void drawPolyline(renderer jr,int[] xpoints,int[] ypoints,int[][] zpoints,int nxpts,int nypts,Color color){
		for(int i=0;i<nypts;i++){
			for(int j=1;j<nxpts;j++){
				jr.addLine3D(xpoints[j-1],ypoints[i],zpoints[j-1][i],xpoints[j],ypoints[i],zpoints[j][i],color);
			}
		}
		for(int i=0;i<nxpts;i++){
			for(int j=1;j<nypts;j++){
				jr.addLine3D(xpoints[i],ypoints[j-1],zpoints[i][j-1],xpoints[i],ypoints[j],zpoints[i][j],color);
			}
		}
	}

	private void drawPolyshape(renderer jr,int[] xpoints,int[] ypoints,int[][] zpoints,int nxpts,int nypts,int shape,Color color){
		for(int i=0;i<nypts;i++){
			for(int j=0;j<nxpts;j++){
				// jr.addPoint3D(xpoints[j],ypoints[i],zpoints[j][i],shapesize,color,Point3D.CIRCLE);
				jr.addPoint3D(xpoints[j],ypoints[i],zpoints[j][i],shape-1,color);
			}
		}
	}

	Color getColor(int index){
		int temp=index;
		if(temp>=8){
			temp=index%8;
		}
		Color[] temp2={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};
		return temp2[temp];
	}

	public void selectSeries(int series){
		selected=series;
		if(selected>=nseries){
			selected=-1;
		}
	}

	public int getSelected(){
		return selected;
	}

	public float[][] getXValues(){
		return xValues;
	}

	public float[] getXValues(int series){
		return xValues[series];
	}

	public float[][] getYValues(){
		return yValues;
	}

	public float[] getYValues(int series){
		return yValues[series];
	}

	public float[][][] getZValues(){
		return zValues;
	}

	public float[][] getZValues(int series){
		return zValues[series];
	}

	public String getxLabel(){
		return xLabel;
	}

	public void setxLabel(String xLabel1){
		xLabel=xLabel1;
	}

	public String getyLabel(){
		return yLabel;
	}

	public void setyLabel(String yLabel1){
		yLabel=yLabel1;
	}

	public String getzLabel(){
		return zLabel;
	}

	public void setzLabel(String zLabel1){
		zLabel=zLabel1;
	}

	public int[][] getNpts(){
		return npts;
	}

	public int getNSeries(){
		return nseries;
	}

	public int getmaxxpts(){
		return maxxpts;
	}

	public int getmaxypts(){
		return maxypts;
	}

	public float[] getLimits(){
		float[] temp={xMin,xMax,yMin,yMax,zMin,zMax};
		return temp;
	}

	public boolean[] getLogAxes(){
		boolean[] temp={logx,logy,logz};
		return temp;
	}

	public int[] getShapes(){
		return shapes;
	}

	public int[] getColors(){
		return colors;
	}

	public ColorProcessor getProcessor(){
		return new ColorProcessor(getImage());
	}

	public Image getImage(){
		renderer jr=new renderer(getWidth(),getHeight());
		drawPlot(jr);
		return jr.renderimage();
	}

	public void saveAsEMF(String path){
		byte[] binaryEMF=getEMFBinary();
		jdataio jdio=new jdataio();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			jdio.writebytearray(os,binaryEMF);
			os.close();
		}catch(IOException e){
			return;
		}
	}
	
	public void saveAsPS(String path){
		byte[] binaryPS=getPSBinary();
		jdataio jdio=new jdataio();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			jdio.writebytearray(os,binaryPS);
			os.close();
		}catch(IOException e){
			return;
		}
	}
	
	public byte[] getPSBinary(){
		renderer jr=new renderer(getWidth(),getHeight());
		drawPlot(jr);
		return jr.renderPS();
	}

	public byte[] getEMFBinary(){
		renderer jr=new renderer(getWidth(),getHeight());
		drawPlot(jr);
		return jr.renderEMF();
	}

	public void saveplot2file(String filename){
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(filename));
			saveplot2os(os);
			os.close();
		}catch(IOException e){
			return;
		}
	}

	public void saveplot2os(OutputStream os){
		jdataio jdio=new jdataio();
		// start with unique identifier for a 3D plot
		jdio.writestring(os,"pw2_file_type");
		jdio.writeintelint(os,1);
		jdio.writestring(os,getxLabel());
		jdio.writestring(os,getyLabel());
		jdio.writestring(os,getzLabel());
		jdio.writeintelint(os,nseries); // number of series'
		jdio.writeintelint(os,getmaxxpts()); // max number of pts
		jdio.writeintelint(os,getmaxypts()); // max number of pts
		jdio.writeintelfloat(os,xMin); // min x axis
		jdio.writeintelfloat(os,xMax); // max x axis
		jdio.writeintelfloat(os,yMin); // min y axis
		jdio.writeintelfloat(os,yMax); // max y axis
		jdio.writeintelfloat(os,zMin); // min z axis
		jdio.writeintelfloat(os,zMax); // max z axis
		jdio.writeintelint(os,logx?1:0); // logx?
		jdio.writeintelint(os,logy?1:0); // logy?
		jdio.writeintelint(os,logz?1:0); // logz?
		for(int l=0;l<nseries;l++){
			jdio.writeintelint(os,npts[0][l]); // number of points in this
			// series
			jdio.writeintelint(os,npts[1][l]); // number of points in this
			// series
			jdio.writeintelint(os,shapes[l]); // shape index
			jdio.writeintelint(os,colors[l]); // color index
			jdio.writeintelfloatarray(os,xValues[l],npts[0][l]); // x values
			jdio.writeintelfloatarray(os,yValues[l],npts[1][l]); // y values
			for(int i=0;i<npts[0][l];i++)
				jdio.writeintelfloatarray(os,zValues[l][i],npts[1][l]); // y
			// values
		}
		// save the errors if they exist
	}

}
