/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import jalgs.jdataio;
import jalgs.jsim.rngs;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Image;
import java.awt.Rectangle;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class Plot2DHist{
	// this is a class that draws a plot with a 2D Histogram image in the
	// background
	private float[] xValues;
	private float[] yValues;
	private float[][] histogram;
	private int npts;
	private float xMin,xMax,yMin,yMax,intMin,intMax;
	private float logxmin,logxmax,logymin,logymax;
	private int binSize;
	private int lutindex;
	private int[] lut;
	private int[][] luts;
	private String xLabel,yLabel;
	private boolean logx,logy;
	public static final int histSize=256;
	public static final int WIDTH=256;
	public static final int HEIGHT=256;
	public static final int TICK_LENGTH=3; // length of ticks
	public final Color gridColor=new Color(0xc0c0c0); // light gray
	// public final Color gridColor = new Color(0x000000); //white
	public static final int LEFT_MARGIN=90;
	public static final int RIGHT_MARGIN=18;
	public static final int TOP_MARGIN=20;
	public static final int BOTTOM_MARGIN=50;
	public static final int fontsize=14;
	private float magnification;
	public static final String[] lutnames={"Greys","Red","Green","Blue","NICE","NICE_whiteback","User"};

	public Plot2DHist(String xLabel1,String yLabel1,float[] xValues1,float[] yValues1,int[] lut1){
		initLuts(lut1);
		xValues=xValues1;
		yValues=yValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		lutindex=6;
		lut=luts[lutindex];
		npts=yValues.length;
		logx=false;
		logy=false;
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		binSize=1;
		updateHistogram();
		temp=findminmax(histogram);
		intMin=temp[0];
		intMax=temp[1];
		magnification=1.0f;
	}

	public Plot2DHist(String xLabel1,String yLabel1,float[] hist2d,int width,int height,float startx,float endx,float starty,float endy,int[] lut1){
		// here we get the arrays of x and y points from an existing 2d
		// histogram
		initLuts(lut1);
		npts=0;
		for(int i=0;i<hist2d.length;i++){
			npts+=(int)hist2d[i];
		}
		logx=false;
		logy=false;
		xValues=new float[npts];
		yValues=new float[npts];
		float pwidth=(endx-startx)/(width-1);
		float pheight=(endy-starty)/(height-1);
		int counter=0;
		for(int i=0;i<height;i++){
			float yval=endy-i*pheight;
			for(int j=0;j<width;j++){
				float xval=startx+j*pwidth;
				int pixelpts=(int)hist2d[j+i*width];
				for(int k=0;k<pixelpts;k++){
					xValues[counter]=xval;
					yValues[counter]=yval;
					counter++;
				}
			}
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		lutindex=6;
		lut=luts[lutindex];
		npts=yValues.length;
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		binSize=1;
		updateHistogram();
		magnification=1.0f;
	}

	public Plot2DHist(InputStream is){
		init_from_is(is);
	}

	public Plot2DHist(String filename){
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
		xMin=jdio.readintelfloat(is);
		xMax=jdio.readintelfloat(is);
		yMin=jdio.readintelfloat(is);
		yMax=jdio.readintelfloat(is);
		intMin=jdio.readintelfloat(is);
		intMax=jdio.readintelfloat(is);
		logx=jdio.readintelint(is)==1;
		logy=jdio.readintelint(is)==1;
		lutindex=jdio.readintelint(is);
		binSize=jdio.readintelint(is);
		int[] templut=new int[256];
		jdio.readintelintfile(is,templut);
		initLuts(templut);
		lut=luts[lutindex];
		npts=jdio.readintelint(is);
		xValues=new float[npts];
		jdio.readintelfloatfile(is,npts,xValues);
		yValues=new float[npts];
		jdio.readintelfloatfile(is,npts,yValues);
		updateHistogram();
		magnification=1.0f;
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
		if(temp==4)
			return true;
		else
			return false;
	}
	
	public static boolean is_this(InputStream is) throws IOException{
		jdataio jdio=new jdataio();
		jdio.readstring(is); // read the label
		int temp=jdio.readintelint(is); // now the identifier
		return (temp==4);
	}

	public void initLuts(int[] userlut){
		luts=new int[7][256];
		int[] r={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,
				200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
		int[] g={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,
				144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,127,123,119,115,111,107,103,99,95,91,87,83,79,
				75,71,67,63,59,55,51,47,43,39,35,31,27,23,19,15,11,7,3,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,255};
		int[] b={0,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,127,123,119,115,111,107,103,99,95,
				91,87,83,79,75,71,67,63,59,55,51,47,43,39,35,31,27,23,19,15,11,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,255};
		for(int i=0;i<256;i++){
			luts[0][i]=0xff000000+(i<<16)+(i<<8)+i;
			luts[1][i]=0xff000000+(i<<16);
			luts[2][i]=0xff000000+(i<<8);
			luts[3][i]=0xff000000+i;
			luts[4][i]=0xff000000+(r[i]<<16)+(g[i]<<8)+b[i];
			luts[5][i]=0xff000000+(r[i]<<16)+(g[i]<<8)+b[i];
		}
		luts[5][0]=0xffffffff;
		luts[5][255]=0xffff7878;
		if(userlut!=null){
			luts[6]=userlut;
		}else{
			luts[6]=luts[5];
		}
	}

	public void setLut(int newlut){
		lut=luts[newlut];
		lutindex=newlut;
	}

	public int getLut(){
		return lutindex;
	}

	public void setLutArray(int[] newlut){
		luts[6]=newlut;
		lut=luts[6];
		lutindex=6;
	}

	public int[] getLutArray(){
		return lut;
	}

	public void setBinSize(int newbinsize){
		binSize=newbinsize;
		updateHistogram();
	}

	public int getBinSize(){
		return binSize;
	}

	public void setmagnification(float newmag){
		magnification=newmag;
	}

	public float getmagnification(){
		return magnification;
	}

	private float[] findminmax(float[] arr,int npts1){
		float[] temp=new float[2];
		temp[0]=arr[0];
		temp[1]=arr[0];
		for(int j=0;j<npts1;j++){
			if(arr[j]<temp[0]){
				temp[0]=arr[j];
			}
			if(arr[j]>temp[1]){
				temp[1]=arr[j];
			}
		}
		return temp;
	}

	private float[] findminmax(float[][] arr){
		float[] temp=new float[2];
		temp[0]=arr[0][0];
		temp[1]=arr[0][0];
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<arr[0].length;j++){
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

	private float findmingt0(float[] arr,float max){
		float temp=max;
		if(max<=0.0f){
			return 0.0f;
		}
		for(int j=0;j<arr.length;j++){
			if(arr[j]<temp&&arr[j]>0.0f){
				temp=arr[j];
			}
		}
		return temp;
	}

	private void updateHistogram(){
		int newhistsize=histSize/binSize;
		histogram=new float[newhistsize][newhistsize];
		float tbinsizex=((float)binSize/(float)histSize)*(xMax-xMin);
		float tbinsizey=((float)binSize/(float)histSize)*(yMax-yMin);
		if((!logx)&&(!logy)){
			for(int i=0;i<npts;i++){
				int dumint1=(int)Math.floor((xValues[i]-xMin)/tbinsizex);
				int dumint2=(int)Math.floor((yValues[i]-yMin)/tbinsizey);
				if((dumint1>=0&&dumint1<newhistsize)&&(dumint2>=0&&dumint2<newhistsize)){
					histogram[dumint1][dumint2]+=1.0f;
				}
			}
		}else{
			if(logx){
				if(xMin<=0.0f){
					logxmin=(float)Math.log(findmingt0(xValues,xMax));
				}else{
					logxmin=(float)Math.log(xMin);
				}
				logxmax=(float)Math.log(xMax);
				tbinsizex=((float)binSize/(float)histSize)*(logxmax-logxmin);
			}
			if(logy){
				if(yMin<=0.0f){
					logymin=(float)Math.log(findmingt0(yValues,yMax));
				}else{
					logymin=(float)Math.log(yMin);
				}
				logymax=(float)Math.log(yMax);
				tbinsizey=((float)binSize/(float)histSize)*(logymax-logymin);
			}
			for(int i=0;i<npts;i++){
				int dumint1=(int)Math.floor((xValues[i]-xMin)/tbinsizex);
				int dumint2=(int)Math.floor((yValues[i]-yMin)/tbinsizey);
				if(logx){
					float temp=0.0f;
					if(xValues[i]>0.0f){
						temp=(float)Math.log(xValues[i]);
					}else{
						temp=logxmin;
					}
					dumint1=(int)Math.floor((temp-logxmin)/tbinsizex);
				}
				if(logy){
					float temp=0.0f;
					if(yValues[i]>0.0f){
						temp=(float)Math.log(yValues[i]);
					}else{
						temp=logymin;
					}
					dumint2=(int)Math.floor((temp-logymin)/tbinsizey);
				}
				if(dumint1>=0&&dumint1<newhistsize&&dumint2>=0&&dumint2<newhistsize){
					histogram[dumint1][dumint2]+=1.0f;
				}
			}
		}
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1,double intMin1,double intMax1){
		xMin=(float)xMin1;
		xMax=(float)xMax1;
		yMin=(float)yMin1;
		yMax=(float)yMax1;
		intMin=(float)intMin1;
		intMax=(float)intMax1;
		updateHistogram();
	}

	public void setLimits(float[] limits){
		xMin=limits[0];
		xMax=limits[1];
		yMin=limits[2];
		yMax=limits[3];
		intMin=limits[4];
		intMax=limits[5];
		updateHistogram();
	}

	public void setLogAxes(boolean logx1,boolean logy1){
		logx=logx1;
		logy=logy1;
		updateHistogram();
	}

	public void autoscale(){
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax(histogram);
		intMin=temp[0];
		intMax=temp[1];
		updateHistogram();
	}

	public void xautoscale(){
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		updateHistogram();
	}

	public void yautoscale(){
		float[] temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		updateHistogram();
	}

	public void intautoscale(){
		float[] temp=findminmax(histogram);
		intMin=temp[0];
		intMax=temp[1];
	}

	public void scalerect(Rectangle rect){
		if(rect!=null){
			float tempxmin=((rect.x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymax=((HEIGHT*magnification+TOP_MARGIN*magnification-rect.y)/(HEIGHT*magnification))*(yMax-yMin)+yMin;
			float tempxmax=((rect.x+rect.width-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymin=((HEIGHT*magnification+TOP_MARGIN*magnification-rect.y-rect.height)/(HEIGHT*magnification))*(yMax-yMin)+yMin;
			if(logx){
				float templogxmin=((rect.x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmin=(float)Math.exp(templogxmin);
				float templogxmax=((rect.x+rect.width-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmax=(float)Math.exp(templogxmax);
			}
			if(logy){
				float templogymax=((HEIGHT*magnification+TOP_MARGIN*magnification-rect.y)/(HEIGHT*magnification))*(logymax-logymin)+logymin;
				tempymax=(float)Math.exp(templogymax);
				float templogymin=((HEIGHT*magnification+TOP_MARGIN*magnification-rect.y-rect.height)/(HEIGHT*magnification))*(logymax-logymin)+logymin;
				tempymin=(float)Math.exp(templogymin);
			}
			xMin=tempxmin;
			xMax=tempxmax;
			yMin=tempymin;
			yMax=tempymax;
			updateHistogram();
		}
	}

	public int[] getindices(Roi roi){
		if(roi!=null){
			int[] indices=new int[npts];
			int counter=0;
			for(int i=0;i<npts;i++){
				int[] location=getPlotPosition(xValues[i],yValues[i]);
				if(roi.contains(location[0],location[1])){
					indices[counter]=i;
					counter++;
				}
			}
			if(counter==0)
				return null;
			int[] indices2=new int[counter];
			System.arraycopy(indices,0,indices2,0,counter);
			return indices2;
		}else{
			return null;
		}
	}

	public float get_score(int x,int y){
		return getPlotCoords(x,y)[2];
	}

	public float get_score(Roi roi){
		if(roi!=null){
			int counter=0;
			for(int i=0;i<npts;i++){
				int[] location=getPlotPosition(xValues[i],yValues[i]);
				if(roi.contains(location[0],location[1])){
					counter++;
				}
			}
			return (float)counter/(float)npts;
		}
		return 0.0f;
	}

	public float[] quadrant_analysis(float[] lims){
		float[] scores=new float[4];
		for(int i=0;i<npts;i++){
			if(xValues[i]>lims[0]){
				if(yValues[i]>lims[1]){
					scores[3]+=1.0f;
				}else{
					scores[1]+=1.0f;
				}
			}else{
				if(yValues[i]>lims[1]){
					scores[2]+=1.0f;
				}else{
					scores[0]+=1.0f;
				}
			}
		}
		return scores;
	}

	public float get_pearson(Roi roi){
		return pearson(roi,xValues,yValues);
	}

	public float pearson(Roi roi,float[] xvals,float[] yvals){
		int counter=0;
		double sumxy=0.0;
		double sumx=0.0;
		double sumxx=0.0;
		double sumy=0.0;
		double sumyy=0.0;
		for(int i=0;i<yvals.length;i++){
			if(roi==null){
				counter++;
				sumxy+=xvals[i]*yvals[i];
				sumx+=xvals[i];
				sumxx+=xvals[i]*xvals[i];
				sumy+=yvals[i];
				sumyy+=yvals[i]*yvals[i];
			}else{
				int[] location=getPlotPosition(xvals[i],yvals[i]);
				if(roi.contains(location[0],location[1])){
					counter++;
					sumxy+=xvals[i]*yvals[i];
					sumx+=xvals[i];
					sumxx+=xvals[i]*xvals[i];
					sumy+=yvals[i];
					sumyy+=yvals[i]*yvals[i];
				}
			}
		}
		sumxy/=counter;
		sumx/=counter;
		sumxx/=counter;
		sumy/=counter;
		sumyy/=counter;
		double temp=(double)counter/(double)(counter-1);
		double covar=temp*(sumxy-sumx*sumy);
		double xdev=Math.sqrt(temp*(sumxx-sumx*sumx));
		double ydev=Math.sqrt(temp*(sumyy-sumy*sumy));
		return (float)(covar/(xdev*ydev));
	}

	public float[] pearson_analysis(int ntrials,Roi roi,PlotWindow2DHist parent){
		float[] newvals=new float[ntrials];
		rngs random=new rngs();
		// if(parent!=null) parent.showLog("getting indices");
		int[] indices=getindices(roi);
		if(indices==null){
			indices=new int[npts];
			for(int i=0;i<indices.length;i++)
				indices[i]=i;
		}else{
			/*
			 * for(int i=0;i<indices.length;i++){ if(parent!=null)
			 * parent.showLog(""+i+" , "+indices[i]); }
			 */
		}
		float[] roixvals=new float[indices.length];
		float[] roiyvals=new float[indices.length];
		for(int i=0;i<indices.length;i++){
			roixvals[i]=xValues[indices[i]];
			roiyvals[i]=yValues[indices[i]];
		}
		// if(parent!=null) parent.showLog("starting trials");
		for(int i=0;i<ntrials;i++){
			int[] order=random.random_order(indices.length);
			float[] newyvals=new float[indices.length];
			for(int j=0;j<indices.length;j++){
				newyvals[j]=roiyvals[order[j]];
			}
			newvals[i]=pearson(null,roixvals,newyvals);
			if(parent!=null){
				// parent.showLog("trial"+i);
				parent.showProgress(i,ntrials);
			}
		}
		return newvals;
	}

	public float[] getPlotCoords(int x,int y){
		// this returns an image position in plot coordinates
		float tempx=((x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
		float tempy=((HEIGHT*magnification+TOP_MARGIN*magnification-y)/(HEIGHT*magnification))*(yMax-yMin)+yMin;
		if(logx){
			float templogx=((x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
			tempx=(float)Math.exp(templogx);
		}
		if(logy){
			float templogy=((HEIGHT*magnification+TOP_MARGIN*magnification-y)/(HEIGHT*magnification))*(logymax-logymin)+logymin;
			tempy=(float)Math.exp(templogy);
		}
		float tempcx=(x-LEFT_MARGIN*magnification)/(WIDTH*magnification*binSize);
		float tempcy=(HEIGHT*magnification+TOP_MARGIN*magnification-y)/(HEIGHT*magnification*binSize);
		float density=histogram[(int)tempcx][(int)tempcy]/npts;
		float[] temp={tempx,tempy,density};
		return temp;
	}

	public int[] getPlotPosition(float x,float y){
		// this maps a value to the image
		float posx=(int)((x-xMin)/(xMax-xMin)*histSize);
		float posy=(int)((y-yMin)/(yMax-yMin)*histSize);
		if(logx){
			// assume that logxmin has already been updated
			float temp=0.0f;
			if(x>0.0f)
				temp=(float)Math.log(x);
			else
				temp=logxmin;
			posx=((temp-logxmin)/(logxmax-logxmin)*histSize);
		}
		if(logy){
			// assume that logymin has already been updated
			float temp=0.0f;
			if(y>0.0f)
				temp=(float)Math.log(y);
			else
				temp=logymin;
			posy=((temp-logymin)/(logymax-logymin)*histSize);
		}
		posx*=magnification;
		posy*=magnification;
		posx+=LEFT_MARGIN*magnification;
		posy=(TOP_MARGIN+HEIGHT)*magnification-posy;
		return new int[]{(int)posx,(int)posy};
	}

	public void updateData(float[] xValues1,float[] yValues1,boolean rescale){
		int length=yValues1.length;
		npts=length;
		xValues=xValues1;
		yValues=yValues1;
		if(rescale){
			autoscale();
		}else{
			updateHistogram();
		}
	}

	public int getWidth(){
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newrightmargin=(int)(magnification*RIGHT_MARGIN);
		return newwidth+newleftmargin+newrightmargin;
	}

	public int getHeight(){
		int newtopmargin=(int)(magnification*TOP_MARGIN);
		int newheight=(int)(magnification*HEIGHT);
		int newbottommargin=(int)(magnification*BOTTOM_MARGIN);
		return newheight+newtopmargin+newbottommargin;
	}

	private void drawPlot(Plotter pr){
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(magnification*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(magnification*HEIGHT);
		int newrightmargin=(int)(magnification*RIGHT_MARGIN);
		int newbottommargin=(int)(magnification*BOTTOM_MARGIN);
		int subsize=(int)magnification*binSize;
		int newhistsize=histSize/binSize;

		/*
		 * int width=newwidth+newleftmargin+newrightmargin; int
		 * height=newheight+newtopmargin+newbottommargin; int[] temp=new
		 * int[width*height]; for(int
		 * i=0;i<width*height;i++){temp[i]=0xffffffff;} for(int
		 * i=0;i<newhistsize;i++){ for(int j=0;j<newhistsize;j++){ int
		 * value=(int)((histogram[j][i]-intMin)/(intMax-intMin)*(float)255);
		 * if(value>255){value=255;} if(value<0){value=0;} for(int
		 * k=0;k<subsize;k++){ for(int l=0;l<subsize;l++){
		 * temp[j*subsize+l+newleftmargin
		 * +(newtopmargin+newheight-i*subsize-k)*width]=lut[value]; } } } }
		 * ColorProcessor cp=new ColorProcessor(width,height,temp);
		 * drawAxisLabels(cp); return cp;
		 */

		int[] temp2=new int[newwidth*newheight];
		int fillval=(int)((0.0f-intMin)/(intMax-intMin)*255);
		if(fillval>255)
			fillval=255;
		if(fillval<0)
			fillval=0;
		for(int i=0;i<newwidth*newheight;i++)
			temp2[i]=lut[fillval];
		for(int i=0;i<newhistsize;i++){
			for(int j=0;j<newhistsize;j++){
				int value=(int)((histogram[j][i]-intMin)/(intMax-intMin)*255);
				if(value>255){
					value=255;
				}
				if(value<0){
					value=0;
				}
				for(int k=0;k<subsize;k++){
					int ypos=newheight-i*subsize-k-1;
					if(ypos>=0){
						for(int l=0;l<subsize;l++){
							int xpos=j*subsize+l;
							if(xpos<newwidth)
								temp2[xpos+ypos*newwidth]=lut[value];
						}
					}
				}
			}
		}
		pr.drawImage(new ColorProcessor(newwidth,newheight,temp2).createImage(),newleftmargin,newtopmargin);
		drawAxisLabels(pr);
		pr.setColor(Color.black);
		pr.unsetClip();
		pr.endPlotting();
		// IJ.showMessage("testdraw3");
	}
	
	public FloatProcessor getHistImage(){
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(magnification*HEIGHT);
		int subsize=(int)magnification*binSize;
		int newhistsize=histSize/binSize;

		float[] temp2=new float[newwidth*newheight];
		for(int i=0;i<newhistsize;i++){
			for(int j=0;j<newhistsize;j++){
				for(int k=0;k<subsize;k++){
					int ypos=newheight-i*subsize-k-1;
					if(ypos>=0){
						for(int l=0;l<subsize;l++){
							int xpos=j*subsize+l;
							if(xpos<newwidth)
								temp2[xpos+ypos*newwidth]=histogram[j][i];
						}
					}
				}
			}
		}
		return new FloatProcessor(newwidth,newheight,temp2,null);
	}

	void drawAxisLabels(Plotter pr){
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(magnification*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(magnification*HEIGHT);
		int newbottommargin=(int)(magnification*BOTTOM_MARGIN);
		int newticklength=(int)(magnification*TICK_LENGTH);
		int newfontsize=(int)(magnification*fontsize);

		// calculate the appropriate label numbers
		float[] xticklabels=new float[4];
		float[] yticklabels=new float[4];
		for(int i=0;i<4;i++){
			if(logx){
				float tempx=logxmin+(i/3.0f)*(logxmax-logxmin);
				xticklabels[i]=(float)Math.exp(tempx);
			}else{
				xticklabels[i]=xMin+(i/3.0f)*(xMax-xMin);
			}
			if(logy){
				float tempy=logymin+(i/3.0f)*(logymax-logymin);
				yticklabels[i]=(float)Math.exp(tempy);
			}else{
				yticklabels[i]=yMin+(i/3.0f)*(yMax-yMin);
			}
		}

		pr.setAntiAliasedText(true);
		Font tempfont=new Font("SansSerif",Font.PLAIN,newfontsize);
		pr.setFont(tempfont);
		FontMetrics fm=pr.getFontMetrics();
		int fontheight=fm.getAscent()+fm.getDescent();

		// draw the y axis labels

		String s=jutils.formatted_string(yticklabels[0]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+newheight+fm.getDescent());
		s=jutils.formatted_string(yticklabels[1]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+(2*newheight)/3+fontheight/2);
		s=jutils.formatted_string(yticklabels[2]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+newheight/3+fontheight/2);
		s=jutils.formatted_string(yticklabels[3]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+fontheight/2);

		// draw the x axis labels
		int y=newtopmargin+newheight+fontheight+newticklength+4;
		s=jutils.formatted_string(xticklabels[0]);
		pr.drawString(s,newleftmargin-8,y);
		s=jutils.formatted_string(xticklabels[1]);
		pr.drawString(s,newleftmargin+newwidth/3-pr.getStringWidth(s)/2,y);
		s=jutils.formatted_string(xticklabels[2]);
		pr.drawString(s,newleftmargin+(2*newwidth)/3-pr.getStringWidth(s)/2,y);
		s=jutils.formatted_string(xticklabels[3]);
		pr.drawString(s,newleftmargin+newwidth-pr.getStringWidth(s)+8,y);

		// now draw the axis labels
		pr.drawString(xLabel,newleftmargin+(newwidth-pr.getStringWidth(xLabel))/2,newtopmargin+newheight+fm.getAscent()+newbottommargin/2);
		pr.drawVerticalString(yLabel,newleftmargin-(4*newleftmargin/5),newtopmargin+newheight/2);
		// now draw the tick marks and grid lines
		pr.drawRect(newleftmargin,newtopmargin,newwidth,newheight);
		pr.setColor(gridColor);
		pr.drawLine(newleftmargin+newwidth/3,newtopmargin,newleftmargin+newwidth/3,newtopmargin+newheight);
		pr.drawLine(newleftmargin+(2*newwidth)/3,newtopmargin,newleftmargin+(2*newwidth)/3,newtopmargin+newheight);
		pr.drawLine(newleftmargin,newtopmargin+newheight/3,newleftmargin+newwidth,newtopmargin+newheight/3);
		pr.drawLine(newleftmargin,newtopmargin+(2*newheight)/3,newleftmargin+newwidth,newtopmargin+(2*newheight)/3);
		pr.setColor(Color.black);
		pr.drawLine(newleftmargin,newtopmargin,newleftmargin-newticklength,newtopmargin);
		pr.drawLine(newleftmargin,newtopmargin+newheight/3,newleftmargin-newticklength,newtopmargin+newheight/3);
		pr.drawLine(newleftmargin,newtopmargin+(2*newheight)/3,newleftmargin-newticklength,newtopmargin+(2*newheight)/3);
		pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin-newticklength,newtopmargin+newheight);

		pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin,newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+newwidth/3,newtopmargin+newheight,newleftmargin+newwidth/3,newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+(2*newwidth)/3,newtopmargin+newheight,newleftmargin+(2*newwidth)/3,newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+newwidth,newtopmargin+newheight,newleftmargin+newwidth,newtopmargin+newheight+newticklength);
	}

	Color getColor(int index){
		int temp=index;
		if(temp>=8){
			temp=index%8;
		}
		Color[] temp2={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};
		return temp2[temp];
	}

	void drawPolyline(ImageProcessor ip,int[] x,int[] y,int n){
		ip.moveTo(x[0],y[0]);
		for(int i=0;i<n;i++)
			ip.lineTo(x[i],y[i]);
	}

	public float[] getXValues(){
		return xValues;
	}

	public float[] getYValues(){
		return yValues;
	}

	public String getxLabel(){
		return xLabel;
	}

	public void setxLabel(String xlabel1){
		xLabel=xlabel1;
	}

	public String getyLabel(){
		return yLabel;
	}

	public void setyLabel(String ylabel1){
		yLabel=ylabel1;
	}

	public int getNpts(){
		return npts;
	}

	public float[] getLimits(){
		float[] temp={xMin,xMax,yMin,yMax,intMin,intMax};
		return temp;
	}

	public boolean[] getLogAxes(){
		boolean[] temp={logx,logy};
		return temp;
	}

	public float[][] getHistogram(){
		float[][] temp=new float[histogram.length][histogram[0].length];
		for(int i=0;i<histogram.length;i++){
			for(int j=0;j<histogram[0].length;j++){
				temp[i][j]=histogram[i][j];
			}
		}
		return temp;
	}

	public void saveAsEMF(String path){
		Plotter pr=new EMFPlotter(getWidth(),getHeight(),path);
		drawPlot(pr);
	}
	
	public void saveAsPS(String path){
		Plotter pr=new PSPlotter(getWidth(),getHeight(),path);
		drawPlot(pr);
	}
	
	public void saveAsPDF(String path){
		Plotter pr=new PDFPlotter(getWidth(),getHeight(),path);
		drawPlot(pr);
	}

	public byte[] getEMFBinary(){
		ByteArrayOutputStream os=new ByteArrayOutputStream();
		Plotter pr=new EMFPlotter(getWidth(),getHeight(),os);
		drawPlot(pr);
		return os.toByteArray();
	}

	public ColorProcessor getProcessor(){
		Plotter pr=new CPPlotter(getWidth(),getHeight());
		drawPlot(pr);
		return ((CPPlotter)pr).get_output();
	}

	public Image getImage(){
		return getProcessor().createImage();
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
		jdio.writeintelint(os,4); //specifies the 2D histogram plot type
		jdio.writestring(os,getxLabel()); //x axis label
		jdio.writestring(os,getyLabel()); //y axis label
		jdio.writeintelfloat(os,xMin); // min x axis
		jdio.writeintelfloat(os,xMax); // max x axis
		jdio.writeintelfloat(os,yMin); // min y axis
		jdio.writeintelfloat(os,yMax); // max y axis
		jdio.writeintelfloat(os,intMin); // min z axis
		jdio.writeintelfloat(os,intMax); // max z axis
		jdio.writeintelint(os,logx?1:0); // logx?
		jdio.writeintelint(os,logy?1:0); // logy?
		jdio.writeintelint(os,lutindex); //the selected LUT for visualization
		jdio.writeintelint(os,binSize); //the bin size (pixels)
		jdio.writeintelintarray(os,luts[6]); // the custom lut (256 length integer)
		jdio.writeintelint(os,xValues.length); //the number of data points
		jdio.writeintelfloatarray(os,xValues,xValues.length);
		jdio.writeintelfloatarray(os,yValues,yValues.length);
		// save the errors if they exist
	}

}
