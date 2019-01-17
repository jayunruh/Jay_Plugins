/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import jalgs.jdataio;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.image.IndexColorModel;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DecimalFormat;

public class Plot4{
	// this is a class that draws a plot usually associated with a PlotWindow4
	// window
	protected float[][] xValues;
	protected float[][] yValues;
	protected float[][][] errors;
	protected boolean showerrors;
	protected int[] npts;
	protected float xMin,xMax,yMin,yMax,xScale,yScale;
	protected float logxmin,logymin,logxscale,logyscale,logxmax,logymax;
	protected int maxpts,nseries;
	protected int selected;
	protected int[] shapes,colors,shapesizes;
	protected String xLabel,yLabel;
	protected boolean logx,logy;
	protected String[] annotations;
	public static final int WIDTH=400;
	public static final int HEIGHT=200;
	public static final int TICK_LENGTH=3; // length of ticks
	public Color gridColor=new Color(0xc0c0c0); // light gray
	public static final int LEFT_MARGIN=90;
	public static final int RIGHT_MARGIN=18;
	public static final int TOP_MARGIN=20;
	public static final int BOTTOM_MARGIN=50;
	public int shapesize=8;
	public int fontsize=14;
	protected float magnification,magratio;
	public static final String[] color_names={"black","blue","green","red","magenta","cyan","yellow","orange"};
	public static final Color[] java_colors={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};
	public static final String[] shape_names={"line","square","+","x","triangle","column"};

	public Plot4(){
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
		nseries=1;
		showerrors=false;
	}

	public Plot4(String xLabel1,String yLabel1,float[][] xValues1,float[][] yValues1,Object npts1){
		xValues=xValues1;
		yValues=yValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		nseries=yValues.length;
		if(npts1 instanceof int[]){
			npts=(int[])npts1;
			maxpts=npts[0];
			for(int i=1;i<nseries;i++){
				if(npts[i]>maxpts){
					maxpts=npts[i];
				}
			}
		}else{
			if(npts1!=null){
				npts=new int[nseries];
				npts[0]=((Integer)npts1).intValue();
				for(int i=1;i<nseries;i++)
					npts[i]=npts[0];
				maxpts=npts[0];
			}else{
				npts=new int[nseries];
				npts[0]=yValues[0].length;
				maxpts=npts[0];
				for(int i=1;i<nseries;i++)
					npts[i]=npts[0];
			}
		}
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		shapesizes=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
			shapesizes[i]=shapesize;
		}
		selected=-1;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		errors=null;
		showerrors=false;
	}

	public Plot4(String xLabel1,String yLabel1,float[][] yValues1,Object npts1){
		yValues=yValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		nseries=yValues.length;
		if(npts1 instanceof int[]){
			npts=(int[])npts1;
			maxpts=npts[0];
			for(int i=1;i<nseries;i++){
				if(npts[i]>maxpts){
					maxpts=npts[i];
				}
			}
		}else{
			if(npts1!=null){
				npts=new int[nseries];
				npts[0]=((Integer)npts1).intValue();
				for(int i=1;i<nseries;i++)
					npts[i]=npts[0];
				maxpts=npts[0];
			}else{
				npts=new int[nseries];
				npts[0]=yValues[0].length;
				maxpts=npts[0];
				for(int i=1;i<nseries;i++)
					npts[i]=npts[0];
			}
		}
		xValues=new float[nseries][maxpts];
		for(int i=0;i<nseries;i++){
			for(int j=0;j<npts[i];j++){
				xValues[i][j]=(float)j+1;
			}
		}
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		shapesizes=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
			shapesizes[i]=shapesize;
		}
		selected=-1;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		errors=null;
		showerrors=false;
	}

	public Plot4(String xLabel1,String yLabel1,float[] xValues1,float[] yValues1){
		xValues=new float[1][yValues1.length];
		for(int i=0;i<yValues1.length;i++){
			xValues[0][i]=xValues1[i];
		}
		yValues=new float[1][yValues1.length];
		for(int i=0;i<yValues1.length;i++){
			yValues[0][i]=yValues1[i];
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		maxpts=yValues[0].length;
		npts=new int[1];
		npts[0]=maxpts;
		nseries=yValues.length;
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		shapesizes=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
			shapesizes[i]=shapesize;
		}
		selected=-1;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		errors=null;
		showerrors=false;
	}

	public Plot4(String xLabel1,String yLabel1,float[] yValues1){
		xValues=new float[1][yValues1.length];
		for(int i=0;i<yValues1.length;i++){
			xValues[0][i]=i+1;
		}
		yValues=new float[1][yValues1.length];
		for(int i=0;i<yValues1.length;i++){
			yValues[0][i]=yValues1[i];
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		maxpts=yValues[0].length;
		npts=new int[1];
		npts[0]=maxpts;
		nseries=yValues.length;
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		shapes=new int[nseries];
		colors=new int[nseries];
		shapesizes=new int[nseries];
		for(int i=0;i<nseries;i++){
			colors[i]=i;
			shapesizes[i]=shapesize;
		}
		selected=-1;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		errors=null;
		showerrors=false;
	}

	public Plot4(InputStream is){
		init_from_is(is);
	}

	public Plot4(String filename){
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
		xLabel=jdio.readstring(is);
		if(xLabel.equals("pw2_file_type")){
			int id=jdio.readintelint(is);
			xLabel=jdio.readstring(is);
		}
		yLabel=jdio.readstring(is);
		nseries=jdio.readintelint(is);
		maxpts=jdio.readintelint(is);
		npts=new int[nseries];
		xValues=new float[nseries][maxpts];
		yValues=new float[nseries][maxpts];
		shapes=new int[nseries];
		colors=new int[nseries];
		shapesizes=new int[nseries];
		xMin=jdio.readintelfloat(is);
		xMax=jdio.readintelfloat(is);
		yMin=jdio.readintelfloat(is);
		yMax=jdio.readintelfloat(is);
		logx=jdio.readintelint(is)==1;
		logy=jdio.readintelint(is)==1;
		for(int l=0;l<nseries;l++){
			shapesizes[l]=shapesize;
			npts[l]=jdio.readintelint(is);
			shapes[l]=jdio.readintelint(is);
			colors[l]=jdio.readintelint(is);
			jdio.readintelfloatfile(is,npts[l],xValues[l]);
			jdio.readintelfloatfile(is,npts[l],yValues[l]);
		}
		showerrors=jdio.readintelint(is)==1;
		if(showerrors){
			errors=new float[2][nseries][maxpts];
			for(int l=0;l<nseries;l++){
				jdio.readintelfloatfile(is,npts[l],errors[0][l]);
				jdio.readintelfloatfile(is,npts[l],errors[1][l]);
			}
		}
		boolean annotated=jdio.readintelint(is)==1;
		if(annotated){
			annotations=new String[nseries];
			for(int i=0;i<nseries;i++){
				annotations[i]=jdio.readstring(is);
			}
		}
		selected=-1;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
	}

	public static boolean is_this(String filename){
		int temp=-1;
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(filename));
			jdataio jdio=new jdataio();
			String label=jdio.readstring(is); // read the label
			if(label.equals("pw2_file_type"))
				temp=jdio.readintelint(is); // now the identifier
			else if(filename.endsWith("pw"))
				temp=0; // if no label exists, it may be a pw file
			is.close();
		}catch(IOException e){
			return false;
		}
		if(temp==0)
			return true;
		else
			return false;
	}
	
	public static boolean is_this(InputStream is) throws IOException{
		//here we assume that this is some sort of plot family
		jdataio jdio=new jdataio();
		String label=jdio.readstring(is); // read the label
		if(label.equals("pw2_file_type")){
			int temp=jdio.readintelint(is); // now the identifier
			if(temp==0) return true;
		}
		return false;
	}

	public void setmagnification(float newmag){
		magnification=newmag;
	}

	public float getmagnification(){
		return magnification;
	}

	public void setmagratio(float newmag){
		magratio=newmag;
	}

	public float getmagratio(){
		return magratio;
	}

	public void setFontSize(int fontsize){
		this.fontsize=fontsize;
	}

	public void setGridWhiteness(int whiteness){
		int val=whiteness;
		if(whiteness>255)
			val=255;
		if(whiteness<0)
			val=0;
		gridColor=new Color(val,val,val);
	}

	public void setShapeSize(int shapesize){
		this.shapesize=shapesize;
	}
	
	public void setAnnotations(String[] annotations){
		this.annotations=annotations;
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

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		xMin=(float)xMin1;
		xMax=(float)xMax1;
		yMin=(float)yMin1;
		yMax=(float)yMax1;
	}

	public void setLimits(float[] limits){
		xMin=limits[0];
		xMax=limits[1];
		yMin=limits[2];
		yMax=limits[3];
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1){
		logx=logx1;
		logy=logy1;
	}

	public void setShowErrors(boolean showerrors){
		this.showerrors=showerrors;
	}

	public void setBinSize(float newbinsize){
	}

	public float getBinSize(){
		return 0.0f;
	}

	public float getBinSizeUnits(){
		return 0.0f;
	}

	public void setBinSizeUnits(float newbinsize){
	}

	public int[] getrectindices(Rectangle rect){
		return null;
	}

	public void updateSeries(float[] xValues1,boolean rescale){
	}

	public void updateSeries(float[] hist,float startx,float endx,boolean rescale){
	}

	public float[][] getHistogram(){
		return null;
	}

	public void autoscale(){
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
	}

	public void xautoscale(){
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
	}

	public void yautoscale(){
		float[] temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
	}

	public void scalerect(Rectangle rect){
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		if(rect!=null){
			float tempxmin=((rect.x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymax=((HEIGHT*ymag+TOP_MARGIN*ymag-rect.y)/(HEIGHT*ymag))*(yMax-yMin)+yMin;
			float tempxmax=((rect.x+rect.width-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymin=((HEIGHT*ymag+TOP_MARGIN*ymag-rect.y-rect.height)/(HEIGHT*ymag))*(yMax-yMin)+yMin;
			if(logx){
				float templogxmin=((rect.x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmin=(float)Math.exp(templogxmin);
				float templogxmax=((rect.x+rect.width-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmax=(float)Math.exp(templogxmax);
			}
			if(logy){
				float templogymax=((HEIGHT*ymag+TOP_MARGIN*ymag-rect.y)/(HEIGHT*ymag))*(logymax-logymin)+logymin;
				tempymax=(float)Math.exp(templogymax);
				float templogymin=((HEIGHT*ymag+TOP_MARGIN*ymag-rect.y-rect.height)/(HEIGHT*ymag))*(logymax-logymin)+logymin;
				tempymin=(float)Math.exp(templogymin);
			}
			xMin=tempxmin;
			xMax=tempxmax;
			yMin=tempymin;
			yMax=tempymax;
		}
	}

	public float[] getPlotCoords(int x,int y){
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		float tempxmin=((x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(xMax-xMin)+xMin;
		float tempymax=((HEIGHT*ymag+TOP_MARGIN*ymag-y)/(HEIGHT*ymag))*(yMax-yMin)+yMin;
		if(logx){
			float templogxmin=((x-LEFT_MARGIN*magnification)/(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
			tempxmin=(float)Math.exp(templogxmin);
		}
		if(logy){
			float templogymax=((HEIGHT*ymag+TOP_MARGIN*ymag-y)/(HEIGHT*ymag))*(logymax-logymin)+logymin;
			tempymax=(float)Math.exp(templogymax);
		}
		float[] temp={tempxmin,tempymax};
		return temp;
	}
	
	/*****************
	 * this resets values like logxmin, logymin without drawing the plot (normally where these would be set)
	 */
	public void resetInternalScales() {
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(ymag*HEIGHT);
		//we can assume that the logy limits and all scales have been updated recently
		logxmin=0;
		logymin=0;
		logxscale=0;
		logyscale=0;
		logxmax=0;
		logymax=0;
		if(logx){
			if(xMin<=0.0f){
				logxmin=(float)Math.log(findmingt0(xValues,npts,xMax));
			}else{
				logxmin=(float)Math.log(xMin);
			}
			logxmax=(float)Math.log(xMax);
			logxscale=newwidth/(logxmax-logxmin);
		}
		if(logy){
			if(yMin<=0.0f){
				logymin=(float)Math.log(findmingt0(yValues,npts,yMax));
			}else{
				logymin=(float)Math.log(yMin);
			}
			logymax=(float)Math.log(yMax);
			logyscale=newheight/(logymax-logymin);
		}
		xScale=newwidth/(xMax-xMin);
		yScale=newheight/(yMax-yMin);
	}
	
	/*************
	 * this returns the pixel position on the plot image for the indicated coordinates, must draw plot (or resetInternalScales before doing this
	 * @param xpos
	 * @param ypos
	 * @return
	 */
	public int[] getCoordsPosition(float xpos,float ypos) {
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newheight=(int)(ymag*HEIGHT);
		int xpos2=newleftmargin+(int)((xpos-xMin)*xScale);
		int ypos2=newtopmargin+newheight-(int)((ypos-yMin)*yScale);
		if(logx) {
			float xtemp;
			if(xpos>0.0f) xtemp=(float)Math.log(xpos);
			else xtemp=logxmin;
			xpos2=newleftmargin+(int)((xtemp-logxmin)*logxscale);
		}
		if(logy) {
			float ytemp;
			if(ypos>0.0f) ytemp=(float)Math.log(ypos);
			else ytemp=logymin;
			ypos2=newtopmargin+newheight-(int)((ytemp-logymin)*logyscale);
		}
		return new int[]{xpos2,ypos2};
	}

	public Calibration getCalibration(){
		Calibration cal=new Calibration();
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(ymag*HEIGHT);
		float[] origin=getPlotCoords(newleftmargin,newtopmargin+newheight);
		cal.xOrigin=origin[0];
		cal.yOrigin=origin[1];
		cal.pixelWidth=1.0/xScale;
		cal.pixelHeight=1.0/yScale;
		cal.setInvertY(true);
		return cal;
	}

	public void updateSeries(float[] xValues1,float[] yValues1,int series,boolean rescale){
		int length=yValues1.length;
		npts[series]=length;
		if(length>maxpts){
			float[][] newxValues=new float[nseries][length];
			float[][] newyValues=new float[nseries][length];
			for(int i=0;i<series;i++){
				for(int j=0;j<maxpts;j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
				}
			}
			for(int j=0;j<length;j++){
				newxValues[series][j]=xValues1[j];
				newyValues[series][j]=yValues1[j];
			}
			for(int i=series+1;i<nseries;i++){
				for(int j=0;j<maxpts;j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
				}
			}
			maxpts=length;
			xValues=newxValues;
			yValues=newyValues;
			if(rescale){
				autoscale();
			}
		}else{
			for(int i=0;i<length;i++){
				xValues[series][i]=xValues1[i];
				yValues[series][i]=yValues1[i];
			}
			for(int i=length;i<maxpts;i++){
				xValues[series][i]=0.0f;
				yValues[series][i]=0.0f;
			}
			if(rescale){
				autoscale();
			}
		}
	}

	public void updateSeries(float[] yValues1,int series,boolean rescale){
		float[] xValues1=getXValues(series);
		updateSeries(xValues1,yValues1,series,rescale);
	}

	public void deleteSeries(int series,boolean rescale){
		nseries-=1;
		float[][] newxValues;
		float[][] newyValues;
		int[] newnpts=new int[nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		int newmaxpts=0;
		if(npts[series]==maxpts){
			for(int i=0;i<=nseries;i++){
				if(i!=series){
					if(npts[i]>newmaxpts){
						newmaxpts=npts[i];
					}
				}
			}
		}else{
			newmaxpts=maxpts;
		}
		newxValues=new float[nseries][newmaxpts];
		newyValues=new float[nseries][newmaxpts];
		for(int i=0;i<series;i++){
			newnpts[i]=npts[i];
			newshapes[i]=shapes[i];
			newcolors[i]=colors[i];
			for(int j=0;j<newmaxpts;j++){
				newxValues[i][j]=xValues[i][j];
				newyValues[i][j]=yValues[i][j];
			}
		}
		for(int i=series+1;i<=nseries;i++){
			newnpts[i-1]=npts[i];
			newshapes[i-1]=shapes[i];
			newcolors[i-1]=colors[i];
			for(int j=0;j<newmaxpts;j++){
				newxValues[i-1][j]=xValues[i][j];
				newyValues[i-1][j]=yValues[i][j];
			}
		}
		maxpts=newmaxpts;
		npts=newnpts;
		xValues=newxValues;
		yValues=newyValues;
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

	public void deleteMultiSeries(boolean[] fate,boolean rescale){
		int ndelete=0;
		for(int i=0;i<nseries;i++){
			if(fate[i])
				ndelete++;
		}
		// nseries-=1;
		float[][] newxValues;
		float[][] newyValues;
		int[] newnpts=new int[nseries-ndelete];
		int[] newshapes=new int[nseries-ndelete];
		int[] newcolors=new int[nseries-ndelete];
		int newmaxpts=0;
		for(int i=0;i<nseries;i++){
			if(!fate[i]){
				if(npts[i]>newmaxpts){
					newmaxpts=npts[i];
				}
			}
		}
		newxValues=new float[nseries-ndelete][newmaxpts];
		newyValues=new float[nseries-ndelete][newmaxpts];
		int counter=0;
		for(int i=0;i<nseries;i++){
			if(!fate[i]){
				newnpts[counter]=npts[i];
				newshapes[counter]=shapes[i];
				newcolors[counter]=colors[i];
				System.arraycopy(xValues[i],0,newxValues[counter],0,newmaxpts);
				System.arraycopy(yValues[i],0,newyValues[counter],0,newmaxpts);
				counter++;
			}
		}
		String[] newannot=null;
		if(annotations!=null){
			newannot=new String[nseries-ndelete];
			counter=0;
			for(int i=0;i<nseries;i++){
				if(!fate[i]){newannot[counter]=annotations[i]; counter++;}
			}
		}
		nseries-=ndelete;
		maxpts=newmaxpts;
		npts=newnpts;
		xValues=newxValues;
		yValues=newyValues;
		shapes=newshapes;
		colors=newcolors;
		annotations=newannot;
		if(rescale){
			autoscale();
		}
		selected=-1;
	}

	public void addPoints(float[] xValues1,float[] yValues1,boolean rescale){
		nseries++;
		float[][] newxValues;
		float[][] newyValues;
		int[] newnpts=new int[nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		int[] newshapesizes=new int[nseries];
		if(yValues1.length>maxpts){
			newxValues=new float[nseries][yValues1.length];
			newyValues=new float[nseries][yValues1.length];
			for(int i=0;i<(nseries-1);i++){
				newnpts[i]=npts[i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				newshapesizes[i]=shapesizes[i];
				for(int j=0;j<maxpts;j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
				}
			}
			maxpts=yValues1.length;
			newnpts[nseries-1]=maxpts;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			newshapesizes[nseries-1]=shapesize;
			for(int j=0;j<maxpts;j++){
				newxValues[nseries-1][j]=xValues1[j];
				newyValues[nseries-1][j]=yValues1[j];
			}
		}else{
			newxValues=new float[nseries][maxpts];
			newyValues=new float[nseries][maxpts];
			for(int i=0;i<(nseries-1);i++){
				newnpts[i]=npts[i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				newshapesizes[i]=shapesizes[i];
				for(int j=0;j<maxpts;j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
				}
			}
			newnpts[nseries-1]=yValues1.length;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			newshapesizes[nseries-1]=shapesize;
			for(int j=0;j<yValues1.length;j++){
				newxValues[nseries-1][j]=xValues1[j];
				newyValues[nseries-1][j]=yValues1[j];
			}
		}
		npts=newnpts;
		shapes=newshapes;
		colors=newcolors;
		shapesizes=newshapesizes;
		xValues=newxValues;
		yValues=newyValues;
		if(selected>=nseries){
			selected=-1;
		}
		if(errors!=null){
			float[][][] newerrs=new float[2][nseries][maxpts];
			for(int i=0;i<errors[0].length;i++){
				System.arraycopy(errors[0][i],0,newerrs[0][i],0,errors[0][i].length);
				System.arraycopy(errors[1][i],0,newerrs[1][i],0,errors[1][i].length);
			}
			errors=newerrs;
		}
		if(rescale){
			autoscale();
		}
	}

	public void addPoints(float[] yValues1,boolean rescale){
		float[] xValues1=new float[yValues1.length];
		for(int i=0;i<yValues1.length;i++){
			xValues1[i]=i+1;
		}
		addPoints(xValues1,yValues1,rescale);
	}

	public void addErrors(float[][][] errors){
		this.errors=errors;
	}

	public void addErrors(float[][] errors){
		// here the upper and lower error bars are the same
		this.errors=new float[2][nseries][];
		this.errors[0]=errors;
		this.errors[1]=errors;
	}
	
	public void addSeriesErrors(int series,float[][] errors){
		if(this.errors==null){
			this.errors=new float[2][nseries][maxpts];
		}
		System.arraycopy(errors[0],0,this.errors[0][series],0,npts[series]);
		System.arraycopy(errors[1],0,this.errors[1][series],0,npts[series]);
	}
	
	public void addSeriesErrors(int series,float[] errors){
		addSeriesErrors(series,new float[][]{errors,errors});
	}

	public int getHeight(){
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newheight=(int)(ymag*HEIGHT);
		int newbottommargin=(int)(ymag*BOTTOM_MARGIN);
		return newheight+newtopmargin+newbottommargin;
	}

	public int getWidth(){
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newrightmargin=(int)(magnification*RIGHT_MARGIN);
		return newwidth+newleftmargin+newrightmargin;
	}

	private void drawPlot(Plotter pr){
		// this draws the plot on the current plotter
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(ymag*HEIGHT);
		Rectangle frame=new Rectangle(newleftmargin,newtopmargin,newwidth,newheight);
		// pr.clear_plot();
		// int[] temp=new int[width*height];
		// for(int i=0;i<width*height;i++){temp[i]=0xffffffff;}
		// ColorProcessor cp=new ColorProcessor(width,height,temp);

		logxmin=0;
		logymin=0;
		logxscale=0;
		logyscale=0;
		logxmax=0;
		logymax=0;
		if(logx){
			if(xMin<=0.0f){
				logxmin=(float)Math.log(findmingt0(xValues,npts,xMax));
			}else{
				logxmin=(float)Math.log(xMin);
			}
			logxmax=(float)Math.log(xMax);
			logxscale=newwidth/(logxmax-logxmin);
		}
		if(logy){
			if(yMin<=0.0f){
				logymin=(float)Math.log(findmingt0(yValues,npts,yMax));
			}else{
				logymin=(float)Math.log(yMin);
			}
			logymax=(float)Math.log(yMax);
			logyscale=newheight/(logymax-logymin);
		}
		// IJ.showMessage("testdraw1");

		drawAxisLabels(pr);
		pr.setClipRect(frame.x,frame.y,frame.width,frame.height);
		xScale=newwidth/(xMax-xMin);
		yScale=newheight/(yMax-yMin);

		// IJ.showMessage("testdraw2");
		if(!logx&&!logy){
			for(int j=0;j<nseries;j++){
				pr.shapesize=(int)(shapesizes[j]*magnification);
				pr.setColor(getColor(colors[j]));
				int xpoints[]=new int[npts[j]];
				int ypoints[]=new int[npts[j]];
				for(int i=0;i<npts[j];i++){
					xpoints[i]=newleftmargin+(int)((xValues[j][i]-xMin)*xScale);
					ypoints[i]=newtopmargin+frame.height-(int)((yValues[j][i]-yMin)*yScale);
				}
				if(j!=selected){
					if(shapes[j]==0){
						pr.drawPolyline(xpoints,ypoints,npts[j]);
					}else if(shapes[j]==5){
						pr.drawPolyColumns(xpoints,ypoints,newtopmargin+frame.height,npts[j],Color.black);
					}else{
						pr.drawPolyshape(xpoints,ypoints,shapes[j],npts[j]);
					}
				}else{
					if(shapes[j]==0){
						pr.drawPolyshape(xpoints,ypoints,1,npts[j]);
					}else{
						pr.drawPolyline(xpoints,ypoints,npts[j]);
					}
				}
				if(showerrors&&errors!=null){
					int[] yerrptsu=new int[npts[j]];
					int[] yerrptsl=new int[npts[j]];
					for(int i=0;i<npts[j];i++){
						yerrptsu[i]=newtopmargin+frame.height-(int)((yValues[j][i]+errors[1][j][i]-yMin)*yScale);
						yerrptsl[i]=newtopmargin+frame.height-(int)((yValues[j][i]-errors[0][j][i]-yMin)*yScale);
					}
					pr.drawPolyerrors(xpoints,yerrptsu,yerrptsl,npts[j]);
				}
			}
		}else{
			for(int j=0;j<nseries;j++){
				pr.shapesize=(int)(shapesizes[j]*magnification);
				pr.setColor(getColor(colors[j]));
				int xpoints[]=new int[npts[j]];
				int ypoints[]=new int[npts[j]];
				for(int i=0;i<npts[j];i++){
					if(logx){
						float xtemp;
						if(xValues[j][i]>0.0f){
							xtemp=(float)Math.log(xValues[j][i]);
						}else{
							xtemp=logxmin;
						}
						xpoints[i]=newleftmargin+(int)((xtemp-logxmin)*logxscale);
					}else{
						xpoints[i]=newleftmargin+(int)((xValues[j][i]-xMin)*xScale);
					}
				}
				for(int i=0;i<npts[j];i++){
					if(logy){
						float ytemp;
						if(yValues[j][i]>0.0f){
							ytemp=(float)Math.log(yValues[j][i]);
						}else{
							ytemp=logymin;
						}
						ypoints[i]=newtopmargin+frame.height-(int)((ytemp-logymin)*logyscale);
					}else{
						ypoints[i]=newtopmargin+frame.height-(int)((yValues[j][i]-yMin)*yScale);
					}
				}
				if(j!=selected){
					if(shapes[j]==0){
						pr.drawPolyline(xpoints,ypoints,npts[j]);
					}else if(shapes[j]==5){
						pr.drawPolyColumns(xpoints,ypoints,newtopmargin+frame.height,npts[j],Color.black);
					}else{
						pr.drawPolyshape(xpoints,ypoints,shapes[j],npts[j]);
					}
				}else{
					if(shapes[j]==0){
						pr.drawPolyshape(xpoints,ypoints,1,npts[j]);
					}else{
						pr.drawPolyline(xpoints,ypoints,npts[j]);
					}
				}
				if(showerrors&&errors!=null){
					int[] yerrptsu=new int[npts[j]];
					int[] yerrptsl=new int[npts[j]];
					if(logy){
						for(int i=0;i<npts[j];i++){
							float ytemp=yValues[j][i]+errors[1][j][i];
							if(ytemp>0.0f){
								ytemp=(float)Math.log(ytemp);
							}else{
								ytemp=logymin;
							}
							yerrptsu[i]=newtopmargin+frame.height-(int)((ytemp-logymin)*logyscale);
							ytemp=yValues[j][i]-errors[1][j][i];
							if(ytemp>0.0f){
								ytemp=(float)Math.log(ytemp);
							}else{
								ytemp=logymin;
							}
							yerrptsl[i]=newtopmargin+frame.height-(int)((ytemp-logymin)*yScale);
						}
					}else{
						for(int i=0;i<npts[j];i++){
							yerrptsu[i]=newtopmargin+frame.height-(int)((yValues[j][i]+errors[1][j][i]-yMin)*yScale);
							yerrptsl[i]=newtopmargin+frame.height-(int)((yValues[j][i]-errors[0][j][i]-yMin)*yScale);
						}
					}
					pr.drawPolyerrors(xpoints,yerrptsu,yerrptsl,npts[j]);
				}
			}
		}
		pr.setColor(Color.black);
		pr.unsetClip();
		pr.endPlotting();
		// IJ.showMessage("testdraw3");
	}

	void drawAxisLabels(Plotter pr){
		// set up formating
		DecimalFormat expformat=new DecimalFormat("0.00E0");
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(ymag*HEIGHT);
		int newrightmargin=(int)(magnification*RIGHT_MARGIN);
		int newbottommargin=(int)(ymag*BOTTOM_MARGIN);
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
		if(gridColor.getRed()<=254){
			pr.setColor(gridColor);
			pr.drawLine(newleftmargin+newwidth/3,newtopmargin,newleftmargin+newwidth/3,newtopmargin+newheight);
			pr.drawLine(newleftmargin+(2*newwidth)/3,newtopmargin,newleftmargin+(2*newwidth)/3,newtopmargin+newheight);
			pr.drawLine(newleftmargin,newtopmargin+newheight/3,newleftmargin+newwidth,newtopmargin+newheight/3);
			pr.drawLine(newleftmargin,newtopmargin+(2*newheight)/3,newleftmargin+newwidth,newtopmargin+(2*newheight)/3);
		}
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
		return java_colors[temp];
	}

	public void selectSeries(int series){
		selected=series;
		if(selected>=nseries){
			selected=-1;
		}
		if(selected<-1){
			selected=nseries-1;
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

	public float[][][] getErrors(){
		return errors;
	}

	public float[] getErrors(int series,boolean upper){
		if(errors==null){
			return null;
		}
		if(upper){
			return errors[1][series];
		}else{
			return errors[0][series];
		}
	}

	public boolean getShowErrors(){
		return showerrors;
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

	public int[] getNpts(){
		return npts;
	}

	public int getNSeries(){
		return nseries;
	}

	public int getmaxpts(){
		return maxpts;
	}

	public float[] getLimits(){
		float[] temp={xMin,xMax,yMin,yMax};
		return temp;
	}

	public boolean[] getLogAxes(){
		boolean[] temp={logx,logy};
		return temp;
	}

	public int[] getShapes(){
		return shapes;
	}

	public int[] getColors(){
		return colors;
	}

	public int[] getShapeSizes(){
		return shapesizes;
	}
	
	public String[] getAnnotations(){
		return annotations;
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
	
	public byte[] getPSBinary(){
		ByteArrayOutputStream os=new ByteArrayOutputStream();
		Plotter pr=new PSPlotter(getWidth(),getHeight(),os);
		drawPlot(pr);
		return os.toByteArray();
	}
	
	public byte[] getPDFBinary(){
		ByteArrayOutputStream os=new ByteArrayOutputStream();
		Plotter pr=new PDFPlotter(getWidth(),getHeight(),os);
		drawPlot(pr);
		return os.toByteArray();
	}
	
	public byte[] getSVGBinary(){
		ByteArrayOutputStream os=new ByteArrayOutputStream();
		Plotter pr=new SVGPlotter(getWidth(),getHeight(),os);
		drawPlot(pr);
		return os.toByteArray();
	}

	public ColorProcessor getProcessor(){
		Plotter pr=new CPPlotter(getWidth(),getHeight());
		drawPlot(pr);
		return ((CPPlotter)pr).get_output();
	}

	public ImageProcessor get8Processor(){
		//some colors are swapping here for some reason.
		ColorProcessor cp=getProcessor();
		if(cp==null)
			return null;
		int[] pixels=(int[])cp.getPixels();
		byte[] pixels8=new byte[pixels.length];
		//lut colors are 0white, 1black, 2blue, 3green, 4red, 5magenta, 6cyan, 7yellow, 8orange
		byte[] rLUT=new byte[256];
		rLUT[0]=(byte)255;
		rLUT[4]=(byte)255;
		rLUT[5]=(byte)255;
		rLUT[7]=(byte)255;
		rLUT[8]=(byte)255;
		byte[] gLUT=new byte[256];
		gLUT[0]=(byte)255;
		gLUT[3]=(byte)255;
		gLUT[6]=(byte)255;
		gLUT[7]=(byte)255;
		gLUT[8]=(byte)200;
		byte[] bLUT=new byte[256];
		bLUT[0]=(byte)255;
		bLUT[2]=(byte)255;
		bLUT[5]=(byte)255;
		bLUT[6]=(byte)255;
		for(int i=9;i<256;i++){
			rLUT[i]=(byte)(i-9);
			gLUT[i]=(byte)(i-9);
			bLUT[i]=(byte)(i-9);
		}
		IndexColorModel cm=new IndexColorModel(8,256,rLUT,bLUT,gLUT);
		for(int i=0;i<pixels.length;i++){
			if(jutils.isgray(pixels[i])){
				int[] temp=jutils.intval2rgb(pixels[i]);
				if(temp[0]<247){
					pixels8[i]=(byte)(temp[0]+9);
				}
			}
			switch(pixels[i]){
			case 0xffffffff:
				pixels8[i]=(byte)0;
				break; // white
			case 0xff000000:
				pixels8[i]=(byte)1;
				break; // black
			case 0xff0000ff:
				pixels8[i]=(byte)2;
				break; // blue
			case 0xff00ff00:
				pixels8[i]=(byte)3;
				break; // green
			case 0xffff0000:
				pixels8[i]=(byte)4;
				break; // red
			case 0xffff00ff:
				pixels8[i]=(byte)5;
				break; // magenta
			case 0xff00ffff:
				pixels8[i]=(byte)6;
				break; // cyan
			case 0xffffff00:
				pixels8[i]=(byte)7;
				break; // yellow
			case 0xffffc800:
				pixels8[i]=(byte)8;
				break; // orange
			}
		}
		return new ByteProcessor(cp.getWidth(),cp.getHeight(),pixels8,cm);
	}

	public Image getImage(){
		ColorProcessor cp=getProcessor();
		if(cp!=null)
			return cp.createImage();
		else
			return null;
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
		jdio.writestring(os,"pw2_file_type");
		jdio.writeintelint(os,0); //identifies this as a Plot4
		jdio.writestring(os,getxLabel()); //x label
		jdio.writestring(os,getyLabel()); //y label
		jdio.writeintelint(os,nseries); // number of series'
		jdio.writeintelint(os,getmaxpts()); // max number of pts
		jdio.writeintelfloat(os,xMin); // min x axis
		jdio.writeintelfloat(os,xMax); // max x axis
		jdio.writeintelfloat(os,yMin); // min y axis
		jdio.writeintelfloat(os,yMax); // max y axis
		jdio.writeintelint(os,logx?1:0); // logx?
		jdio.writeintelint(os,logy?1:0); // logy?
		for(int l=0;l<nseries;l++){
			jdio.writeintelint(os,npts[l]); // number of points in this series
			jdio.writeintelint(os,shapes[l]); // shape index
			jdio.writeintelint(os,colors[l]); // color index
			jdio.writeintelfloatarray(os,xValues[l],npts[l]); // x values
			jdio.writeintelfloatarray(os,yValues[l],npts[l]); // y values
		}
		// save the errors if they exist
		if(errors==null)
			showerrors=false;
		jdio.writeintelint(os,showerrors?1:0);
		if(showerrors){
			for(int l=0;l<nseries;l++){
				jdio.writeintelfloatarray(os,errors[0][l],npts[l]); // lower
				// errors
				jdio.writeintelfloatarray(os,errors[1][l],npts[l]); // upper
				// errors
			}
		}
		if(annotations!=null && annotations.length==nseries){
			jdio.writeintelint(os,1);
			for(int i=0;i<nseries;i++) jdio.writestring(os,annotations[i]);
		}
	}

	public Plot4 getCopy(){
		Plot4 temp=new Plot4(xLabel,yLabel,xValues,yValues,npts);
		temp.colors=colors;
		temp.shapes=shapes;
		temp.xMin=xMin;
		temp.yMin=yMin;
		temp.logx=logx;
		temp.logy=logy;
		temp.magnification=magnification;
		temp.magratio=magratio;
		return temp;
	}

	public Plot4 getSelFitCopy(){
		int selindex=selected;
		if(selindex<0)
			selindex=0;
		float[] tempy=new float[npts[selindex]];
		float[] tempx=new float[npts[selindex]];
		System.arraycopy(xValues[selindex],0,tempx,0,npts[selindex]);
		System.arraycopy(yValues[selindex],0,tempy,0,npts[selindex]);
		Plot4 temp=new Plot4(xLabel,yLabel,tempx,tempy);
		temp.xMin=xMin;
		temp.yMin=yMin;
		temp.xMax=xMax;
		temp.yMax=yMax;
		temp.logx=logx;
		temp.logy=logy;
		temp.magnification=magnification;
		temp.magratio=magratio;
		return temp;
	}

}
