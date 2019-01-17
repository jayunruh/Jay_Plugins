/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.process.ColorProcessor;
import jalgs.algutils;
import jalgs.jdataio;
import jalgs.jstatistics;
import jalgs.jfit.fit_kde_gaus;
import jalgs.jsim.rngs;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Rectangle;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DecimalFormat;

public class PlotColumn{
	//this is a class that plots column and whisker plots with optional scatter overlay
	//basically equivalent to a multihistogram plot
	//may eventually support violin plots
	private float[][] histogram;
	private float[][] stats;
	public String uwhiskerstat,uboxstat,lwhiskerstat,lboxstat;
	public float[] uwhiskerextras,uboxextras,lwhiskerextras,lboxextras;
	private int color;
	private float binSize;
	public static final int histSize=256;
	protected float[][] xValues;
	protected float[][] yValues;
	protected float[][][] errors;
	protected boolean showerrors;
	protected int[] npts;
	protected float xMin,xMax,yMin,yMax,xScale,yScale;
	protected float logxmin,logymin,logxscale,logyscale,logxmax,logymax;
	protected int maxpts,nseries;
	protected int selected;
	protected int[] shapes,colors,types;
	protected String xLabel,yLabel;
	protected boolean logx,logy;
	protected String[] annotations;
	public static final int WIDTH=400;
	public static final int HEIGHT=200;
	public static final int TICK_LENGTH=3; // length of ticks
	public static final int LEFT_MARGIN=90;
	public static final int RIGHT_MARGIN=18;
	public static final int TOP_MARGIN=20;
	public static final int BOTTOM_MARGIN=50;
	public int shapesize=8;
	public int fontsize=14;
	public boolean centered_data=true;
	protected float magnification,magratio;
	public static final String[] color_names={"black","blue","green","red","magenta","cyan","yellow","orange"};
	public static final String[] shape_names={"line","square","+","x","triangle"};
	public static final String[] type_names={"column","box_whisker","box_whisker_data","violin","violin_data"};

	public PlotColumn(String xLabel1,String yLabel1,float[] xValues1,int color1){
		xValues=new float[][]{xValues1};
		npts=new int[]{xValues1.length};
		nseries=1;
		maxpts=xValues1.length;
		xLabel=xLabel1;
		yLabel=yLabel1;
		float[] temp=findminmax(xValues[0]);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		color=color1;
		binSize=4;
		updateHistogram();
		//temp=findminmax(histogram[0]);
		//yMin=temp[0];
		//yMax=temp[1];
		xMin=0.0f;
		xMax=1.0f;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
	}
	
	public PlotColumn(String xLabel1,String yLabel1,float[][] xValues1,int[] npts,int color1){
		xValues=xValues1;
		nseries=xValues1.length;
		maxpts=xValues1[0].length;
		if(npts!=null) this.npts=npts;
		else {
			for(int i=0;i<nseries;i++) npts[i]=maxpts;
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		float[] temp=findminmax(xValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		logx=false;
		logy=false;
		color=color1;
		binSize=4;
		updateHistogram();
		//temp=findminmax(histogram);
		//yMin=temp[0];
		//yMax=temp[1];
		xMin=0.0f;
		xMax=nseries;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
	}

	public PlotColumn(InputStream is){
		init_from_is(is);
	}

	public PlotColumn(String filename){
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
		logx=jdio.readintelint(is)==1;
		logy=jdio.readintelint(is)==1;
		color=jdio.readintelint(is);
		binSize=jdio.readintelfloat(is);
		//int npts=jdio.readintelint(is);
		nseries=jdio.readintelint(is);
		maxpts=jdio.readintelint(is);
		centered_data=jdio.readintelint(is)==1;
		npts=new int[nseries];
		xValues=new float[nseries][maxpts];
		types=new int[nseries];
		shapes=new int[nseries];
		colors=new int[nseries];
		for(int i=0;i<nseries;i++){
			types[i]=jdio.readintelint(is);
			shapes[i]=jdio.readintelint(is);
			colors[i]=jdio.readintelint(is);
			npts[i]=jdio.readintelint(is);
			jdio.readintelfloatfile(is,npts[i],xValues[i]);
		}
		boolean annotated=jdio.readintelint(is)==1;
		if(annotated){
			annotations=new String[nseries];
			for(int i=0;i<nseries;i++){
				annotations[i]=jdio.readstring(is);
			}
		}
		updateHistogram();
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
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
		if(temp==5)
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
	
	protected float[] findminmax(float[][] arr){
		int[] tempnpts=new int[arr.length];
		for(int i=0;i<arr.length;i++) tempnpts[i]=arr[i].length;
		return findminmax(arr,tempnpts);
	}

	private float[] findminmax(float[] arr){
		float[] temp=new float[2];
		temp[0]=arr[0];
		temp[1]=arr[0];
		for(int i=0;i<arr.length;i++){
			if(arr[i]<temp[0]){
				temp[0]=arr[i];
			}
			if(arr[i]>temp[1]){
				temp[1]=arr[i];
			}
		}
		return temp;
	}

	private float findmingt0(float[] arr,float max){
		float temp=max;
		if(max<=0.0f){
			return 0.0f;
		}
		for(int i=0;i<arr.length;i++){
			if(arr[i]<temp&&arr[i]>0.0f){
				temp=arr[i];
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
	
	protected float findmingt0(float[][] arr,float max){
		int[] tempnpts=new int[arr.length];
		for(int i=0;i<arr.length;i++) tempnpts[i]=arr[i].length;
		return findmingt0(arr,tempnpts,max);
	}

	public void setBinSize(float newbinsize){
		binSize=newbinsize;
		updateHistogram();
	}

	public float getBinSize(){
		return binSize;
	}

	public float getBinSizeUnits(){
		return (binSize/histSize)*(yMax-yMin);
	}

	public void setBinSizeUnits(float newbinsize){
		binSize=newbinsize*histSize/(yMax-yMin);
		updateHistogram();
		//yautoscale();
	}
	
	public void setFontSize(int fontsize){
		this.fontsize=fontsize;
	}
	
	public void setShapeSize(int shapesize){
		this.shapesize=shapesize;
	}
	
	public void setAnnotations(String[] annotations){
		this.annotations=annotations;
	}

	private void updateHistogram(){
		/*int newhistsize=(int)(histSize/binSize);
		if(!logy){
			float tbinsize=(binSize/histSize)*(yMax-yMin);
			histogram=new float[nseries][newhistsize];
			for(int j=0;j<nseries;j++){
    			for(int i=0;i<npts[j];i++){
    				int dumint1=(int)Math.floor((xValues[j][i]-yMin)/tbinsize);
    				if(dumint1>=0&&dumint1<newhistsize){
    					histogram[j][dumint1]++;
    				}
    			}
			}
		}else{
			// need to implement logarithmic histograming here
			histogram=new float[nseries][newhistsize];
			if(yMin<=0.0f){
				logymin=(float)Math.log((double)findmingt0(xValues,npts,yMax));
			}else{
				logymin=(float)Math.log((double)yMin);
			}
			logxmax=(float)Math.log((double)yMax);
			float tbinsize=(binSize/histSize)*(logymax-logymin);
			for(int j=0;j<nseries;j++){
    			for(int i=0;i<npts[j];i++){
    				float temp=0.0f;
    				if(xValues[j][i]>0.0f){
    					temp=(float)Math.log(xValues[j][i]);
    				}else{
    					temp=logymin;
    				}
    				int dumint1=(int)Math.floor((temp-logymin)/tbinsize);
    				if(dumint1>=0&&dumint1<newhistsize){
    					histogram[j][dumint1]++;
    				}
    			}
			}
		}*/
		//the statistics are the mean, median, box stat lower, upper (sterr), and whisker stat lower, upper (max,min)
		stats=new float[nseries][6];
		if(uwhiskerstat==null || uwhiskerstat.length()<2) uwhiskerstat="Max";
		if(uboxstat==null || uboxstat.length()<2) uboxstat="StErr";
		if(lwhiskerstat==null || lwhiskerstat.length()<2) lwhiskerstat="Min";
		if(lboxstat==null || lboxstat.length()<2) lboxstat="StErr";
		for(int i=0;i<nseries;i++){
			Rectangle r=new Rectangle(0,0,npts[i],1);
			stats[i][0]=jstatistics.getstatistic("Avg",xValues[i],npts[i],1,r,null);
			stats[i][1]=jstatistics.getstatistic("Median",xValues[i],npts[i],1,r,null);
			stats[i][2]=jstatistics.getstatistic(lboxstat,xValues[i],npts[i],1,r,lboxextras);
			stats[i][3]=jstatistics.getstatistic(uboxstat,xValues[i],npts[i],1,r,uboxextras);
			stats[i][4]=jstatistics.getstatistic(lwhiskerstat,xValues[i],npts[i],1,r,lwhiskerextras);
			stats[i][5]=jstatistics.getstatistic(uwhiskerstat,xValues[i],npts[i],1,r,uwhiskerextras);
			if(!lwhiskerstat.equals("Min")) stats[i][4]=stats[i][0]-stats[i][4];
			if(!uwhiskerstat.equals("Max")) stats[i][5]=stats[i][0]+stats[i][5];
			if(!lboxstat.equals("Min")) stats[i][2]=stats[i][0]-stats[i][2];
			if(!uboxstat.equals("Max")) stats[i][3]=stats[i][0]+stats[i][3];
		}
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		yMin=(float)yMin1;
		yMax=(float)yMax1;
		updateHistogram();
		xMin=(float)xMin1;
		xMax=(float)xMax1;
	}

	public void setLimits(float[] limits){
		yMin=limits[2];
		yMax=limits[3];
		updateHistogram();
		xMin=limits[0];
		xMax=limits[1];
		//IJ.log(""+yMin+" , "+yMax);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1){
		logx=logx1;
		logy=logy1;
		updateHistogram();
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
	
	public void setAllColors(int color){
		if(colors==null) colors=new int[nseries];
		for(int i=0;i<nseries;i++) colors[i]=color;
	}
	
	public void setAllShapes(int shape){
		if(shapes==null) shapes=new int[nseries];
		for(int i=0;i<nseries;i++) shapes[i]=shape;
	}
	
	public void setAllTypes(int type){
		if(types==null) types=new int[nseries];
		for(int i=0;i<nseries;i++) types[i]=type;
	}

	public void autoscale(){
		xMin=0.0f;
		xMax=nseries;
		float[] temp=findminmax(xValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		updateHistogram();
	}

	public void xautoscale(){
		xMin=0.0f;
		xMax=nseries;
	}

	public void yautoscale(){
		float[] temp=findminmax(xValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		updateHistogram();
	}

	/*public void scalerect(Rectangle rect){
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		if(rect!=null){
			float tempxmin=((float)(rect.x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymax=((float)(HEIGHT*ymag+TOP_MARGIN*ymag-rect.y)/(float)(HEIGHT*ymag))*(yMax-yMin)+yMin;
			float tempxmax=((float)(rect.x+rect.width-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempymin=((float)(HEIGHT*ymag+TOP_MARGIN*ymag-rect.y-rect.height)/(float)(HEIGHT*ymag))*(yMax-yMin)+yMin;
			if(logx){
				float templogxmin=((float)(rect.x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmin=(float)Math.exp(templogxmin);
				float templogxmax=((float)(rect.x+rect.width-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmax=(float)Math.exp(templogxmax);
			}
			if(logy){
				float templogymax=((float)(HEIGHT*ymag+TOP_MARGIN*ymag-rect.y)/(float)(HEIGHT*ymag))*(logymax-logymin)+logymin;
				tempymax=(float)Math.exp(templogymax);
				float templogymin=((float)(HEIGHT*ymag+TOP_MARGIN*ymag-rect.y-rect.height)/(float)(HEIGHT*ymag))*(logymax-logymin)+logymin;
				tempymin=(float)Math.exp(templogymin);
			}
			xMin=tempxmin;
			xMax=tempxmax;
			yMin=tempymin;
			yMax=tempymax;
			updateHistogram();
		}
	}

	public int[] getrectindices(Rectangle rect){
		if(rect!=null){
			float tempxmin=((float)(rect.x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(xMax-xMin)+xMin;
			float tempxmax=((float)(rect.x+rect.width-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(xMax-xMin)+xMin;
			if(logx){
				float templogxmin=((float)(rect.x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
				tempxmin=(float)Math.exp(templogxmin);
			}
			int[] temp=new int[xValues.length];
			int counter=0;
			for(int i=0;i<xValues.length;i++){
				if(xValues[0][i]>=tempxmin&&xValues[0][i]<tempxmax){
					temp[counter]=i;
					counter++;
				}
			}
			int[] indices=new int[counter];
			System.arraycopy(temp,0,indices,0,counter);
			return indices;
		}else{
			return null;
		}
	}
	
	public float getPlotCoords(int x){
		float tempx=((float)(x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(xMax-xMin)+xMin;
		if(logx){
			float templogx=((float)(x-LEFT_MARGIN*magnification)/(float)(WIDTH*magnification))*(logxmax-logxmin)+logxmin;
			tempx=(float)Math.exp(templogx);
		}
		return tempx;
	}
	
	public int getPlotPosition(float x){
		// this maps a value to the image
		float posx=(int)((x-xMin)/(xMax-xMin)*(float)histSize);
		if(logx){
			// assume that logxmin has already been updated
			float temp=0.0f;
			if(x>0.0f)
				temp=(float)Math.log(x);
			else
				temp=logxmin;
			posx=((temp-logxmin)/(logxmax-logxmin)*(float)histSize);
		}
		posx*=magnification;
		posx+=LEFT_MARGIN*magnification;
		return (int)posx;
	}
	
	public float get_score(Roi roi){
		if(roi!=null){
			int counter=0;
			for(int i=0;i<npts[0];i++){
				int location=getPlotPosition(xValues[0][i]);
				Rectangle r=roi.getBounds();
				if(location>=r.x && location<(r.x+r.width)){
					counter++;
				}
			}
			return (float)counter/(float)npts[0];
		}
		return 0.0f;
	}*/

	public void updateSeries(float[] xValues1,int series,boolean rescale){
		if(xValues1.length>maxpts){
			float[][] newxvals=new float[nseries][xValues1.length];
			for(int i=0;i<nseries;i++){
				if(i==series) System.arraycopy(xValues1,0,newxvals[i],0,xValues1.length);
				else System.arraycopy(xValues[i],0,newxvals[i],0,maxpts);
			}
			maxpts=xValues1.length;
			xValues=newxvals;
		}
		xValues[series]=xValues1;
		if(rescale){
			autoscale();
		}else{
			updateHistogram();
		}
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
		//plot types are column/error, box/whisker, box/whisker with data, and violin
		//note that logx is really logy.  There is no such thing as logx
		float ymag=magnification*magratio/((float)HEIGHT/(float)WIDTH);
		int newleftmargin=(int)(magnification*LEFT_MARGIN);
		int newtopmargin=(int)(ymag*TOP_MARGIN);
		int newwidth=(int)(magnification*WIDTH);
		int newheight=(int)(ymag*HEIGHT);
		Rectangle frame=new Rectangle(newleftmargin,newtopmargin,newwidth,newheight);
		logymin=0;
		logyscale=0;
		logymax=0;
		if(logy){
			if(yMin<=0.0f){
				logymin=(float)Math.log(findmingt0(xValues,npts,yMax));
			}else{
				logymin=(float)Math.log(yMin);
			}
			logymax=(float)Math.log(yMax);
			logyscale=newheight/(logymax-logymin);
		}
		drawAxisLabels(pr);
		pr.setClipRect(frame.x,frame.y,frame.width,frame.height);
		xScale=newwidth/(xMax-xMin);
		yScale=newheight/(yMax-yMin);
		//IJ.log(""+yMin+" , "+yMax);

		// IJ.showMessage("testdraw2");
		//Color fillcolor=getColor(color);
		float binwidth=2.0f*magnification*binSize;
		float widthfraction=0.75f;
		int colwidth=(int)(widthfraction*xScale);
		if(types==null){
			types=new int[nseries];
			for(int i=0;i<nseries;i++) types[i]=2; //default is box + data
		}
		if(colors==null){
			colors=new int[nseries];
			for(int i=0;i<nseries;i++) colors[i]=color; //default is black
		}
		if(shapes==null){
			shapes=new int[nseries];
			for(int i=0;i<nseries;i++) shapes[i]=3; //default is x
		}

		rngs random=new rngs();
		for(int i=0;i<nseries;i++){
			//types[i]=2;
			//shapes[i]=3;
			if(types[i]==0){ //column type
    			float xval=i+0.5f;
    			int xstart=(int)(xval*xScale+frame.x-0.5f*colwidth);
    			int xend=xstart+colwidth;
    			int temp2=(int)((stats[i][0]-yMin)*yScale);
    			if(logy) {
    				temp2=(int)(((float)Math.log(stats[i][0])-logymin)*logyscale);
    			}
    			if(temp2<0) temp2=0;
    			if(temp2>frame.height) temp2=frame.height;
    			int ystart=newtopmargin+frame.height-temp2;
    			int yend=newtopmargin+frame.height;
    			// first fill the rectangle
    			Color fillcolor=getColor(colors[i]);
    			//IJ.log(fillcolor.toString());
    			pr.setColor(fillcolor);
    			pr.fillRect(xstart,ystart,xend-xstart,yend-ystart);
    			// then outline it
    			pr.setColor(Color.black);
    			pr.drawRect(xstart,ystart,xend-xstart,yend-ystart);
    			//now add the error bars
    			int yl=newtopmargin+frame.height-(int)((stats[i][2]-yMin)*yScale);
    			int yu=newtopmargin+frame.height-(int)((stats[i][3]-yMin)*yScale);
    			if(logy) {
    				yl=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][2])-logymin)*logyscale);
    				yu=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][3])-logymin)*logyscale);
    			}
    			int xc=xstart+(int)(0.5f*colwidth);
    			pr.drawErrors(xc,yu,yl);
			}
			if(types[i]>0 && types[i]<3){ //box and whisker type
    			float xval=i+0.5f;
    			int xstart=(int)(xval*xScale+frame.x-colwidth/2);
    			int xend=xstart+colwidth;
    			int temp2=(int)((stats[i][0]-yMin)*yScale);
    			if(logy) {
    				temp2=(int)(((float)Math.log(stats[i][0])-logymin)*logyscale);
    			}
    			if(temp2<0) temp2=0;
    			if(temp2>frame.height) temp2=frame.height;
    			int ystart=newtopmargin+frame.height-temp2;
    			int yend=newtopmargin+frame.height;
    			int xc=xstart+(int)(0.5f*colwidth);
    			pr.setColor(getColor(colors[i]));
    			//first draw the mean
    			pr.drawSquare(xc,ystart);
    			//now the median
    			int temp3=(int)((stats[i][1]-yMin)*yScale);
    			if(logy) {
    				temp3=(int)(((float)Math.log(stats[i][1])-logymin)*logyscale);
    			}
    			if(temp3<0) temp3=0;
    			if(temp3>frame.height) temp3=frame.height;
    			int medstart=newtopmargin+frame.height-temp3;
    			pr.drawLine(xstart,medstart,xend,medstart);
    			//now add the box
    			int yl=newtopmargin+frame.height-(int)((stats[i][2]-yMin)*yScale);
    			int yu=newtopmargin+frame.height-(int)((stats[i][3]-yMin)*yScale);
    			if(logy) {
    				yl=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][2])-logymin)*logyscale);
    				yu=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][3])-logymin)*logyscale);
    			}
    			pr.drawRect(xstart,yl,xend-xstart,yu-yl);
    			//and the whiskers
    			yl=newtopmargin+frame.height-(int)((stats[i][4]-yMin)*yScale);
    			yu=newtopmargin+frame.height-(int)((stats[i][5]-yMin)*yScale);
    			if(logy) {
    				yl=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][4])-logymin)*logyscale);
    				yu=newtopmargin+frame.height-(int)(((float)Math.log(stats[i][5])-logymin)*logyscale);
    			}
    			pr.drawErrors(xc,yu,yl);
    			if(types[i]==2){ //add the data points on top of the plot
    				int[] newxvals=new int[npts[i]];
    				int[] newyvals=new int[npts[i]];
    				for(int j=0;j<npts[i];j++){
    					newyvals[j]=newtopmargin+frame.height-(int)((xValues[i][j]-yMin)*yScale);
    					if(logy) newyvals[j]=newtopmargin+frame.height-(int)(((float)Math.log(xValues[i][j])-logymin)*logyscale);
    					if(!centered_data) {
    						if(j%2==0) newxvals[j]=(int)(random.unidev(xend,xc)+0.5);
    						else newxvals[j]=(int)(random.unidev(xc,xstart)+0.5);
    					}
    					else newxvals[j]=(int)(0.5f*(float)(xend+xstart)+0.5f);
    				}
    				pr.drawPolyshape(newxvals,newyvals,shapes[i],npts[i]);
    			}
			}
			if(types[i]>2) { //this is a violin plot type
				float xval=i+0.5f;
    			int xstart=(int)(xval*xScale+frame.x-colwidth/2);
    			int xend=xstart+colwidth;
    			int temp2=(int)((stats[i][0]-yMin)*yScale);
    			if(logy) {
    				temp2=(int)(((float)Math.log(stats[i][0])-logymin)*logyscale);
    			}
    			if(temp2<0) temp2=0;
    			if(temp2>frame.height) temp2=frame.height;
    			int ystart=newtopmargin+frame.height-temp2;
    			int yend=newtopmargin+frame.height;
    			int xc=xstart+(int)(0.5f*colwidth);
    			pr.setColor(getColor(colors[i]));
    			//first draw the mean
    			pr.drawSquare(xc,ystart);
    			//get the kde curve
    			float[] xvals=(float[])algutils.get_subarray(xValues[i],0,npts[i]);
    			float[][] kde=null;
    			if(logy) {
    				float[] tempxvals=xvals.clone();
    				for(int j=0;j<npts[i];j++) {
    					if(tempxvals[j]<=0.0f) tempxvals[j]=logxmin; //would be better to eliminate it perhaps?
    					tempxvals[j]=(float)Math.log(tempxvals[j]);
    				}
    				kde=fit_kde_gaus.get_kde(tempxvals);
    			} else {
    				kde=fit_kde_gaus.get_kde(xvals);
    			}
    			//new PlotWindow4("kde test","x","y",kde[0],kde[1]).draw();
    			float kdexspace=kde[0][1]-kde[0][0];
    			float firstx=kde[0][0]-kdexspace;
    			float lastx=kde[0][kde[0].length-1]+kdexspace;
    			float kdemax=jstatistics.getstatistic("Max",kde[1],null);
    			int[][] kdexvals=new int[2][kde[0].length+2];
    			int[] kdeyvals=new int[kde[0].length+2];
    			if(!logy) {
        			for(int j=0;j<kde[0].length;j++) {
        				kdeyvals[j+1]=newtopmargin+frame.height-(int)((kde[0][j]-yMin)*yScale);
        				kdexvals[0][j+1]=xc+(int)(0.5f*(float)colwidth*kde[1][j]/kdemax);
        				kdexvals[1][j+1]=xc-(int)(0.5f*(float)colwidth*kde[1][j]/kdemax);
        			}
        			kdeyvals[0]=newtopmargin+frame.height-(int)((firstx-yMin)*yScale);
        			kdeyvals[kdeyvals.length-1]=newtopmargin+frame.height-(int)((lastx-yMin)*yScale);
        			kdexvals[0][0]=xc; kdexvals[0][kdeyvals.length-1]=xc; kdexvals[1][0]=xc; kdexvals[1][kdeyvals.length-1]=xc;
    			} else {
        			for(int j=0;j<kde[0].length;j++) {
        				kdeyvals[j+1]=newtopmargin+frame.height-(int)((kde[0][j]-logymin)*logyscale);
        				kdexvals[0][j+1]=xc+(int)(0.5f*(float)colwidth*kde[1][j]/kdemax);
        				kdexvals[1][j+1]=xc-(int)(0.5f*(float)colwidth*kde[1][j]/kdemax);
        			}
        			kdeyvals[0]=newtopmargin+frame.height-(int)((firstx-logymin)*logyscale);
        			kdeyvals[kdeyvals.length-1]=newtopmargin+frame.height-(int)((lastx-logymin)*logyscale);
        			kdexvals[0][0]=xc; kdexvals[0][kdeyvals.length-1]=xc; kdexvals[1][0]=xc; kdexvals[1][kdeyvals.length-1]=xc;
    			}
    			pr.drawPolyLine(kdexvals[0],kdeyvals,false);
    			pr.drawPolyLine(kdexvals[1],kdeyvals,false);
    			//now interpolate along the kde curve to draw the plot
    			if(types[i]==4){ //add the data points on top of the plot
    				int[] newxvals=new int[npts[i]];
    				int[] newyvals=new int[npts[i]];
    				for(int j=0;j<npts[i];j++){
    					newyvals[j]=newtopmargin+frame.height-(int)((xValues[i][j]-yMin)*yScale);
    					if(logy) newyvals[j]=newtopmargin+frame.height-(int)(((float)Math.log(xValues[i][j])-logymin)*logyscale);
    					if(!centered_data) {
    						if(j%2==0) newxvals[j]=(int)(random.unidev(xend,xc)+0.5);
    						else newxvals[j]=(int)(random.unidev(xc,xstart)+0.5);
    					}
    					else newxvals[j]=(int)(0.5f*(float)(xend+xstart)+0.5f);
    				}
    				pr.drawPolyshape(newxvals,newyvals,shapes[i],npts[i]);
    			}
			}
		}
		pr.setColor(Color.black);
		pr.unsetClip();
		pr.endPlotting();
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
		float[] yticklabels=new float[4];
		for(int i=0;i<4;i++){
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
		for(int i=0;i<nseries;i++){
			if(annotations!=null) s=annotations[i];
			else s=""+(i+1);
			pr.drawString(s,newleftmargin+(int)((0.5f+i)*newwidth/nseries-0.5f*pr.getStringWidth(s)),y);
		}

		// now draw the axis labels
		pr.drawString(xLabel,newleftmargin+(newwidth-pr.getStringWidth(xLabel))/2,newtopmargin+newheight+fm.getAscent()+newbottommargin/2);
		pr.drawVerticalString(yLabel,newleftmargin-(4*newleftmargin/5),newtopmargin+newheight/2);
		// now draw the plot outline
		pr.drawRect(newleftmargin,newtopmargin,newwidth,newheight);
		pr.setColor(Color.black);
		//and now the left tick marks
		pr.drawLine(newleftmargin,newtopmargin,newleftmargin-newticklength,newtopmargin);
		pr.drawLine(newleftmargin,newtopmargin+newheight/3,newleftmargin-newticklength,newtopmargin+newheight/3);
		pr.drawLine(newleftmargin,newtopmargin+(2*newheight)/3,newleftmargin-newticklength,newtopmargin+(2*newheight)/3);
		pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin-newticklength,newtopmargin+newheight);
		//and the bottom
		/*pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin,newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+(int)(newwidth/3),newtopmargin+newheight,newleftmargin+(int)(newwidth/3),newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+(int)((2*newwidth)/3),newtopmargin+newheight,newleftmargin+(int)((2*newwidth)/3),newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+newwidth,newtopmargin+newheight,newleftmargin+newwidth,newtopmargin+newheight+newticklength);*/
	}

	Color getColor(int index){
		int temp=index;
		if(temp>=8){
			temp=index%8;
		}
		Color[] temp2={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};
		return temp2[temp];
	}
	
	public int[] getColors(){
		return colors;
	}
	
	public int[] getTypes(){
		return colors;
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

	public float[] getYValues(int series){
		return null;
	}

	public float[] getErrors(int series,boolean upper){
		return null;
	}

	public void selectSeries(int series){
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
	
	public int getmaxpts(){
		return maxpts;
	}
	
	public int getNSeries(){
		return nseries;
	}

	public float[] getLimits(){
		float[] temp={xMin,xMax,yMin,yMax};
		return temp;
	}

	public boolean[] getLogAxes(){
		boolean[] temp={logx,logy};
		return temp;
	}
	
	public String[] getAnnotations(){
		return annotations;
	}

	public float[][] getHistogram(){
		float[][] temp=new float[2*nseries][histogram.length];
		float tbinsize=(binSize/histSize)*(yMax-yMin);
		if(!logy){
			for(int j=0;j<nseries;j++){
    			for(int i=0;i<histogram[0].length;i++){
    				temp[2*j][i]=yMin+tbinsize*i;
    				temp[2*j+1][i]=histogram[j][i];
    			}
			}
		}else{
			for(int j=0;j<nseries;j++){
    			for(int i=0;i<histogram.length;i++){
    				temp[2*j][i]=logymin+((float)i/(float)(histogram[0].length-1))*(logymax-logymin);
    				temp[2*j][i]=(float)Math.exp(temp[2*j][i]);
    				temp[2*j+1][i]=histogram[j][i];
    			}
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
	
	public byte[] getPDFBinary(){
		ByteArrayOutputStream os=new ByteArrayOutputStream();
		Plotter pr=new PDFPlotter(getWidth(),getHeight(),os);
		drawPlot(pr);
		return os.toByteArray();
	}

	public ColorProcessor getProcessor(){
		Plotter pr=new CPPlotter(getWidth(),getHeight());
		drawPlot(pr);
		return ((CPPlotter)pr).get_output();
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
		jdio.writeintelint(os,5);
		jdio.writestring(os,getxLabel());
		jdio.writestring(os,getyLabel());
		jdio.writeintelfloat(os,xMin); // min x axis
		jdio.writeintelfloat(os,xMax); // max x axis
		jdio.writeintelfloat(os,yMin); // min y axis
		jdio.writeintelfloat(os,yMax); // max y axis
		jdio.writeintelint(os,logx?1:0); // logx?
		jdio.writeintelint(os,logy?1:0); // logy?
		jdio.writeintelint(os,color);
		jdio.writeintelfloat(os,binSize);
		jdio.writeintelint(os,nseries);
		jdio.writeintelint(os,maxpts);
		jdio.writeintelint(os,(centered_data?1:0));
		for(int i=0;i<nseries;i++){
			jdio.writeintelint(os,types[i]);
			jdio.writeintelint(os,shapes[i]);
			jdio.writeintelint(os,colors[i]);
    		jdio.writeintelint(os,npts[i]);
    		jdio.writeintelfloatarray(os,xValues[i],npts[i]);
		}
		if(annotations!=null && annotations.length==nseries){
			jdio.writeintelint(os,1);
			for(int i=0;i<nseries;i++) jdio.writestring(os,annotations[i]);
		}
	}

	/*public Plot4 getSelFitCopy(){
		float[][] hist=getHistogram();
		Plot4 temp=new Plot4(xLabel,yLabel,hist[0],hist[1]);
		temp.setLimits(xMin,xMax,yMin,yMax);
		temp.setLogAxes(logx,logy);
		temp.setmagnification(magnification);
		temp.setmagratio((float)HEIGHT/(float)WIDTH);
		return temp;
	}*/

}
