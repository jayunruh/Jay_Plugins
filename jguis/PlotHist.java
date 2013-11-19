/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import jalgs.jdataio;

import java.awt.*;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.text.*;

public class PlotHist{
	// this is a class that draws a plot usually associated with a PlotWindow3
	// window
	private float[] histogram;
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
	protected int[] shapes,colors;
	protected String xLabel,yLabel;
	protected boolean logx,logy;
	public static final int WIDTH=512;
	public static final int HEIGHT=200;
	public static final int TICK_LENGTH=3; // length of ticks
	public Color gridColor=new Color(0xc0c0c0); // light gray
	public static final int LEFT_MARGIN=90;
	public static final int RIGHT_MARGIN=18;
	public static final int TOP_MARGIN=20;
	public static final int BOTTOM_MARGIN=50;
	public static final int shapesize=8;
	public static final int fontsize=12;
	protected float magnification,magratio;
	public static final String[] color_names={"black","blue","green","red","magenta","cyan","yellow","orange"};
	public static final String[] shape_names={"line","square","+","x","triangle"};

	public PlotHist(String xLabel1,String yLabel1,float[] xValues1,int color1){
		xValues=new float[][]{xValues1};
		npts=new int[]{xValues1.length};
		xLabel=xLabel1;
		yLabel=yLabel1;
		float[] temp=findminmax(xValues[0]);
		xMin=temp[0];
		xMax=temp[1];
		logx=false;
		logy=false;
		color=color1;
		binSize=4;
		updateHistogram();
		temp=findminmax(histogram);
		yMin=temp[0];
		yMax=temp[1];
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
	}

	public PlotHist(String xLabel1,String yLabel1,float[] hist,float startx,float endx,int color1){
		// here we get the x and y values from an existing histogram
		int sum=0;
		for(int i=0;i<hist.length;i++){
			sum+=(int)hist[i];
		}
		xValues=new float[1][sum];
		int counter=0;
		for(int i=0;i<hist.length;i++){
			float tempxval=startx+((float)i/(float)(hist.length-1))*(endx-startx);
			for(int j=0;j<(int)hist[i];j++){
				xValues[0][counter]=tempxval;
				counter++;
			}
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		xMin=startx;
		xMax=endx;
		logx=false;
		logy=false;
		color=color1;
		binSize=4;
		updateHistogram();
		float[] temp=findminmax(histogram);
		yMin=temp[0];
		yMax=temp[1];
		magnification=1.0f;
		magnification=1.0f;
		magratio=(float)HEIGHT/(float)WIDTH;
		selected=-1;
	}

	public PlotHist(InputStream is){
		init_from_is(is);
	}

	public PlotHist(String filename){
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
		int npts=jdio.readintelint(is);
		xValues=new float[1][npts];
		jdio.readintelfloatfile(is,npts,xValues[0]);
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
		if(temp==3)
			return true;
		else
			return false;
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

	public void setBinSize(float newbinsize){
		binSize=newbinsize;
		updateHistogram();
	}

	public float getBinSize(){
		return binSize;
	}

	public float getBinSizeUnits(){
		return (binSize/histSize)*(xMax-xMin);
	}

	public void setBinSizeUnits(float newbinsize){
		binSize=newbinsize*histSize/(xMax-xMin);
		updateHistogram();
		yautoscale();
	}

	private void updateHistogram(){
		int newhistsize=(int)(histSize/binSize);
		if(!logx){
			float tbinsize=(binSize/histSize)*(xMax-xMin);
			histogram=new float[newhistsize];
			for(int i=0;i<xValues[0].length;i++){
				int dumint1=(int)Math.floor((xValues[0][i]-xMin)/tbinsize);
				if(dumint1>=0&&dumint1<newhistsize){
					histogram[dumint1]++;
				}
			}
		}else{
			// need to implement logarithmic histograming here
			histogram=new float[newhistsize];
			if(xMin<=0.0f){
				logxmin=(float)Math.log((double)findmingt0(xValues[0],xMax));
			}else{
				logxmin=(float)Math.log((double)xMin);
			}
			logxmax=(float)Math.log((double)xMax);
			float tbinsize=(binSize/histSize)*(logxmax-logxmin);
			for(int i=0;i<xValues[0].length;i++){
				float temp=0.0f;
				if(xValues[0][i]>0.0f){
					temp=(float)Math.log(xValues[0][i]);
				}else{
					temp=logxmin;
				}
				int dumint1=(int)Math.floor((temp-logxmin)/tbinsize);
				if(dumint1>=0&&dumint1<newhistsize){
					histogram[dumint1]++;
				}
			}
		}
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		xMin=(float)xMin1;
		xMax=(float)xMax1;
		updateHistogram();
		yMin=(float)yMin1;
		yMax=(float)yMax1;
	}

	public void setLimits(float[] limits){
		xMin=limits[0];
		xMax=limits[1];
		updateHistogram();
		yMin=limits[2];
		yMax=limits[3];
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

	public void autoscale(){
		float[] temp=findminmax(xValues[0]);
		xMin=temp[0];
		xMax=temp[1];
		updateHistogram();
		temp=findminmax(histogram);
		yMin=temp[0];
		yMax=temp[1];
	}

	public void xautoscale(){
		float[] temp=findminmax(xValues[0]);
		xMin=temp[0];
		xMax=temp[1];
		updateHistogram();
	}

	public void yautoscale(){
		float[] temp=findminmax(histogram);
		yMin=temp[0];
		yMax=temp[1];
	}

	public void scalerect(Rectangle rect){
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

	public void updateSeries(float[] xValues1,boolean rescale){
		xValues[0]=xValues1;
		if(rescale){
			autoscale();
		}else{
			updateHistogram();
		}
	}

	public void updateSeries(float[] hist,float startx,float endx,boolean rescale){
		int sum=0;
		for(int i=0;i<hist.length;i++){
			sum+=(int)hist[i];
		}
		xValues[0]=new float[sum];
		int counter=0;
		for(int i=0;i<hist.length;i++){
			float tempxval=startx+((float)i/(float)(hist.length-1))*(endx-startx);
			for(int j=0;j<(int)hist[i];j++){
				xValues[0][counter]=tempxval;
				counter++;
			}
		}
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
				logymin=(float)Math.log((double)findmingt0(histogram,yMax));
			}else{
				logymin=(float)Math.log((double)yMin);
			}
			logymax=(float)Math.log((double)yMax);
			logyscale=(float)newheight/(logymax-logymin);
		}
		// IJ.showMessage("testdraw1");

		drawAxisLabels(pr);
		pr.setClipRect(frame.x,frame.y,frame.width,frame.height);
		xScale=(float)newwidth/(xMax-xMin);
		yScale=(float)newheight/(yMax-yMin);

		// IJ.showMessage("testdraw2");
		Color fillcolor=getColor(color);
		float binwidth=2.0f*magnification*binSize;
		if(!logy){
			for(int i=0;i<histogram.length;i++){
				int temp2=(int)((histogram[i]-yMin)*yScale);
				if(temp2<0)
					temp2=0;
				if(temp2>frame.height)
					temp2=frame.height;
				int ystart=newtopmargin+frame.height-temp2;
				int yend=newtopmargin+frame.height;
				int xstart=(int)((float)i*binwidth+(float)frame.x);
				int xend=(int)((float)i*binwidth+binwidth+(float)frame.x);
				// first fill the rectangle
				pr.setColor(fillcolor);
				pr.fillRect(xstart,ystart,xend-xstart,yend-ystart);
				// then outline it
				pr.setColor(Color.black);
				pr.drawRect(xstart,ystart,xend-xstart,yend-ystart);
			}
		}else{
			pr.setColor(getColor(color));
			for(int i=0;i<histogram.length;i++){
				float ytemp;
				if(histogram[i]>0.0f){
					ytemp=(float)Math.log((double)histogram[i]);
				}else{
					ytemp=logymin;
				}
				int ystart=newtopmargin+frame.height-(int)((ytemp-logymin)*logyscale);
				int yend=newtopmargin+frame.height;
				int xstart=(int)((float)i*binwidth+(float)frame.x);
				int xend=(int)((float)i*binwidth+binwidth+(float)frame.x);
				// first fill the rectangle
				pr.setColor(fillcolor);
				pr.fillRect(xstart,ystart,xend-xstart,yend-ystart);
				// then outline it
				pr.setColor(Color.black);
				pr.drawRect(xstart,ystart,xend-xstart,yend-ystart);
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
		float[] xticklabels=new float[4];
		float[] yticklabels=new float[4];
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
		}

		pr.setAntiAliasedText(true);
		Font tempfont=new Font("SansSerif",Font.PLAIN,newfontsize);
		pr.setFont(tempfont);
		FontMetrics fm=pr.getFontMetrics();
		int fontheight=fm.getAscent()+fm.getDescent();

		// draw the y axis labels

		String s=jutils.formatted_string((double)yticklabels[0]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+newheight+fm.getDescent());
		s=jutils.formatted_string((double)yticklabels[1]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+(int)((2*newheight)/3)+fontheight/2);
		s=jutils.formatted_string((double)yticklabels[2]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+(int)(newheight/3)+fontheight/2);
		s=jutils.formatted_string((double)yticklabels[3]);
		pr.drawString(s,newleftmargin-pr.getStringWidth(s)-newticklength-2,newtopmargin+fontheight/2);

		// draw the x axis labels
		int y=newtopmargin+newheight+fontheight+newticklength+4;
		s=jutils.formatted_string((double)xticklabels[0]);
		pr.drawString(s,newleftmargin-8,y);
		s=jutils.formatted_string((double)xticklabels[1]);
		pr.drawString(s,newleftmargin+(int)(newwidth/3-pr.getStringWidth(s)/2),y);
		s=jutils.formatted_string((double)xticklabels[2]);
		pr.drawString(s,newleftmargin+(int)((2*newwidth)/3-pr.getStringWidth(s)/2),y);
		s=jutils.formatted_string((double)xticklabels[3]);
		pr.drawString(s,newleftmargin+newwidth-pr.getStringWidth(s)+8,y);

		// now draw the axis labels
		pr.drawString(xLabel,newleftmargin+(newwidth-pr.getStringWidth(xLabel))/2,newtopmargin+newheight+fm.getAscent()+newbottommargin/2);
		pr.drawVerticalString(yLabel,newleftmargin-(4*newleftmargin/5),newtopmargin+(int)(newheight/2));
		// now draw the tick marks and grid lines
		pr.drawRect(newleftmargin,newtopmargin,newwidth,newheight);
		pr.setColor(gridColor);
		pr.drawLine(newleftmargin+(int)(newwidth/3),newtopmargin,newleftmargin+(int)(newwidth/3),newtopmargin+newheight);
		pr.drawLine(newleftmargin+(int)((2*newwidth)/3),newtopmargin,newleftmargin+(int)((2*newwidth)/3),newtopmargin+newheight);
		pr.drawLine(newleftmargin,newtopmargin+(int)(newheight/3),newleftmargin+newwidth,newtopmargin+(int)(newheight/3));
		pr.drawLine(newleftmargin,newtopmargin+(int)((2*newheight)/3),newleftmargin+newwidth,newtopmargin+(int)((2*newheight)/3));
		pr.setColor(Color.black);
		pr.drawLine(newleftmargin,newtopmargin,newleftmargin-newticklength,newtopmargin);
		pr.drawLine(newleftmargin,newtopmargin+(int)(newheight/3),newleftmargin-newticklength,newtopmargin+(int)(newheight/3));
		pr.drawLine(newleftmargin,newtopmargin+(int)((2*newheight)/3),newleftmargin-newticklength,newtopmargin+(int)((2*newheight)/3));
		pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin-newticklength,newtopmargin+newheight);

		pr.drawLine(newleftmargin,newtopmargin+newheight,newleftmargin,newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+(int)(newwidth/3),newtopmargin+newheight,newleftmargin+(int)(newwidth/3),newtopmargin+newheight+newticklength);
		pr.drawLine(newleftmargin+(int)((2*newwidth)/3),newtopmargin+newheight,newleftmargin+(int)((2*newwidth)/3),newtopmargin+newheight+newticklength);
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

	public float[] getLimits(){
		float[] temp={xMin,xMax,yMin,yMax};
		return temp;
	}

	public boolean[] getLogAxes(){
		boolean[] temp={logx,logy};
		return temp;
	}

	public float[][] getHistogram(){
		float[][] temp=new float[2][histogram.length];
		float tbinsize=(binSize/histSize)*(xMax-xMin);
		if(!logx){
			for(int i=0;i<histogram.length;i++){
				temp[0][i]=xMin+tbinsize*(float)i;
				temp[1][i]=histogram[i];
			}
		}else{
			for(int i=0;i<histogram.length;i++){
				temp[0][i]=logxmin+((float)i/(float)(histogram.length-1))*(logxmax-logxmin);
				temp[0][i]=(float)Math.exp((double)temp[0][i]);
				temp[1][i]=histogram[i];
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
		jdio.writeintelint(os,3);
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
		jdio.writeintelint(os,xValues[0].length);
		jdio.writeintelfloatarray(os,xValues[0],xValues[0].length);
		// save the errors if they exist
	}

	public Plot4 getSelFitCopy(){
		float[][] hist=getHistogram();
		Plot4 temp=new Plot4(xLabel,yLabel,hist[0],hist[1]);
		temp.setLimits(xMin,xMax,yMin,yMax);
		temp.setLogAxes(logx,logy);
		temp.setmagnification(magnification);
		temp.setmagratio((float)HEIGHT/(float)WIDTH);
		return temp;
	}

}
