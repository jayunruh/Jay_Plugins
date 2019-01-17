/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import j3D.renderer;
import jalgs.jdataio;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.List;

public class Traj3D extends Plot3D{
	// here we have a single zValue per xy value
	// private Label coordinates;
	private float[][] xValues;
	private float[][] yValues;
	private float[][] zValues;
	private int[] npts;
	public boolean thick=true;
	public boolean drawIndex;

	public Traj3D(String xLabel1,String yLabel1,String zLabel1,float[][] xValues1,float[][] yValues1,float[][] zValues1,Object npts1){
		super();
		xValues=xValues1;
		yValues=yValues1;
		zValues=zValues1;
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		nseries=zValues.length;
		if(npts1 instanceof int[]){
			npts=(int[])npts1;
			maxxpts=npts[0];
			for(int i=1;i<nseries;i++){
				if(npts[i]>maxxpts){
					maxxpts=npts[i];
				}
			}
		}else{
			npts=new int[nseries];
			npts[0]=zValues[0].length;
			maxxpts=npts[0];
			for(int i=1;i<nseries;i++){
				npts[i]=npts[0];
			}
		}
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax(zValues,npts);
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

	public Traj3D(String xLabel1,String yLabel1,String zLabel1,float[] xValues1,float[] yValues1,float[] zValues1){
		super();
		xValues=new float[1][zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			xValues[0][i]=xValues1[i];
		}
		yValues=new float[1][zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			yValues[0][i]=yValues1[i];
		}
		zValues=new float[1][zValues1.length];
		for(int i=0;i<zValues1.length;i++){
			zValues[0][i]=zValues1[i];
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		maxxpts=zValues1.length;
		npts=new int[1];
		npts[0]=maxxpts;
		nseries=1;
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax(zValues,npts);
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
	
	public Traj3D(String xLabel1,String yLabel1,String zLabel1,List<List<float[]>> points){
		super();
		maxxpts=0;
		npts=new int[points.size()];
		for(int i=0;i<points.size();i++){
			npts[i]=points.get(i).size();
			if(npts[i]>maxxpts) maxxpts=npts[i];
		}
		nseries=npts.length;
		xValues=new float[npts.length][maxxpts];
		yValues=new float[npts.length][maxxpts];
		zValues=new float[npts.length][maxxpts];
		for(int i=0;i<points.size();i++){
			List<float[]> series=points.get(i);
			for(int j=0;j<series.size();j++){
				float[] coords=series.get(j);
				xValues[i][j]=coords[0]; yValues[i][j]=coords[1]; zValues[i][j]=coords[2];
			}
		}
		xLabel=xLabel1;
		yLabel=yLabel1;
		zLabel=zLabel1;
		
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax(zValues,npts);
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
	

	public Traj3D(InputStream is){
		init_from_is(is);
	}

	public Traj3D(String filename){
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
		maxypts=0;
		npts=new int[nseries];
		xValues=new float[nseries][maxxpts];
		yValues=new float[nseries][maxxpts];
		zValues=new float[nseries][maxxpts];
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
			npts[l]=jdio.readintelint(is);
			shapes[l]=jdio.readintelint(is);
			colors[l]=jdio.readintelint(is);
			jdio.readintelfloatfile(is,npts[l],xValues[l]);
			jdio.readintelfloatfile(is,npts[l],yValues[l]);
			jdio.readintelfloatfile(is,npts[l],zValues[l]); // y values
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
		if(temp==2)
			return true;
		else
			return false;
	}
	
	public static boolean is_this(InputStream is) throws IOException{
		jdataio jdio=new jdataio();
		jdio.readstring(is); // read the label
		int temp=jdio.readintelint(is); // now the identifier
		return (temp==2);
	}

	public void autoscale(){
		float[] temp=findminmax(xValues,npts);
		xMin=temp[0];
		xMax=temp[1];
		temp=findminmax(yValues,npts);
		yMin=temp[0];
		yMax=temp[1];
		temp=findminmax(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
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

	public void zautoscale(){
		float[] temp=findminmax(zValues,npts);
		zMin=temp[0];
		zMax=temp[1];
	}

	public void updateSeries(float[] xValues1,float[] yValues1,float[] zValues1,int series,boolean rescale){
		int xlength=zValues1.length;
		npts[series]=xlength;
		if(xlength>maxxpts){
			if(xlength>maxxpts){
				maxxpts=xlength;
			}
			float[][] newxValues=new float[nseries][maxxpts];
			float[][] newyValues=new float[nseries][maxxpts];
			float[][] newzValues=new float[nseries][maxxpts];
			for(int i=0;i<series;i++){
				for(int j=0;j<npts[i];j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
					newzValues[i][j]=zValues[i][j];
				}
			}
			for(int j=0;j<npts[series];j++){
				newxValues[series][j]=xValues1[j];
				newyValues[series][j]=yValues1[j];
				newzValues[series][j]=zValues1[j];
			}
			for(int i=series+1;i<nseries;i++){
				for(int j=0;j<npts[i];j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
					newzValues[i][j]=zValues[i][j];
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
				yValues[series][i]=yValues1[i];
				zValues[series][i]=zValues1[i];
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
		float[][] newzValues;
		int[] newnpts=new int[nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		int newmaxxpts=0;
		if(npts[series]==maxxpts){
			for(int i=0;i<=nseries;i++){
				if(i!=series){
					if(npts[i]>newmaxxpts){
						newmaxxpts=npts[i];
					}
				}
			}
		}else{
			newmaxxpts=maxxpts;
		}
		newxValues=new float[nseries][newmaxxpts];
		newyValues=new float[nseries][newmaxxpts];
		newzValues=new float[nseries][newmaxxpts];
		for(int i=0;i<series;i++){
			newnpts[i]=npts[i];
			newshapes[i]=shapes[i];
			newcolors[i]=colors[i];
			for(int j=0;j<newmaxxpts;j++){
				newxValues[i][j]=xValues[i][j];
				newyValues[i][j]=yValues[i][j];
				newzValues[i][j]=zValues[i][j];
			}
		}
		for(int i=series+1;i<=nseries;i++){
			newnpts[i-1]=npts[i];
			newshapes[i-1]=shapes[i];
			newcolors[i-1]=colors[i];
			for(int j=0;j<newmaxxpts;j++){
				newxValues[i-1][j]=xValues[i][j];
				newyValues[i-1][j]=yValues[i][j];
				newzValues[i-1][j]=zValues[i][j];
			}
		}
		maxxpts=newmaxxpts;
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

	public void addPoints(float[] xValues1,float[] yValues1,float[] zValues1,boolean rescale){
		nseries++;
		float[][] newxValues;
		float[][] newyValues;
		float[][] newzValues;
		int[] newnpts=new int[nseries];
		int[] newshapes=new int[nseries];
		int[] newcolors=new int[nseries];
		if(zValues1.length>maxxpts){
			newxValues=new float[nseries][zValues1.length];
			newyValues=new float[nseries][zValues1.length];
			newzValues=new float[nseries][zValues1.length];
			for(int i=0;i<(nseries-1);i++){
				newnpts[i]=npts[i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				for(int j=0;j<npts[i];j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
					newzValues[i][j]=zValues[i][j];
				}
			}
			maxxpts=xValues1.length;
			newnpts[nseries-1]=maxxpts;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			for(int j=0;j<maxxpts;j++){
				newxValues[nseries-1][j]=xValues1[j];
				newyValues[nseries-1][j]=yValues1[j];
				newzValues[nseries-1][j]=zValues1[j];
			}
		}else{
			newxValues=new float[nseries][maxxpts];
			newyValues=new float[nseries][maxxpts];
			newzValues=new float[nseries][maxxpts];
			for(int i=0;i<(nseries-1);i++){
				newnpts[i]=npts[i];
				newshapes[i]=shapes[i];
				newcolors[i]=colors[i];
				for(int j=0;j<maxxpts;j++){
					newxValues[i][j]=xValues[i][j];
					newyValues[i][j]=yValues[i][j];
					newzValues[i][j]=zValues[i][j];
				}
			}
			newnpts[nseries-1]=xValues1.length;
			newshapes[nseries-1]=0;
			newcolors[nseries-1]=nseries-1;
			for(int j=0;j<xValues1.length;j++){
				newxValues[nseries-1][j]=xValues1[j];
				newyValues[nseries-1][j]=yValues1[j];
				newzValues[nseries-1][j]=zValues1[j];
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
				logxmin=(float)Math.log(findmingt0(xValues,npts,xMax));
			}else{
				logxmin=(float)Math.log(xMin);
			}
			logxmax=(float)Math.log(xMax);
			logxscale=WIDTH/(logxmax-logxmin);
		}
		if(logy){
			if(yMin<=0.0f){
				logymin=(float)Math.log(findmingt0(yValues,npts,yMax));
			}else{
				logymin=(float)Math.log(yMin);
			}
			logymax=(float)Math.log(yMax);
			logyscale=WIDTH/(logymax-logymin);
		}
		if(logz){
			if(zMin<=0.0f){
				logzmin=(float)Math.log(findmingt0(zValues,npts,zMax));
			}else{
				logzmin=(float)Math.log(zMin);
			}
			logzmax=(float)Math.log(zMax);
			logzscale=HEIGHT/(logzmax-logzmin);
		}
		// IJ.showMessage("testdraw1");

		xScale=WIDTH/(xMax-xMin);
		yScale=WIDTH/(yMax-yMin);
		zScale=HEIGHT/(zMax-zMin);

		// jr.setBufferSize(width,height);
		drawAxisLabels(jr);
		if(!logx&&!logy&&!logz){
			for(int j=0;j<nseries;j++){
				Color tempcolor=getColor(colors[j]);
				int xpoints[]=new int[npts[j]];
				int ypoints[]=new int[npts[j]];
				int zpoints[]=new int[npts[j]];
				for(int i=0;i<npts[j];i++){
					xpoints[i]=LEFT_MARGIN+(int)((xValues[j][i]-xMin)*xScale);
					if(xpoints[i]<LEFT_MARGIN){
						xpoints[i]=LEFT_MARGIN;
					}
					if(xpoints[i]>LEFT_MARGIN+WIDTH){
						xpoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[j];i++){
					ypoints[i]=LEFT_MARGIN+(int)((yValues[j][i]-yMin)*yScale);
					if(ypoints[i]<LEFT_MARGIN){
						ypoints[i]=LEFT_MARGIN;
					}
					if(ypoints[i]>LEFT_MARGIN+WIDTH){
						ypoints[i]=LEFT_MARGIN+WIDTH;
					}
				}
				for(int i=0;i<npts[j];i++){
					zpoints[i]=TOP_MARGIN+HEIGHT-(int)((zValues[j][i]-zMin)*zScale);
					if(zpoints[i]<TOP_MARGIN){
						zpoints[i]=TOP_MARGIN;
					}
					if(zpoints[i]>TOP_MARGIN+HEIGHT){
						zpoints[i]=TOP_MARGIN+HEIGHT;
					}
				}
				if(j!=selected){
					if(shapes[j]==0){
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}else{
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[j],shapes[j],tempcolor);
						if(drawAnnotations && annotations!=null) drawPolytext(jr,xpoints,ypoints,zpoints,npts[j],annotations[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}
				}else{
					if(shapes[j]==0){
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[j],1,tempcolor);
						if(drawAnnotations && annotations!=null) drawPolytext(jr,xpoints,ypoints,zpoints,npts[j],annotations[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}else{
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}
				}
			}
		}else{
			for(int j=0;j<nseries;j++){
				Color tempcolor=getColor(colors[j]);
				int xpoints[]=new int[npts[j]];
				int ypoints[]=new int[npts[j]];
				int zpoints[]=new int[npts[j]];
				for(int i=0;i<npts[j];i++){
					if(logx){
						float xtemp;
						if(xValues[j][i]>0.0f){
							xtemp=(float)Math.log(xValues[j][i]);
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
				for(int i=0;i<npts[j];i++){
					if(logy){
						float ytemp;
						if(yValues[j][i]>0.0f){
							ytemp=(float)Math.log(yValues[j][i]);
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
				for(int i=0;i<npts[j];i++){
					if(logz){
						float ztemp;
						if(zValues[j][i]>0.0f){
							ztemp=(float)Math.log(zValues[j][i]);
						}else{
							ztemp=logzmin;
						}
						zpoints[i]=TOP_MARGIN+HEIGHT-(int)((ztemp-logzmin)*logzscale);
					}else{
						zpoints[i]=TOP_MARGIN+HEIGHT-(int)((zValues[j][i]-zMin)*zScale);
					}
					if(zpoints[i]<TOP_MARGIN){
						zpoints[i]=TOP_MARGIN;
					}
					if(zpoints[i]>TOP_MARGIN+HEIGHT){
						zpoints[i]=TOP_MARGIN+HEIGHT;
					}
				}
				if(j!=selected){
					if(shapes[j]==0){
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}else{
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[j],shapes[j],tempcolor);
						if(drawAnnotations && annotations!=null) drawPolytext(jr,xpoints,ypoints,zpoints,npts[j],annotations[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}
				}else{
					if(shapes[j]==0){
						drawPolyshape(jr,xpoints,ypoints,zpoints,npts[j],1,tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
					}else{
						drawPolyline(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
						if(drawAnnotations && annotations!=null) drawPolytext(jr,xpoints,ypoints,zpoints,npts[j],annotations[j],tempcolor);
						if(drawIndex) drawPolyindex(jr,xpoints,ypoints,zpoints,npts[j],tempcolor);
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
			if(logz){
				float tempz=logzmin+(i/3.0f)*(logzmax-logzmin);
				zticklabels[i]=(float)Math.exp(tempz);
			}else{
				zticklabels[i]=zMin+(i/3.0f)*(zMax-zMin);
			}
		}

		// draw the z axis labels
		String s=jutils.formatted_string(zticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+HEIGHT-5,Color.BLACK);
		s=jutils.formatted_string(zticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+(2*HEIGHT)/3-5,Color.BLACK);
		s=jutils.formatted_string(zticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN+HEIGHT/3-5,Color.BLACK);
		s=jutils.formatted_string(zticklabels[3]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+10,LEFT_MARGIN,TOP_MARGIN-5,Color.BLACK);

		jr.addText3D(zLabel,LEFT_MARGIN+WIDTH+25,LEFT_MARGIN,TOP_MARGIN+HEIGHT/2-15,Color.BLACK);

		// now the x axis labels
		s=jutils.formatted_string(xticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(xticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH/3-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(xticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+2*WIDTH/3-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(xticklabels[3]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH-10,LEFT_MARGIN+WIDTH+40,TOP_MARGIN+HEIGHT,Color.BLACK);

		jr.addText3D(xLabel,LEFT_MARGIN+WIDTH/2,LEFT_MARGIN+WIDTH+60,TOP_MARGIN+HEIGHT,Color.BLACK);

		// and the y axis labels
		s=jutils.formatted_string(yticklabels[0]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(yticklabels[1]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+WIDTH/3+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(yticklabels[2]);
		jr.addText3D(s,LEFT_MARGIN+WIDTH+20,LEFT_MARGIN+2*WIDTH/3+10,TOP_MARGIN+HEIGHT,Color.BLACK);
		s=jutils.formatted_string(yticklabels[3]);
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

	private void drawPolyline(renderer jr,int[] xpoints,int[] ypoints,int[] zpoints,int npts,Color color){
		for(int i=1;i<npts;i++){
			jr.addLine3D(xpoints[i-1],ypoints[i-1],zpoints[i-1],xpoints[i],ypoints[i],zpoints[i],color);
			jr.elements[jr.elements.length-1].thick=thick;
		}
	}

	private void drawPolyshape(renderer jr,int[] xpoints,int[] ypoints,int[] zpoints,int npts,int shape,Color color){
		for(int i=0;i<npts;i++){
			jr.addPoint3D(xpoints[i],ypoints[i],zpoints[i],shape-1,color);
		}
	}
	
	private void drawPolytext(renderer jr,int[] xpoints,int[] ypoints,int[] zpoints,int npts,String text,Color color){
		for(int i=0;i<npts;i++){
				jr.addText3D(text,xpoints[i]+shapesize,ypoints[i],zpoints[i],color);
		}
	}
	
	private void drawPolyindex(renderer jr,int[] xpoints,int[] ypoints,int[] zpoints,int npts,Color color){
		for(int i=0;i<npts;i++){
				jr.addText3D(""+i,xpoints[i]+shapesize,ypoints[i],zpoints[i],color);
		}
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
		return new float[][][]{zValues};
	}

	public float[][] getZValues(int series){
		return new float[][]{zValues[series]};
	}

	public int[][] getNpts(){
		return new int[][]{npts};
	}
	
	public void setThick(boolean thick){this.thick=thick;}
	
	public void showIndex(boolean drawIndex){this.drawIndex=drawIndex;}
	
	/*public ColorProcessor getProcessor(){
		return new ColorProcessor(getImage());
	}

	public Image getImage(){
		renderer jr=new renderer(getWidth(),getHeight());
		drawPlot(jr);
		return jr.renderimage();
	}*/

	public void saveplot2os(OutputStream os){
		jdataio jdio=new jdataio();
		// start with unique identifier for a 3D plot
		jdio.writestring(os,"pw2_file_type");
		jdio.writeintelint(os,2);
		jdio.writestring(os,getxLabel());
		jdio.writestring(os,getyLabel());
		jdio.writestring(os,getzLabel());
		jdio.writeintelint(os,nseries); // number of series'
		jdio.writeintelint(os,getmaxxpts()); // max number of pts
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
			jdio.writeintelint(os,npts[l]); // number of points in this series
			jdio.writeintelint(os,shapes[l]); // shape index
			jdio.writeintelint(os,colors[l]); // color index
			jdio.writeintelfloatarray(os,xValues[l],npts[l]); // x values
			jdio.writeintelfloatarray(os,yValues[l],npts[l]); // y values
			jdio.writeintelfloatarray(os,zValues[l],npts[l]); // z values
		}
		// save the annotations if they exist
		if(annotations!=null && annotations.length==nseries){
			jdio.writeintelint(os,1);
			for(int i=0;i<nseries;i++) jdio.writestring(os,annotations[i]);
		}
	}

}
