/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Image;

public abstract class Plotter{
	// this is an abstract superclass for different plotters
	public int shapesize;
	public int currx,curry;
	public static final String[] shape_names={"line","square","+","x","triangle"};
	public int width,height;

	public abstract void setColor(Color color);

	public abstract Color getColor();

	public abstract void setFont(Font font);

	public abstract Font getFont();

	public abstract void drawLine(int x1,int y1,int x2,int y2);
	
	public abstract void drawPolyLine(int[] xcoords,int[] ycoords,boolean close);

	public abstract void drawString(String string,int x,int y);

	public abstract void drawVerticalString(String string,int x,int y);

	public abstract void fillRect(int x,int y,int width,int height);

	public abstract void clear_plot();

	public abstract void setClipRect(int x,int y,int width,int height);

	public abstract void unsetClip();

	public abstract void setAntiAliasedText(boolean aliased);

	public abstract boolean getAntiAliasedText();

	public abstract int getStringWidth(String string);

	public abstract FontMetrics getFontMetrics();

	public abstract void drawImage(Image img,int x,int y);

	public abstract void endPlotting();

	public void setFontSize(int fontsize){
		Font temp=getFont();
		Font temp2=new Font(temp.getName(),temp.getStyle(),fontsize);
		setFont(temp2);
	}

	public int getFontSize(){
		return getFont().getSize();
	}

	public void drawRect(int x,int y,int width,int height){
		moveTo(x,y);
		lineTo(x+width,y);
		lineTo(x+width,y+height);
		lineTo(x,y+height);
		lineTo(x,y);
	}

	public void fillBoundedRect(int x,int y,int width,int height,Color rectcolor){
		Color temp=getColor();
		fillRect(x,y,width,height);
		setColor(rectcolor);
		drawRect(x,y,width,height);
		setColor(temp);
	}

	public void drawPolyColumns(int[] heights,int xstart,int ystart,float binwidth,int npts,Color rectcolor){
		// here the columns touch each other
		for(int i=0;i<heights.length;i++){
			int ystart2=ystart-heights[i];
			int xstart2=(int)((float)i*binwidth+(float)xstart);
			int xend2=(int)((float)i*binwidth+binwidth+(float)xstart);
			fillBoundedRect(xstart2,ystart2,xend2-xstart2,heights[i],rectcolor);
		}
	}

	public void drawPolyColumns(int[] x,int[] y,int ystart,int npts,Color rectcolor){
		// here the columns are a fixed width (can overlap)
		// columns are centered on x positions
		int hw=shapesize/2;
		for(int i=0;i<npts;i++){
			fillBoundedRect(x[i]-hw,y[i],shapesize,ystart-y[i],rectcolor);
		}
	}

	public void lineTo(int x,int y){
		drawLine(currx,curry,x,y);
		currx=x;
		curry=y;
	}

	public void moveTo(int x,int y){
		currx=x;
		curry=y;
	}

	public void drawPolyline(int[] x,int[] y,int n){
		//moveTo(x[0],y[0]);
		//for(int i=0;i<n;i++)
		//	lineTo(x[i],y[i]);
		drawPolyLine(x,y,false);
	}

	public void drawPolyshape(int[] x,int[] y,int shape1,int n){
		moveTo(x[0],y[0]);
		for(int i=0;i<n;i++){
			if(shape1==1){
				drawSquare(x[i],y[i]);
			}
			if(shape1==2){
				drawPlus(x[i],y[i]);
			}
			if(shape1==3){
				drawX(x[i],y[i]);
			}
			if(shape1==4){
				drawTriangle(x[i],y[i]);
			}
		}
	}

	public void drawPolyerrors(int[] x,int[] yu,int[] yl,int n){
		for(int i=0;i<n;i++){
			if(yu[i]!=yl[i]){
				moveTo(x[i]-shapesize/2,yl[i]);
				lineTo(x[i]+shapesize/2,yl[i]);
				moveTo(x[i],yl[i]);
				lineTo(x[i],yu[i]);
				moveTo(x[i]-shapesize/2,yu[i]);
				lineTo(x[i]+shapesize/2,yu[i]);
			}
		}
	}

	public void drawSquare(int x,int y){
		moveTo(x-shapesize/2,y-shapesize/2);
		lineTo(x+shapesize/2,y-shapesize/2);
		lineTo(x+shapesize/2,y+shapesize/2);
		lineTo(x-shapesize/2,y+shapesize/2);
		lineTo(x-shapesize/2,y-shapesize/2);
	}

	public void drawPlus(int x,int y){
		moveTo(x-shapesize/2,y);
		lineTo(x+shapesize/2,y);
		moveTo(x,y-shapesize/2);
		lineTo(x,y+shapesize/2);
	}

	public void drawX(int x,int y){
		moveTo(x-shapesize/2,y-shapesize/2);
		lineTo(x+shapesize/2,y+shapesize/2);
		moveTo(x-shapesize/2,y+shapesize/2);
		lineTo(x+shapesize/2,y-shapesize/2);
	}

	public void drawTriangle(int x,int y){
		moveTo(x,y-shapesize/2);
		lineTo(x-shapesize/2,y+shapesize/2);
		lineTo(x+shapesize/2,y+shapesize/2);
		lineTo(x,y-shapesize/2);
	}

}
