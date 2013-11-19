/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.ps.PSGraphics2D;

import ij.IJ;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.io.File;
import java.io.OutputStream;
import java.util.Properties;

public class PSPlotter extends Plotter{
	public PSGraphics2D vg2;
	public boolean antialias;

	public PSPlotter(int width,int height,OutputStream os){
		this.width=width;
		this.height=height;
		shapesize=8;
		try{
			vg2=new PSGraphics2D(os,new Dimension(width,height));
			vg2.setDeviceIndependent(true);
			vg2.startExport();
		}catch(Throwable e){
			IJ.log(e.toString());
		}
	}

	public PSPlotter(int width,int height,String path){
		this.width=height;
		this.height=height;
		shapesize=8;
		try{
			File file=new File(path);
			vg2=new PSGraphics2D(file,new Dimension(width,height));
			vg2.setDeviceIndependent(true);
			vg2.startExport();
		}catch(Throwable e){
			IJ.log(e.toString());
		}
	}

	public void endPlotting(){
		vg2.endExport();
	}

	public void clear_plot(){
		vg2.clearRect(0,0,width,height);
	}

	public void setColor(Color color){
		vg2.setColor(color);
	}

	public void setFont(Font font){
		vg2.setFont(font);
	}

	public void drawLine(int x1,int y1,int x2,int y2){
		vg2.drawLine(x1,y1,x2,y2);
	}
	
	public void drawPolyLine(int[] xcoords,int[] ycoords,boolean close){
		if(!close) vg2.drawPolyline(xcoords,ycoords,xcoords.length);
		else vg2.drawPolygon(xcoords,ycoords,xcoords.length);
	}

	public void drawString(String string,int x,int y){
		vg2.drawString(string,x,y);
	}

	public void drawVerticalString(String string,int x,int y){
		// x and y here are the center of the string
		AffineTransform at=vg2.getTransform();
		int halfwidth=getStringWidth(string)/2;
		FontMetrics fm=getFontMetrics();
		int halfheight=(fm.getAscent()+fm.getDescent())/2;
		vg2.rotate(-0.5*Math.PI,x,y);
		vg2.drawString(string,x-halfwidth,y+halfheight);
		vg2.setTransform(at);
	}

	public void fillRect(int x,int y,int width,int height){
		vg2.fillRect(x,y,width,height);
	}

	public void setAntiAliasedText(boolean antialias){
		Properties user=new Properties();
		user.setProperty(ImageGraphics2D.ANTIALIAS_TEXT,"true");
		vg2.setProperties(user);
		this.antialias=true;
	}

	public void setClipRect(int x,int y,int width,int height){
		vg2.setClip(x,y,width,height);
	}

	public int getStringWidth(String string){
		FontMetrics fm=vg2.getFontMetrics();
		java.awt.geom.Rectangle2D r=fm.getStringBounds(string,vg2);
		return (int)r.getWidth();
	}

	public Color getColor(){
		return vg2.getColor();
	}

	public Font getFont(){
		return vg2.getFont();
	}

	public boolean getAntiAliasedText(){
		return antialias;
	}

	public FontMetrics getFontMetrics(){
		return vg2.getFontMetrics();
	}

	public void unsetClip(){
		vg2.setClip(null);
	}

	public void drawImage(Image img,int x,int y){
		vg2.drawImage(img,x,y,null);
	}

}
