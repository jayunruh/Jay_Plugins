/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Image;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferInt;
import java.awt.image.DirectColorModel;
import java.awt.image.PixelGrabber;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.OutputStream;
import java.util.Properties;

import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.svg.SVGGraphics2D;

public class SVGPlotter extends Plotter{
	public SVGGraphics2D vg2;
	public boolean antialias;

	public SVGPlotter(int width,int height,OutputStream os){
		this.width=width;
		this.height=height;
		shapesize=8;
		try{
			vg2=new SVGGraphics2D(os,new Dimension(width,height));
			vg2.setDeviceIndependent(true);
			vg2.startExport();
		}catch(Throwable e){
			IJ.log(e.toString());
		}
	}

	public SVGPlotter(int width,int height,String path){
		this.width=height;
		this.height=height;
		shapesize=8;
		try{
			File file=new File(path);
			vg2=new SVGGraphics2D(file,new Dimension(width,height));
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
		//Image img2=copyImage(img,0);
		vg2.drawImage(img,x,y,null);
	}
	
	public void setTransparency(Image img,int transparency){
		int[] pixels=new int[width*height];
		PixelGrabber pg=new PixelGrabber(img, 0,0,width,height,pixels,0,width);
		try{
			pg.grabPixels();
		}catch(InterruptedException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for(int i=0;i<pixels.length;i++){
			int[] rgb=jutils.intval2rgb(pixels[i]);
			int temp=(transparency<<24)|(rgb[0]<<16)|(rgb[1]<<8)|rgb[2];
			pixels[i]=temp;
		}
	}
	
	public Image copyImage(Image img,int transparency){
		int[] pixels=new int[width*height];
		PixelGrabber pg=new PixelGrabber(img, 0,0,width,height,pixels,0,width);
		try{
			pg.grabPixels();
		}catch(InterruptedException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for(int i=0;i<pixels.length;i++){
			int[] rgb=jutils.intval2rgb(pixels[i]);
			int temp=(transparency<<24)|(rgb[0]<<16)|(rgb[1]<<8)|rgb[2];
			pixels[i]=temp;
		}
		DataBuffer dataBuffer=new DataBufferInt(pixels,width*height,0);
		ColorModel cm=new DirectColorModel(32,0xff0000,0xff00,0xff,0xff000000);
		WritableRaster wr = cm.createCompatibleWritableRaster(1, 1);
		SampleModel rgbSampleModel=wr.getSampleModel();
		rgbSampleModel=rgbSampleModel.createCompatibleSampleModel(width,height);
		WritableRaster rgbRaster=Raster.createWritableRaster(rgbSampleModel,dataBuffer,null);
		Image img2=new BufferedImage(cm,rgbRaster,false,null);
		return img2;
	}

}
