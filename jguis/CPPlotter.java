/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.gui.Roi;
import ij.process.Blitter;
import ij.process.ColorProcessor;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Image;
import java.awt.Rectangle;

public class CPPlotter extends Plotter{
	public ColorProcessor cp;
	public boolean antialias;
	public Color color;
	public Rectangle clip;

	public CPPlotter(int width,int height){
		this.width=width;
		this.height=height;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		shapesize=8;
	}

	public void clear_plot(){
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
	}

	public void setColor(Color color){
		cp.setColor(color);
		this.color=color;
	}

	public void drawLine(int x1,int y1,int x2,int y2){
		cp.drawLine(x1,y1,x2,y2);
	}
	
	public void drawPolyLine(int[] xcoords,int[] ycoords,boolean close){
		for(int i=0;i<(xcoords.length-1);i++){
			drawLine(xcoords[i],ycoords[i],xcoords[i+1],ycoords[i+1]);
		}
		if(close){
			drawLine(xcoords[xcoords.length-1],ycoords[xcoords.length-1],xcoords[0],ycoords[0]);
		}
	}

	public void drawString(String string,int x,int y){
		cp.drawString(string,x,y);
	}

	public void drawVerticalString(String string,int x,int y){
		// x and y here are the center of the string
		FontMetrics fm=cp.getFontMetrics();
		int ascent=fm.getAscent();
		int descent=fm.getDescent();
		int height=ascent+descent;
		int width=cp.getStringWidth(string);
		int[] mask=new int[width*height];
		for(int i=0;i<width*height;i++){
			mask[i]=0xffffffff;
		}
		ColorProcessor cp2=new ColorProcessor(width,height,mask);
		cp2.setAntialiasedText(antialias);
		cp2.setFont(cp.getFont());
		cp2.drawString(string,0,ascent);
		ColorProcessor cp3=(ColorProcessor)cp2.rotateLeft();
		mask=(int[])cp3.getPixels();
		int newx=x-height/2;
		int newy=y-width/2;
		int[] pixels=(int[])cp.getPixels();
		int owidth=cp.getWidth();
		int oheight=cp.getHeight();
		for(int i=0;i<width;i++){
			int b=newy+i;
			if(b>=0&&b<oheight){
				for(int j=0;j<height;j++){
					int a=newx+j;
					if(a>=0&&a<owidth){
						if(mask[j+i*height]!=0xffffffff){
							pixels[a+owidth*b]=mask[j+i*height];
						}
					}
				}
			}
		}
	}

	public void fillRect(int x,int y,int width,int height){
		// have the implement the clip rect ourselves here
		int x1=x;
		int y1=y;
		int x2=x+width;
		int y2=y+height;
		if(x1<clip.x)
			x1=clip.x;
		if(x2>(clip.x+clip.width))
			x2=clip.x+clip.width;
		if(x2<x1)
			return;
		if(y1<clip.y)
			y1=clip.y;
		if(y2>(clip.y+clip.width))
			y2=clip.y+clip.width;
		if(y2<y1)
			return;
		Roi fillroi=new Roi(x1,y1,x2-x1,y2-y1);
		cp.fill(fillroi);
	}

	public void setFont(Font font){
		cp.setFont(font);
	}

	public void setAntiAliasedText(boolean antialias){
		cp.setAntialiasedText(antialias);
		this.antialias=true;
	}

	public void setClipRect(int x,int y,int width,int height){
		clip=new Rectangle(x,y,width,height);
		cp.setClipRect(clip);
	}

	public ColorProcessor get_output(){
		return cp;
	}

	public int getStringWidth(String string){
		return cp.getStringWidth(string);
	}

	public Color getColor(){
		return color;
	}

	public Font getFont(){
		return cp.getFont();
	}

	public boolean getAntiAliasedText(){
		return antialias;
	}

	public FontMetrics getFontMetrics(){
		return cp.getFontMetrics();
	}

	public void unsetClip(){
		cp.setClipRect(null);
		clip=new Rectangle(width,height);
	}

	public void endPlotting(){
	}

	public void drawImage(Image img,int x,int y){
		ColorProcessor cp2=new ColorProcessor(img);
		cp.copyBits(cp2,x,y,Blitter.COPY);
	}

}
