/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import java.awt.Graphics;
import java.awt.Image;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

public class Stack_Printer implements Printable{
	private ImagePlus imp;
	private int nslices;
	private boolean rotate;

	public Stack_Printer(boolean rotate){
		this.rotate=rotate;
	}

	public void print(ImagePlus imp){
		this.imp=imp;
		nslices=imp.getStackSize();
		PrinterJob pj=PrinterJob.getPrinterJob();
		pj.setPrintable(this);
		// pj.pageDialog(pj.defaultPage());
		if(pj.printDialog()){
			imp.startTiming();
			try{
				pj.print();
			}catch(PrinterException e){
				IJ.log(""+e);
			}
		}
	}

	public int print(Graphics g,PageFormat pf,int pageIndex){
		if(pageIndex>=nslices)
			return NO_SUCH_PAGE;
		ImageProcessor ip=imp.getStack().getProcessor(pageIndex+1);
		if(rotate)
			ip=ip.rotateLeft();
		// new ImagePlus("ip", ip.duplicate()).show();
		int width=ip.getWidth();
		int height=ip.getHeight();
		int margin=0;
		double scale=1.0;
		int dstWidth=(int)(width*scale);
		int dstHeight=(int)(height*scale);
		int pageX=(int)pf.getImageableX();
		int pageY=(int)pf.getImageableY();
		int dstX=pageX+margin;
		int dstY=pageY+margin;
		Image img=ip.createImage();
		double pageWidth=pf.getImageableWidth()-2*margin;
		double pageHeight=pf.getImageableHeight()-2*margin;
		if(dstWidth>pageWidth||dstHeight>pageHeight){
			// scale to fit page
			double hscale=pageWidth/dstWidth;
			double vscale=pageHeight/dstHeight;
			double scale2=hscale<=vscale?hscale:vscale;
			dstWidth=(int)(dstWidth*scale2);
			dstHeight=(int)(dstHeight*scale2);
		}
		g.drawImage(img,dstX,dstY,dstX+dstWidth,dstY+dstHeight,0,0,width,height,null);
		return PAGE_EXISTS;
	}

}
