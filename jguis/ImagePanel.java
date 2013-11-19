/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.gui.GenericDialog;
import ij.process.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class ImagePanel extends JPanel implements MouseListener,MouseMotionListener,ActionListener{

	public Image plotimage;
	public float[] imagedata;
	private boolean zooming;
	private int startx,starty,currx,curry,panelwidth,panelheight,imagewidth,imageheight;
	private int xMin,yMin,width,height,maxmag,currmag;
	private float intMin,intMax,aspectratio;

	public void init(int maxmag1,int currmag1,int imagewidth1,int imageheight1,float[] data){
		zooming=false;
		maxmag=maxmag1;
		currmag=currmag1;
		imagewidth=imagewidth1;
		imageheight=imageheight1;
		aspectratio=(float)imagewidth/(float)imageheight;
		panelwidth=imagewidth*maxmag;
		panelheight=imageheight*maxmag;
		imagedata=data;
		if(imagedata==null){
			imagedata=new float[imagewidth*imageheight];
		}
		addMouseMotionListener(this);
		addMouseListener(this);
		autoscale();
	}

	public void setBounds(int x,int y,int mag){
		panelwidth=imagewidth*mag;
		panelheight=imageheight*mag;
		super.setBounds(x,y,panelwidth,panelheight);
		repaint();
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(int xMin1,int yMin1,int Width1,int Height1,float intMin1,float intMax1){
		xMin=xMin1;
		yMin=yMin1;
		width=Width1;
		height=Height1;
		intMin=intMin1;
		intMax=intMax1;
		if(width>height){
			height=(int)(aspectratio/width);
		}else{
			width=(int)(aspectratio*height);
		}
		repaint();
	}

	public void autoscale(){
		float[] minmax=getminmax(imagedata);
		intMin=minmax[0];
		intMax=minmax[1];
		xMin=0;
		yMin=0;
		width=imagewidth;
		height=imageheight;
		repaint();
	}

	float[] getminmax(float[] array){
		float[] minmax={array[0],array[0]};
		for(int i=1;i<imagewidth*imageheight;i++){
			if(array[i]<minmax[0]){
				minmax[0]=array[i];
			}
			if(array[i]>minmax[1]){
				minmax[1]=array[i];
			}
		}
		return minmax;
	}

	public void changeImage(float[] newimage){
		imagedata=newimage;
		autoscale();
	}

	public void paint(Graphics g){
		// start by creating a cropped version of the imagedata
		float[] cropped=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if((xMin+j)<imagewidth&&(yMin+i)<imageheight){
					cropped[j+i*width]=imagedata[xMin+j+(yMin+i)*imagewidth];
				}
			}
		}
		// now scale the image
		FloatProcessor fp=new FloatProcessor(width,height,cropped,null);
		fp.setMinAndMax(intMin,intMax);
		ImageProcessor fpscaled=fp.resize(currmag*imagewidth,currmag*imageheight);
		// now get the image and draw it
		plotimage=fpscaled.createImage();
		g.setColor(Color.white);
		g.fillRect(0,0,panelwidth,panelheight);
		g.setColor(Color.black);
		g.clipRect(0,0,panelwidth,panelheight);
		g.drawImage(plotimage,0,0,this);
		if(zooming){
			// draw the zoom rectangle
			// make sure to maintain aspect ratio
			g.setColor(Color.RED);
			int tempwidth=Math.abs(currx-startx);
			int widthsign=(currx-startx)>0?1:-1;
			int tempheight=Math.abs(curry-starty);
			int heightsign=(curry-starty)>0?1:-1;
			if(tempwidth>tempheight*aspectratio){
				tempheight=(int)(heightsign*tempwidth/aspectratio);
				tempwidth*=widthsign;
			}else{
				tempwidth=(int)(widthsign*tempheight*aspectratio);
				tempheight*=heightsign;
			}
			int tempcurrx=startx+tempwidth;
			int tempcurry=starty+tempheight;
			g.drawLine(tempcurrx,tempcurry,startx,tempcurry);
			g.drawLine(startx,tempcurry,startx,starty);
			g.drawLine(startx,starty,tempcurrx,starty);
			g.drawLine(tempcurrx,starty,tempcurrx,tempcurry);
			g.setColor(Color.BLACK);
		}
	}

	public void update(Graphics g){
		paint(g);
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
		if(e.getButton()==MouseEvent.BUTTON3){
			handlepopup(e);
		}
		if(e.getButton()==MouseEvent.BUTTON1){
			startx=e.getX();
			starty=e.getY();
			currx=startx;
			curry=starty;
			zooming=true;
		}
	}

	public void mouseReleased(MouseEvent e){
		if(zooming){
			int endx=e.getX();
			int endy=e.getY();
			if(endx==startx||endy==starty){
				zooming=false;
				return;
			}

			int tempwidth=Math.abs(endx-startx);
			int widthsign=(endx-startx)>0?1:-1;
			int tempheight=Math.abs(endy-starty);
			int heightsign=(endy-starty)>0?1:-1;
			if(tempwidth>tempheight*aspectratio){
				tempheight=(int)(heightsign*tempwidth/aspectratio);
				tempwidth*=widthsign;
			}else{
				tempwidth=(int)(widthsign*tempheight*aspectratio);
				tempheight*=heightsign;
			}
			int tempcurrx=startx+tempwidth;
			int tempcurry=starty+tempheight;
			int rx=startx;
			if(tempcurrx<startx){
				rx=tempcurrx;
			}
			int ry=starty;
			if(tempcurry<starty){
				ry=tempcurry;
			}
			xMin=rx;
			yMin=ry;
			width=Math.abs(tempwidth);
			height=Math.abs(tempheight);

			zooming=false;
			repaint();
		}
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		if(zooming){
			currx=e.getX();
			curry=e.getY();
			repaint();
		}
	}

	public void actionPerformed(ActionEvent e){
		String cmd=e.getActionCommand();
		if(cmd=="Autoscale"){
			autoscale();
		}
		if(cmd=="Edit"){
			editPlot();
		}
	}

	void handlepopup(MouseEvent e){
		int x=e.getX();
		int y=e.getY();
		JPopupMenu popup=new JPopupMenu("");
		JMenuItem mi;
		mi=popup.add("Autoscale");
		mi.addActionListener(this);
		mi=popup.add("Edit");
		mi.addActionListener(this);
		add(popup);
		popup.show(this,x,y);
	}

	void editPlot(){
		GenericDialog gd=new GenericDialog("Image Options");
		gd.addNumericField("x min",xMin,0);
		gd.addNumericField("x max",xMin+width,0);
		gd.addNumericField("y min",yMin,0);
		gd.addNumericField("y max",yMin+height,0);
		gd.addNumericField("intensity min",intMin,5,10,null);
		gd.addNumericField("intensity max",intMax,5,10,null);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		xMin=(int)gd.getNextNumber();
		width=(int)gd.getNextNumber()-xMin;
		yMin=(int)gd.getNextNumber();
		height=(int)gd.getNextNumber()-yMin;
		intMin=(float)gd.getNextNumber();
		intMax=(float)gd.getNextNumber();
		if(width>height){
			height=(int)(aspectratio/width);
		}else{
			width=(int)(aspectratio*height);
		}
		repaint();
	}
}
