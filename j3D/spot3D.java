/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import java.awt.Color;
import java.awt.Graphics;

public class spot3D extends element3D implements Cloneable{
	public point3D point;
	public int shape; // shapes are square,plus,x,triangle
	public int shapesize=8;
	public int transshapesize;

	public spot3D(int x,int y,int z,int shape1,Color color1){
		point=new point3D(x,y,z);
		shape=shape1;
		color=color1;
		transshapesize=shapesize;
	}
	
	public spot3D(point3D pt,int shape1,Color color1){
		point=pt;
		shape=shape1;
		color=color1;
		transshapesize=shapesize;
	}

	public void moveto(int ptx,int pty,int ptz){
		point.moveto(ptx,pty,ptz);
	}

	public void translate(int transx,int transy,int transz){
		point.translate(transx,transy,transz);
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		point.rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		point.rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		point.rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		point.addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
	}
	
	public void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1){
		point.set_rot_about_vector(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		point.addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		point.addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		point.transform_perspective(horizon_dist,centerx,centery,centerz);
		if(horizon_dist<=0.0){
			transshapesize=shapesize;
		}else{
			double tempz=point.z-centerz;
			double temphordist=(tempz+horizon_dist)/horizon_dist;
			if(temphordist<=0){
				transshapesize=0;
			}else{
				transshapesize=(int)(shapesize*temphordist);
			}
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		point.transform_negative_perspective(horizon_dist,centerx,centery,centerz);
		if(horizon_dist<=0.0){
			transshapesize=shapesize;
		}else{
			double tempz=centerz-point.z;
			double temphordist=(horizon_dist-tempz)/horizon_dist;
			if(temphordist<=0){
				transshapesize=0;
			}else{
				transshapesize=(int)(shapesize*temphordist);
			}
		}
	}
	
	public void transform(double[][] transmat,int centerx,int centery,int centerz){
		point.transform(transmat,centerx,centery,centerz);
	}

	public int getzpos(){
		return point.rz;
	}

	public void drawelement(Graphics g){
		Color tempcolor=g.getColor();
		g.setColor(color);
		if(shape==0){
			drawSquare(g,point.rx,point.ry,transshapesize);
		}
		if(shape==1){
			drawPlus(g,point.rx,point.ry,transshapesize);
		}
		if(shape==2){
			drawX(g,point.rx,point.ry,transshapesize);
		}
		if(shape==3){
			drawTriangle(g,point.rx,point.ry,transshapesize);
		}
		if(shape==4){
			drawCircle(g,point.rx,point.ry,transshapesize);
		}
		if(shape==5){
			drawFilledSquare(g,point.rx,point.ry,transshapesize);
		}
		if(shape==6){
			drawFilledCircle(g,point.rx,point.ry,transshapesize);
		}
		g.setColor(tempcolor);
	}

	void drawSquare(Graphics g,int x,int y,int size){
		g.drawLine(x-size/2,y-size/2,x+size/2,y-size/2);
		g.drawLine(x+size/2,y-size/2,x+size/2,y+size/2);
		g.drawLine(x+size/2,y+size/2,x-size/2,y+size/2);
		g.drawLine(x-size/2,y+size/2,x-size/2,y-size/2);
	}

	void drawCircle(Graphics g,int x,int y,int size){
		g.drawOval(x-size/2,y-size/2,size,size);
	}

	void drawFilledCircle(Graphics g,int x,int y,int size){
		g.fillOval(x-size/2,y-size/2,size,size);
	}

	void drawFilledSquare(Graphics g,int x,int y,int size){
		g.fillRect(x-size/2,y-size/2,size,size);
	}

	void drawPlus(Graphics g,int x,int y,int size){
		g.drawLine(x-size/2,y,x+size/2,y);
		g.drawLine(x,y-size/2,x,y+size/2);
	}

	void drawX(Graphics g,int x,int y,int size){
		g.drawLine(x-size/2,y-size/2,x+size/2,y+size/2);
		g.drawLine(x-size/2,y+size/2,x+size/2,y-size/2);
	}

	void drawTriangle(Graphics g,int x,int y,int size){
		g.drawLine(x,y-size/2,x-size/2,y+size/2);
		g.drawLine(x-size/2,y+size/2,x+size/2,y+size/2);
		g.drawLine(x+size/2,y+size/2,x,y-size/2);
	}
	
	public spot3D clone(){
		spot3D spot=new spot3D(point.clone(),shape,color);
		spot.shapesize=shapesize;
		spot.transshapesize=transshapesize;
		return spot;
	}

}
