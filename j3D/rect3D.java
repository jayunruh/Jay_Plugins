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

public class rect3D extends element3D implements Cloneable{
	public point3D[] pt;
	public int xsize,ysize,zsize;

	// this is really a cube with rectangular sides
	public rect3D(int x,int y,int z,int xsize1,int ysize1,int zsize1,Color color1){
		pt=new point3D[8];
		xsize=xsize1;
		ysize=ysize1;
		zsize=zsize1;
		int halfxsize=xsize/2;
		int halfysize=ysize/2;
		int halfzsize=zsize/2;
		// start with the upper square
		pt[0]=new point3D(x-halfxsize,y-halfysize,z-halfzsize);
		pt[1]=new point3D(x+halfxsize,y-halfysize,z-halfzsize);
		pt[2]=new point3D(x+halfxsize,y+halfysize,z-halfzsize);
		pt[3]=new point3D(x-halfxsize,y+halfysize,z-halfzsize);
		pt[4]=new point3D(x-halfxsize,y-halfysize,z+halfzsize);
		pt[5]=new point3D(x+halfxsize,y-halfysize,z+halfzsize);
		pt[6]=new point3D(x+halfxsize,y+halfysize,z+halfzsize);
		pt[7]=new point3D(x-halfxsize,y+halfysize,z+halfzsize);
		color=color1;
	}
	
	public rect3D(point3D[] pt,int xsize,int ysize,int zsize,Color color){
		this.pt=pt;
		this.xsize=xsize;
		this.ysize=ysize;
		this.zsize=zsize;
		this.color=color;
	}

	public void moveto(int x,int y,int z){
		int halfxsize=xsize/2;
		int halfysize=ysize/2;
		int halfzsize=zsize/2;
		pt[0].moveto(x-halfxsize,y-halfysize,z-halfzsize);
		pt[1].moveto(x+halfxsize,y-halfysize,z-halfzsize);
		pt[2].moveto(x+halfxsize,y+halfysize,z-halfzsize);
		pt[3].moveto(x-halfxsize,y+halfysize,z-halfzsize);
		pt[4].moveto(x-halfxsize,y-halfysize,z+halfzsize);
		pt[5].moveto(x+halfxsize,y-halfysize,z+halfzsize);
		pt[6].moveto(x+halfxsize,y+halfysize,z+halfzsize);
		pt[7].moveto(x-halfxsize,y+halfysize,z+halfzsize);
	}

	public void translate(int transx,int transy,int transz){
		for(int i=0;i<8;i++){
			pt[i].translate(transx,transy,transz);
		}
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1){
		for(int i=0;i<pt.length;i++){
			pt[i].set_rot_about_vector(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
		}
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<8;i++){
			pt[i].transform_negative_perspective(horizon_dist,centerx,centery,centerz);
		}
	}
	
	public void transform(double[][] transmat,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].transform(transmat,centerx,centery,centerz);
		}
	}

	public int getzpos(){
		float sum=0.0f;
		for(int i=0;i<8;i++){
			sum+=pt[i].rz;
		}
		int zpos=(int)(0.125f*sum);
		return zpos;
	}

	public void drawelement(Graphics g){
		Color tempcolor=g.getColor();
		g.setColor(color);
		// draw the top square
		j3Dutils.drawLine(g,pt[0].rx,pt[0].ry,pt[1].rx,pt[1].ry,thick);
		j3Dutils.drawLine(g,pt[1].rx,pt[1].ry,pt[2].rx,pt[2].ry,thick);
		j3Dutils.drawLine(g,pt[2].rx,pt[2].ry,pt[3].rx,pt[3].ry,thick);
		j3Dutils.drawLine(g,pt[3].rx,pt[3].ry,pt[0].rx,pt[0].ry,thick);
		// draw the bottom square
		j3Dutils.drawLine(g,pt[4].rx,pt[4].ry,pt[5].rx,pt[5].ry,thick);
		j3Dutils.drawLine(g,pt[5].rx,pt[5].ry,pt[6].rx,pt[6].ry,thick);
		j3Dutils.drawLine(g,pt[6].rx,pt[6].ry,pt[7].rx,pt[7].ry,thick);
		j3Dutils.drawLine(g,pt[7].rx,pt[7].ry,pt[4].rx,pt[4].ry,thick);
		// draw the vertical lines
		j3Dutils.drawLine(g,pt[0].rx,pt[0].ry,pt[4].rx,pt[4].ry,thick);
		j3Dutils.drawLine(g,pt[1].rx,pt[1].ry,pt[5].rx,pt[5].ry,thick);
		j3Dutils.drawLine(g,pt[2].rx,pt[2].ry,pt[6].rx,pt[6].ry,thick);
		j3Dutils.drawLine(g,pt[3].rx,pt[3].ry,pt[7].rx,pt[7].ry,thick);
		g.setColor(tempcolor);
	}
	
	public rect3D clone(){
		point3D[] temp=new point3D[pt.length];
		for(int i=0;i<temp.length;i++) temp[i]=pt[i].clone();
		return new rect3D(temp,xsize,ysize,zsize,color);
	}

}
