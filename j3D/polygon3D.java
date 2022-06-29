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
import java.awt.Polygon;

public class polygon3D extends element3D implements Cloneable{
	public point3D[] pt;
	public boolean closed;

	// here we have a 3D line
	public polygon3D(int[] xpts,int[] ypts,int[] zpts,Color color1){
		pt=new point3D[xpts.length];
		color=color1;
		for(int i=0;i<xpts.length;i++){
			pt[i]=new point3D(xpts[i],ypts[i],zpts[i]);
		}
		color=color1;
		closed=true;
	}

	public polygon3D(Polygon poly,Color color1){
		pt=new point3D[poly.xpoints.length];
		color=color1;
		for(int i=0;i<poly.xpoints.length;i++){
			pt[i]=new point3D(poly.xpoints[i],poly.ypoints[i],0);
		}
		color=color1;
		closed=true;
	}
	
	public polygon3D(point3D[] pt,boolean closed,Color color){
		this.pt=pt;
		this.closed=closed;
		this.color=color;
	}

	public Polygon getPolygon(){
		int[] xpoints=new int[pt.length];
		int[] ypoints=new int[pt.length];
		for(int i=0;i<pt.length;i++){
			xpoints[i]=pt[i].rx;
			ypoints[i]=pt[i].ry;
		}
		return new Polygon(xpoints,ypoints,pt.length);
	}

	public void moveto(int x,int y,int z){
		translate(x-pt[0].rx,y-pt[0].ry,z-pt[0].rz);
	}

	public void translate(int transx,int transy,int transz){
		for(int i=0;i<pt.length;i++){
			pt[i].translate(transx,transy,transz);
		}
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1){
		for(int i=0;i<pt.length;i++){
			pt[i].set_rot_about_vector(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
		}
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
			pt[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<pt.length;i++){
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
		for(int i=0;i<pt.length;i++){
			sum+=pt[i].rz;
		}
		int zpos=(int)(sum/pt.length);
		return zpos;
	}

	public void drawelement(Graphics g){
		Color tempcolor=g.getColor();
		g.setColor(color);
		for(int i=1;i<pt.length;i++){
			j3Dutils.drawLine(g,pt[i-1].rx,pt[i-1].ry,pt[i].rx,pt[i].ry,thick);
		}
		if(closed){
			j3Dutils.drawLine(g,pt[pt.length-1].rx,pt[pt.length-1].ry,pt[0].rx,pt[0].ry,thick);
		}
		g.setColor(tempcolor);
	}
	
	public polygon3D clone(){
		point3D[] temp=new point3D[pt.length];
		for(int i=0;i<pt.length;i++) temp[i]=pt[i].clone();
		return new polygon3D(temp,closed,color);
	}

}
