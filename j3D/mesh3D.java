/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import jalgs.interpolation;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Polygon;

public class mesh3D extends element3D implements Cloneable{
	//this is essentially a collection of 3D planar polygons that are connected to make a meshwork
	public polygon3D[] polys;
	public line3D[] connections;
	
	/*public static mesh3D make_sphere(int[] center,Color color,int rad){
		for(int i=center[2]-rad;i<center[2]+rad;i++){
			float rad2=0.0f;
		}
	}*/

	public mesh3D(Polygon[] poly,int[] zpos,Color color1){
		color=color1;
		polys=new polygon3D[poly.length];
		for(int i=0;i<polys.length;i++) polys[i]=poly2poly3D(poly[i],zpos[i],color);
		//each neighboring polygon is connected
		//these connections are lines from points in the bigger polygon to closest points in the smaller polygon
		//always want connections from polygon with more points to polygon with less
		line3D[][] conns=new line3D[polys.length-1][];
		int nconns=0;
		for(int i=0;i<(polys.length-1);i++){
			polygon3D src=polys[i];
			polygon3D dest=polys[i+1];
			if(dest.pt.length>src.pt.length){
				src=polys[i+1];
				dest=polys[i];
			}
			conns[i]=get_connections(src,dest);
			nconns+=conns[i].length;
		}
		connections=new line3D[nconns];
		int counter=0;
		for(int i=0;i<conns.length;i++){
			for(int j=0;j<conns[i].length;j++){
				connections[counter]=conns[i][j];
				counter++;
			}
		}
		conns=null;
	}
	
	public mesh3D(polygon3D[] polys,line3D[] connections,Color color){
		this.color=color;
		this.polys=polys;
		this.connections=connections;
	}
	
	public polygon3D poly2poly3D(Polygon poly,int z,Color color){
		int[] zpts=new int[poly.npoints];
		for(int i=0;i<poly.npoints;i++){
			zpts[i]=z;
		}
		return new polygon3D(poly.xpoints,poly.ypoints,zpts,color);
	}
	
	public line3D[] get_connections(polygon3D src,polygon3D dest){
		float[][] fpdest=new float[3][dest.pt.length];
		for(int i=0;i<dest.pt.length;i++){
			fpdest[0][i]=dest.pt[i].rx;
			fpdest[1][i]=dest.pt[i].ry;
			fpdest[2][i]=dest.pt[i].rz;
		}
		line3D[] conns=new line3D[src.pt.length];
		for(int i=0;i<src.pt.length;i++){
			float[] srcpt={src.pt[i].rx,src.pt[i].rx,src.pt[i].rz};
			float[] closest=interpolation.get_closest_index_3D(fpdest,srcpt);
			conns[i]=new line3D(src.pt[i],dest.pt[(int)closest[0]],color);
		}
		return conns;
	}

	public void moveto(int x,int y,int z){
		translate(x-polys[0].pt[0].rx,y-polys[0].pt[0].ry,z-polys[0].pt[0].rz);
	}

	public void translate(int transx,int transy,int transz){
		for(int i=0;i<polys.length;i++){
			polys[i].translate(transx,transy,transz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].translate(transx,transy,transz);
		}
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1){
		for(int i=0;i<polys.length;i++){
			polys[i].rotatecossinabout(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].rotatecossinabout(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
		}
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].transform_negative_perspective(horizon_dist,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].transform_negative_perspective(horizon_dist,centerx,centery,centerz);
		}
	}
	
	public void transform(double[][] transmat,int centerx,int centery,int centerz){
		for(int i=0;i<polys.length;i++){
			polys[i].transform(transmat,centerx,centery,centerz);
		}
		for(int i=0;i<connections.length;i++){
			connections[i].transform(transmat,centerx,centery,centerz);
		}
	}

	public int getzpos(){
		float sum=0.0f;
		for(int i=0;i<polys.length;i++){
			sum+=polys[i].getzpos();
		}
		int zpos=(int)(sum/polys.length);
		return zpos;
	}

	public void drawelement(Graphics g){
		Color tempcolor=g.getColor();
		g.setColor(color);
		//start by drawing all of the polygons
		for(int i=0;i<polys.length;i++) polys[i].drawelement(g);
		//now draw the connecting lines
		/*for(int i=0;i<connections.length;i++){
			connections[i].drawelement(g);
		}*/
		g.setColor(tempcolor);
	}
	
	public mesh3D clone(){
		polygon3D[] temp=new polygon3D[polys.length];
		for(int i=0;i<polys.length;i++) temp[i]=polys[i].clone();
		line3D[] tcons=new line3D[connections.length];
		for(int i=0;i<tcons.length;i++) tcons[i]=connections[i].clone();
		return new mesh3D(temp,tcons,color);
	}

}
