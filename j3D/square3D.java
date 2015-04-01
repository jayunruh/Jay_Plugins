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

public class square3D extends element3D implements Cloneable{
	public point3D[] pt;
	public final int size;
	public final String[] orientations={"xy","xz","yz"};
	public int orientation;

	// here we have a 3D square
	public square3D(int x,int y,int z,int size1,int orientation1,Color color1){
		pt=new point3D[4];
		size=size1;
		orientation=orientation1;
		int halfsize=size/2;
		if(orientation==0){
			pt[0]=new point3D(x-halfsize,y-halfsize,z);
			pt[1]=new point3D(x+halfsize,y-halfsize,z);
			pt[2]=new point3D(x+halfsize,y+halfsize,z);
			pt[3]=new point3D(x-halfsize,y+halfsize,z);
		}else{
			if(orientation==1){
				pt[0]=new point3D(x-halfsize,y,z-halfsize);
				pt[1]=new point3D(x+halfsize,y,z-halfsize);
				pt[2]=new point3D(x+halfsize,y,z+halfsize);
				pt[3]=new point3D(x-halfsize,y,z+halfsize);
			}else{
				pt[0]=new point3D(x,y-halfsize,z-halfsize);
				pt[1]=new point3D(x,y+halfsize,z-halfsize);
				pt[2]=new point3D(x,y+halfsize,z+halfsize);
				pt[3]=new point3D(x,y-halfsize,z+halfsize);
			}
		}
		color=color1;
	}
	
	public square3D(point3D[] pt,int size,int orientation,Color color){
		this.pt=pt;
		this.size=size;
		this.orientation=orientation;
		this.color=color;
	}

	public void moveto(int x,int y,int z){
		int halfsize=size/2;
		if(orientation==0){
			pt[0]=new point3D(x-halfsize,y-halfsize,z);
			pt[1]=new point3D(x+halfsize,y-halfsize,z);
			pt[2]=new point3D(x+halfsize,y+halfsize,z);
			pt[3]=new point3D(x-halfsize,y+halfsize,z);
		}else{
			if(orientation==1){
				pt[0]=new point3D(x-halfsize,y,z-halfsize);
				pt[1]=new point3D(x+halfsize,y,z-halfsize);
				pt[2]=new point3D(x+halfsize,y,z+halfsize);
				pt[3]=new point3D(x-halfsize,y,z+halfsize);
			}else{
				pt[0]=new point3D(x,y-halfsize,z-halfsize);
				pt[1]=new point3D(x,y+halfsize,z-halfsize);
				pt[2]=new point3D(x,y+halfsize,z+halfsize);
				pt[3]=new point3D(x,y-halfsize,z+halfsize);
			}
		}
	}

	public void translate(int transx,int transy,int transz){
		for(int i=0;i<4;i++){
			pt[i].translate(transx,transy,transz);
		}
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].rotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].rotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].rotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].addrotatecossin(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
		}
	}
	
	public void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1){
		for(int i=0;i<pt.length;i++){
			pt[i].set_rot_about_vector(cosval1,sinval1,ux1,uy1,uz1,centerx1,centery1,centerz1);
		}
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].addrotatedeg(degx1,degy1,degz1,centerx,centery,centerz);
		}
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].addrotaterad(radx1,rady1,radz1,centerx,centery,centerz);
		}
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
			pt[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		for(int i=0;i<4;i++){
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
		for(int i=0;i<4;i++){
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
		g.setColor(tempcolor);
	}
	
	public square3D clone(){
		point3D[] temp=new point3D[pt.length];
		for(int i=0;i<pt.length;i++) temp[i]=pt[i].clone();
		return new square3D(temp,size,orientation,color);
	}

}
