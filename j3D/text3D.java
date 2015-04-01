/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;

public class text3D extends element3D implements Cloneable{
	public point3D point;
	public String text; // shapes are square,plus,x,triangle
	public static final int fontsize=10;

	public text3D(String text1,int x,int y,int z,Color color1){
		point=new point3D(x,y,z);
		text=text1;
		color=color1;
	}
	
	public text3D(String text1,point3D pt,Color color1){
		point=pt;
		text=text1;
		color=color1;
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
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		point.transform_negative_perspective(horizon_dist,centerx,centery,centerz);
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
		Font temp=g.getFont().deriveFont((float)fontsize);
		g.setFont(temp);
		g.drawString(text,point.rx,point.ry);
		g.setColor(tempcolor);
	}
	
	public text3D clone(){
		return new text3D(text,point.clone(),color);
	}

}
