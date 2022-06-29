/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import java.awt.Graphics;

public class j3Dutils{

	public static void draw3DLine(Graphics[] g,int x1,int y1,int z1,int x2,int y2,int z2,boolean thick){
		float length=(float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
		float dx=(x2-x1)/length;
		float dy=(y2-y1)/length;
		float dz=(z2-z1)/length;
		float xpos=x1, ypos=y1, zpos=z1;
		for(int i=0;i<(int)length;i++){
			int zind=(int)Math.round(zpos);
			if(zind>0 && zind<(g.length)){
				if(!thick) drawPoint(g[zind],(int)Math.round(xpos),(int)Math.round(ypos));
				else drawDot(g[zind],(int)Math.round(xpos),(int)Math.round(ypos));
			}
			xpos+=dx;
			ypos+=dy;
			zpos+=dz;
		}
		int zind=(int)Math.round(z2);
		if(zind>0 && zind<(g.length)){
			if(!thick) drawPoint(g[zind],(int)Math.round(x2),(int)Math.round(y2));
			else drawDot(g[zind],(int)Math.round(x2),(int)Math.round(y2));
		}
	}
	
	public static void drawLine(Graphics g,int x1,int y1,int x2,int y2,boolean thick){
		if(thick){
			drawThickLine(g,x1,y1,x2,y2);
		}else{
			g.drawLine(x1,y1,x2,y2);
		}
	}

	public static void drawThickLine(Graphics g,int x1,int y1,int x2,int y2){
		float length=(float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		float dx=(x2-x1)/length;
		float dy=(y2-y1)/length;
		float xpos=x1;
		float ypos=y1;
		for(int i=0;i<(int)length;i++){
			drawDot(g,(int)xpos,(int)ypos);
			xpos+=dx;
			ypos+=dy;
		}
		drawDot(g,x2,y2);
	}

	public static void drawDot(Graphics g,int x,int y){
		g.drawLine(x-1,y-1,x,y-1);
		g.drawLine(x-1,y,x,y);
	}
	
	public static void drawPoint(Graphics g,int x,int y){
		g.fillRect(x,y,1,1);
	}
	
	public static double[] norm_vec(double[] vec){
		double norm=Math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
		double[] temp={vec[0]/norm,vec[1]/norm,vec[2]/norm};
		return temp;
	}
	
	public static float get_3D_dist(point3D pt1,point3D pt2){
		return (float)Math.sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y)+(pt1.z-pt2.z)*(pt1.z-pt2.z));
	}
	
	public static float get_2D_dist(point3D pt1,point3D pt2){
		return (float)Math.sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y));
	}

}
