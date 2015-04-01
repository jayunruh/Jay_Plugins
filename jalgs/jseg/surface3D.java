/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.interpolation;

import java.awt.Polygon;
import java.awt.Rectangle;

public class surface3D{
	// this is the equivalent of a polygon in 3D
	// for simplicity, it exists as an ordered array from bottom to top of java
	// polygons with z indices
	// the upper and lower polygons are considered hard boundaries for the 3D
	// surface
	public Polygon[] polys;
	public int[] zpositions;
	public Rectangle xybounds;
	public int[] zbounds;
	public float[] centroid;

	public surface3D(Polygon[] polys,int bottom,int top){
		this.polys=polys;
		zpositions=new int[polys.length];
		float increment=(top-bottom)/(float)(polys.length-1);
		xybounds=polys[0].getBounds();
		for(int i=0;i<polys.length;i++){
			zpositions[i]=Math.round(bottom+increment*i);
			Rectangle temp=polys[i].getBounds();
			if(temp.x<xybounds.x)
				xybounds.x=temp.x;
			if(temp.y<xybounds.y)
				xybounds.y=temp.y;
			if((temp.x+temp.width)>(xybounds.x+xybounds.width))
				xybounds.width=(temp.x+temp.width)-xybounds.x;
			if((temp.y+temp.height)>(xybounds.y+xybounds.height))
				xybounds.height=(temp.y+temp.height)-xybounds.y;
		}
		zbounds=new int[]{bottom,top};
		centroid=null;
	}

	public surface3D(Polygon[] polys,int startz){
		this.polys=polys;
		zpositions=new int[polys.length];
		xybounds=polys[0].getBounds();
		for(int i=0;i<polys.length;i++){
			zpositions[i]=startz+i;
			Rectangle temp=polys[i].getBounds();
			if(temp.x<xybounds.x)
				xybounds.x=temp.x;
			if(temp.y<xybounds.y)
				xybounds.y=temp.y;
			if((temp.x+temp.width)>(xybounds.x+xybounds.width))
				xybounds.width=(temp.x+temp.width)-xybounds.x;
			if((temp.y+temp.height)>(xybounds.y+xybounds.height))
				xybounds.height=(temp.y+temp.height)-xybounds.y;
		}
		zbounds=new int[]{startz,startz+polys.length};
		centroid=null;
	}

	public float[] getcentroid(){
		double xsum=0.0;
		double ysum=0.0;
		double zsum=0.0;
		double count=0.0;
		for(int i=0;i<polys.length;i++){
			float zval=zpositions[i];
			Rectangle temp=polys[i].getBounds();
			for(int j=temp.y;j<(temp.y+temp.height);j++){
				for(int k=temp.x;k<(temp.x+temp.width);k++){
					if(polys[i].contains(k,j)){
						count+=1.0;
						xsum+=k;
						ysum+=j;
						zsum+=zval;
					}
				}
			}
		}
		centroid=new float[3];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		centroid[2]=(float)(zsum/count+0.5);
		return centroid;
	}

	public boolean contains(int[] point){
		return contains(point[0],point[1],point[2]);
	}

	public boolean contains(int x,int y,int z){
		if(z>zbounds[0]&&z<zbounds[1]){
			return polys[z-zbounds[0]].contains(x,y);
		}
		return false;
	}

	public Rectangle getXYBounds(){
		return xybounds;
	}

	public int[] getZBounds(){
		return zbounds;
	}

	public void contract(float pixels){
		if(centroid==null)
			getcentroid();
		float[][][] coords=new float[polys.length][][];
		float[] zcoords=new float[polys.length];
		for(int i=0;i<polys.length;i++){
			float zdist=(zpositions[i]-centroid[2]);
			float distance=Math.abs(zdist);
			float multiplier=(distance-pixels)/distance;
			zcoords[i]=zdist*multiplier+centroid[2];
		}
		float zmin=zcoords[0];
		float zmax=zcoords[polys.length-1];
		for(int i=0;i<polys.length;i++){
			float zdist=zpositions[i]-centroid[2];
			int[] xcoords=polys[i].xpoints;
			int[] ycoords=polys[i].ypoints;
			coords[i]=new float[2][xcoords.length];
			for(int j=0;j<xcoords.length;j++){
				float xdist=xcoords[j]-centroid[0];
				float ydist=ycoords[j]-centroid[1];
				float distance=(float)Math.sqrt(zdist*zdist+ydist*ydist+xdist*xdist);
				float multiplier=(distance-pixels)/distance;
				coords[i][0][j]=xdist*multiplier+centroid[0];
				coords[i][1][j]=ydist*multiplier+centroid[1];
			}
		}
		// now transform the new coordinates into the integer coordinate system
		// find the first z position that has data
		int zstart=(int)Math.ceil(zmin);
		int zend=(int)Math.floor(zmax);
		zbounds[0]=zstart;
		zbounds[1]=zend;
		zpositions=new int[zend-zstart+1];
		for(int i=0;i<=(zend-zstart);i++){
			zpositions[i]=zstart+i;
		}
		polys=new Polygon[zend-zstart+1];
		for(int i=0;i<=(zend-zstart);i++){
			float index=interpolation.get_float_index(zcoords,zpositions[i]);
			int prev=(int)index;
			float rem=index-prev;
			float[][] temp=coords[prev];
			if(rem>0.0f)
				temp=interpolate_polygon(coords[prev],coords[prev+1],rem);
			polys[i]=fp2poly(temp);
		}
		getcentroid();
	}

	public void expand(float pixels){
		contract(-pixels);
	}

	public float[][] interpolate_polygon(float[][] poly1,float[][] poly2,float fraction){
		int length1=poly1[0].length;
		int length2=poly2[0].length;
		float[][] tp1=poly1;
		float[][] tp2=poly2;
		float f=fraction;
		if(length1>length2){
			tp1=poly2;
			tp2=poly1;
			f=1.0f-fraction;
		}
		float[][] newpoly=new float[2][tp1[0].length];
		for(int i=0;i<tp1[0].length;i++){
			// find the closest point in the longer polygon
			float[] temp=interpolation.get_closest_point(tp2,new float[]{tp1[0][i],tp1[1][i]});
			// now interpolate between our current point and the closest point
			newpoly[0][i]=tp1[0][i]+f*(temp[0]-tp1[0][i]);
			newpoly[1][i]=tp1[1][i]+f*(temp[1]-tp1[1][i]);
		}
		return newpoly;
	}

	public float calcdist(float x1,float y1,float x2,float y2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}

	public Polygon fp2poly(float[][] poly){
		float[] bounds=getfpbounds(poly);
		int width=(int)(bounds[2]-bounds[0]);
		int height=(int)(bounds[3]-bounds[1]);
		float[] temp=new float[(width+4)*(height+4)];
		int xoff=(int)bounds[0]-2;
		int yoff=(int)bounds[2]-2;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(contains(j+xoff,i+yoff,poly)){
					temp[j+i*width]=1.0f;
				}
			}
		}
		findblobs3 fb=new findblobs3(width+4,height+4);
		Polygon poly2=fb.get_object_outline(temp,1);
		poly2.translate(-xoff,-yoff);
		return poly2;
	}

	/**
	 * Returns 'true' if the point (x,y) is inside this polygon. This is a Java
	 * version of the remarkably small C program by W. Randolph Franklin at
	 * http:
	 * //www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html#The
	 * %20C%20Code taken from ImageJ's FloatPolygon
	 */
	public boolean contains(float x,float y,float[][] poly){
		return contains(x,y,poly[0],poly[1]);
	}

	public boolean contains(float x,float y,float[] xpoints,float[] ypoints){
		boolean inside=false;
		int npoints=xpoints.length;
		for(int i=0,j=npoints-1;i<npoints;j=i++){
			if(((ypoints[i]>y)!=(ypoints[j]>y))&&(x<(xpoints[j]-xpoints[i])*(y-ypoints[i])/(ypoints[j]-ypoints[i])+xpoints[i]))
				inside=!inside;
		}
		return inside;
	}

	public float[] getfpbounds(float[][] poly){
		float xmin=poly[0][0];
		float xmax=poly[0][0];
		float ymin=poly[1][0];
		float ymax=poly[1][0];
		for(int i=1;i<poly[0].length;i++){
			if(poly[0][i]<xmin)
				xmin=poly[0][i];
			if(poly[0][i]>xmax)
				xmax=poly[0][i];
			if(poly[1][i]<ymin)
				ymin=poly[1][i];
			if(poly[1][i]>ymax)
				ymax=poly[1][i];
		}
		return new float[]{xmin,ymin,xmax,ymax};
	}

	public void draw(float[][] image,int width,int height,int xoff,int yoff,int zoff){
		byte[][] temp=new byte[image.length][image[0].length];
		for(int i=0;i<image.length;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(contains(j-xoff,k-yoff,zoff)){
						temp[i][k+j*width]=(byte)255;
					}
				}
			}
		}
	}

}
