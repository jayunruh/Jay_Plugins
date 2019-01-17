/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.interpolation;
import quickhull3d.Point3d;
import quickhull3d.QuickHull3D;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

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
	public static boolean contains(float x,float y,float[][] poly){
		return contains(x,y,poly[0],poly[1]);
	}

	public static boolean contains(float x,float y,float[] xpoints,float[] ypoints){
		boolean inside=false;
		int npoints=xpoints.length;
		for(int i=0,j=npoints-1;i<npoints;j=i++){
			if(((ypoints[i]>y)!=(ypoints[j]>y))&&(x<(xpoints[j]-xpoints[i])*(y-ypoints[i])/(ypoints[j]-ypoints[i])+xpoints[i]))
				inside=!inside;
		}
		return inside;
	}

	public static float[] getfpbounds(float[][] poly){
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
	
	/********************
	 * here we use QuickHull3D (https://github.com/Quickhull3d/quickhull3d) to do a 3D convex hull contour set
	 * @param points: the points to wrap the convex hull around
	 * @return
	 */
	public static float[][][] constructFromConvHull(float[][] points,float zres) {
		Point3d[] points2=new Point3d[points.length];
		//the limits of the convex hull are the same as those of the original point cloud
		float[] lims={points[0][0],points[0][0],points[0][1],points[0][1],points[0][2],points[0][2]};
		for(int i=0;i<points.length;i++) {
			points2[i]=new Point3d(points[i][0],points[i][1],points[i][2]);
			if(points[i][0]<lims[0]) lims[0]=points[i][0];
			if(points[i][0]>lims[1]) lims[1]=points[i][0];
			if(points[i][1]<lims[2]) lims[2]=points[i][1];
			if(points[i][1]>lims[3]) lims[3]=points[i][1];
			if(points[i][2]<lims[4]) lims[4]=points[i][2];
			if(points[i][2]>lims[5]) lims[5]=points[i][2];
		}
		QuickHull3D hull=new QuickHull3D();
		hull.build(points2); //builds the hull
		hull.triangulate(); //ensures all of the faces are triangles
		//each xy cross section of the convex hull should be a single (convex) polygon
		//the vertices lie along the edges of the faces
		int[][] faceindices=hull.getFaces();
		float[][][] faces2=new float[faceindices.length][3][];
		for(int i=0;i<faceindices.length;i++) {
			for(int j=0;j<3;j++) faces2[i][j]=points[faceindices[i][j]];
		}
		float[][][] edges=getEdges(faces2);
		//now scan through z and find all of the cross section vertices
		int nzs=1+(int)((lims[5]-lims[4])/zres);
		float[][][] crossverts=new float[nzs][][];
		for(int i=0;i<nzs;i++) {
			float zpos=lims[4]+i*zres;
			List<float[]> crossverts2=new ArrayList<float[]>();
			for(int j=0;j<edges.length;j++) {
				float[][] temp=get3DLineXYCrossing(edges[j],zpos);
				if(temp!=null) {
					//this line crosses the current plane, add it if it's unique
					addIfUnique(crossverts2,temp[0]);
					if(temp.length>1) addIfUnique(crossverts2,temp[1]);
				}
			}
			crossverts[i]=new float[crossverts2.size()][];
			for(int j=0;j<crossverts2.size();j++) {
				crossverts[i][j]=crossverts2.get(j);
			}
			//crossverts2.toArray(crossverts[i]);
		}
		//now scan through z again, running 2D convex hull to get the polygons ordered right
		for(int i=0;i<nzs;i++) {
			float zpos=lims[4]+i*zres;
			float[][] temp=getConvHull2D(crossverts[i]);
			if(temp!=null) {
				crossverts[i]=new float[temp[0].length][3];
				for(int j=0;j<temp[0].length;j++) {
					crossverts[i][j][0]=temp[0][j];
					crossverts[i][j][1]=temp[1][j];
					crossverts[i][j][2]=zpos;
				}
			} else {
				crossverts[i]=null;
			}
		}
		return crossverts;
	}
	
	/*************
	 * adds a point to a list if unique
	 * @param arr
	 * @param point
	 */
	public static void addIfUnique(List<float[]> arr,float[] point) {
		for(int i=0;i<arr.size();i++) {
			float[] temp=arr.get(i);
			if(isSamePoint(temp,point)) return;
		}
		arr.add(point);
	}
	
	/**************
	 * here we identify all of the shared lines (edges) in the faces and return only unique lines
	 * @param faces
	 * @return
	 */
	public static float[][][] getEdges(float[][][] faces) {
		//note that a line should belong to no more than two faces
		boolean[][] elim=new boolean[faces.length][2];
		for(int i=0;i<faces.length;i++){
			if(!elim[i][0]) {
				float[][] line1=new float[][]{faces[i][0],faces[i][1]};
				for(int j=i+1 ;j<faces.length;j++){
					float[][] line3={faces[j][0],faces[j][1]};
					if(!elim[j][0] && isSameLine(line1,line3)) {elim[j][0]=true; break;}
					float[][] line4={faces[j][1],faces[j][2]};
					if(!elim[j][1] && isSameLine(line1,line4)) {elim[j][1]=true; break;}
				}
			}
			if(!elim[i][1]) {
				float[][] line2=new float[][]{faces[i][1],faces[i][2]};
				for(int j=i+1 ;j<faces.length;j++){
					float[][] line3={faces[j][0],faces[j][1]};
					if(!elim[j][0] && isSameLine(line2,line3)) {elim[j][0]=true; break;}
					float[][] line4={faces[j][1],faces[j][2]};
					if(!elim[j][1] && isSameLine(line2,line4)) {elim[j][1]=true; break;}
				}
			}
		}
		int nedges=0;
		for(int i=0;i<faces.length;i++) {
			if(!elim[i][0]) nedges++;
			if(!elim[i][1]) nedges++;
		}
		float[][][] edges=new float[nedges][][];
		int counter=0;
		for(int i=0;i<faces.length;i++) {
			if(!elim[i][0]) {
				edges[counter]=new float[][]{faces[i][0],faces[i][1]};
				counter++;
			}
			if(!elim[i][1]) {
				edges[counter]=new float[][]{faces[i][1],faces[i][2]};
				counter++;
			}
		}
		return edges;
	}
	
	public static boolean isSameLine(float[][] line1, float[][] line2){
		  if(isSamePoint(line1[0],line2[0]) && isSamePoint(line1[1],line2[1])) return true;
		  if(isSamePoint(line1[1],line2[0]) && isSamePoint(line1[0],line2[1])) return true;
		  return false;
	}

	public static boolean isSamePoint(float[] pt1, float[] pt2){
	  if(pt1[0]==pt2[0] && pt1[1]==pt2[1] && pt1[2]==pt2[2]) return true;
	  return false;
	}
	
	public static float[][] get3DLineXYCrossing(float[][] line,float zpos) {
		if(line[0][2]<zpos && line[1][2]<zpos) return null; //line is completely below the plane
		if(line[0][2]>zpos && line[1][2]>zpos) return null; //line is completely above the plane
		if(line[0][2]==zpos) {
			if(line[1][2]==zpos) return line; //line is in this xy plane
			else return new float[][]{line[0]}; //first point is in the plane
		} else if(line[1][2]==zpos) {
			return new float[][]{line[1]}; //second point is in the plane
		}
		//float len=(float)Math.sqrt((line[1][0]-line[0][0])*(line[1][0]-line[0][0])+(line[1][1]-line[0][1])*(line[1][1]-line[0][1])+(line[1][2]-line[0][2])*(line[1][2]-line[0][2]));
		float zfrac=(zpos-line[0][2])/(line[1][2]-line[0][2]);
		float xinterp=line[0][0]+zfrac*(line[1][0]-line[0][0]);
		float yinterp=line[0][1]+zfrac*(line[1][1]-line[0][1]);
		return new float[][] {{xinterp,yinterp,zpos}};
	}
	
	/****************
	 * returns the convex hull face triangles
	 * @param points
	 * @return
	 */
	public static float[][][] getConvHull3D(float[][] points) {
		Point3d[] points2=new Point3d[points.length];
		for(int i=0;i<points.length;i++) {
			points2[i]=new Point3d(points[i][0],points[i][1],points[i][2]);
		}
		QuickHull3D hull=new QuickHull3D();
		hull.build(points2); //builds the hull
		hull.triangulate(); //ensures all of the faces are triangles
		int[][] faceindices=hull.getFaces();
		float[][][] faces2=new float[faceindices.length][3][];
		for(int i=0;i<faceindices.length;i++) {
			for(int j=0;j<3;j++) faces2[i][j]=points[faceindices[i][j]];
		}
		return faces2;
	}
	
	/***********
	 * a 2D giftwrap convex hull adapted from ImageJ for floating point data sets
	 * @param points
	 * @return
	 */
	public static float[][] getConvHull2D(float[][] points){
		if(points==null) return null;
		if(points.length<4) {
			float[][] temp=new float[2][points.length];
			for(int i=0;i<points.length;i++) {
				temp[0][i]=points[i][0];
				temp[1][i]=points[i][1];
			}
			return temp; //these are by definition convex hulls
		}
		int n = points.length;
		float[] xCoordinates=new float[n];
		float[] yCoordinates=new float[n];
		float[] lims= {points[0][0],points[0][0],points[0][1],points[0][1]};
		for(int i=0;i<n;i++) {
			xCoordinates[i]=points[i][0];
			yCoordinates[i]=points[i][1];
			if(xCoordinates[i]<lims[0]) lims[0]=xCoordinates[i];
			if(xCoordinates[i]>lims[1]) lims[1]=xCoordinates[i];
			if(yCoordinates[i]<lims[2]) lims[2]=yCoordinates[i];
			if(yCoordinates[i]>lims[3]) lims[3]=yCoordinates[i];
		}
        float[] xx = new float[n];
        float[] yy = new float[n];
        int n2 = 0;
        float smallestY = lims[2];
        float x, y;
        float smallestX = lims[1]; //initialize to greatest x
        //find the smallest x point at the smallest y value
        int p1 = 0;
        for (int i=0; i<n; i++) {
            x = xCoordinates[i];
            y = yCoordinates[i];
            if (y==smallestY && x<smallestX) {
                smallestX = x;
                p1 = i;
            }
        }
        int pstart = p1;
        int p2,p3;
        float x1, y1, x2, y2, x3, y3;
        float determinate;
        int count = 0;
        do {
            x1 = xCoordinates[p1];
            y1 = yCoordinates[p1];
            p2 = p1+1; if (p2==n) p2=0;
            x2 = xCoordinates[p2];
            y2 = yCoordinates[p2];
            p3 = p2+1; if (p3==n) p3=0;
            do {
                x3 = xCoordinates[p3];
                y3 = yCoordinates[p3];
                determinate = x1*(y2-y3)-y1*(x2-x3)+(y3*x2-y2*x3);
                if (determinate>0.0f)
                    {x2=x3; y2=y3; p2=p3;}
                p3 += 1;
                if (p3==n) p3 = 0;
            } while (p3!=p1);
            if (n2<n) { 
                xx[n2] = x1;
                yy[n2] = y1;
                n2++;
            } else {
                count++;
                if (count>10) return null;
            }
            p1 = p2;
        } while (p1!=pstart);
        return new float[][] {(float[])algutils.get_subarray(xx,0,n2),(float[])algutils.get_subarray(yy,0,n2)};
	}

}
