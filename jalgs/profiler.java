/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.awt.Polygon;

import jalgs.jseg.measure_object;

public class profiler{
	// here we have static methods to generate line and polyline thick profiles

	public static int getPolygonLength(Polygon polyroi,boolean connected){
		int[] xvals=polyroi.xpoints;
		int[] yvals=polyroi.ypoints;
		int nlines=polyroi.npoints-1;
		int length=0;
		for(int i=0;i<nlines;i++){
			int linelength=(int)Math.sqrt((xvals[i+1]-xvals[i])*(xvals[i+1]-xvals[i])+(yvals[i+1]-yvals[i])*(yvals[i+1]-yvals[i]))-1;
			if(linelength==0) linelength=1;
			length+=linelength;
		}
		if(connected){
			length+=(int)Math.sqrt((xvals[nlines]-xvals[0])*(xvals[nlines]-xvals[0])+(yvals[nlines]-yvals[0])*(yvals[nlines]-yvals[0]))-1;
		}else{
			length+=1;
		}
		return length;
	}
	
	public static int get3DPolygonLength(float[] xvals,float[] yvals,float[] zvals,boolean connected){
		int nlines=xvals.length-1;
		int length=0;
		for(int i=0;i<nlines;i++){
			int linelength=(int)get3DLength(new float[]{xvals[i],yvals[i],zvals[i],xvals[i+1],yvals[i+1],zvals[i+1]})-1;
			if(linelength==0) linelength=1;
			length+=linelength;
		}
		if(connected){
			length+=(int)get3DLength(new float[]{xvals[nlines],yvals[nlines],zvals[nlines],xvals[0],yvals[0],zvals[0]})-1;
		}else{
			length+=1;
		}
		return length;
	}

	/*************
	 * This outputs a straightened 2D profile for a polyline
	 * @param pixels
	 * @param width
	 * @param height
	 * @param polyroi
	 * @param connected
	 * @param linewidth
	 * @param or_index
	 * @return
	 */
	public static float[] getStraightened(Object pixels,int width,int height,Polygon polyroi,boolean connected,int linewidth,int or_index){
		if(linewidth<2){
			return getProfile(pixels,width,height,polyroi,connected,linewidth,or_index);
		}
		int length=getPolygonLength(polyroi,connected);
		int[] xvals=polyroi.xpoints;
		int[] yvals=polyroi.ypoints;
		int nlines=polyroi.npoints-1;
		float[] profile=new float[length*linewidth];
		int counter=0;
		for(int i=0;i<nlines;i++){
			float[] coords={xvals[i],yvals[i],xvals[i+1],yvals[i+1]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth+0.5f;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[(counter+k)*linewidth+j]=tempfloat[k];
				}
				if(i==(nlines-1)&&!connected){
					profile[(counter+templength-1)*linewidth+j]=tempfloat[templength-1];
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines],yvals[nlines],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[(counter+k)*linewidth+j]+=tempfloat[k];
				}
			}
		}
		return profile;
	}
	
	/*************
	 * This outputs a straightened 3D profile
	 * @param image
	 * @param width
	 * @param height
	 * @param xvals
	 * @param yvals
	 * @param zvals
	 * @param connected
	 * @param linewidth
	 * @param or_index
	 * @param zratio
	 * @return
	 */
	public static float[] get3DStraightened(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		if(linewidth<2){
			return get3DProfile(image,width,height,xvals,yvals,zvals,connected,linewidth,or_index,zratio);
		}
		int length=get3DPolygonLength(xvals,yvals,zvals,connected);
		int nlines=xvals.length-1;
		float[] profile=new float[length*linewidth];
		int counter=0;
		for(int i=0;i<nlines;i++){
			float[] coords={xvals[i],yvals[i],xvals[i+1],yvals[i+1]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth+0.5f;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[i],newcoords[2],newcoords[3],zvals[i+1]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				if(templength==1) templength=2;
				for(int k=0;k<templength-1;k++){
					profile[(counter+k)*linewidth+j]=tempfloat[k];
				}
				if(i==(nlines-1)&&!connected){
					if(tempfloat.length>1) profile[(counter+templength-1)*linewidth+j]=tempfloat[templength-1];
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines],yvals[nlines],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[nlines-1],newcoords[2],newcoords[3],zvals[0]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[(counter+k)*linewidth+j]+=tempfloat[k];
				}
			}
		}
		return profile;
	}

	/*************
	 * This outputs a 2D thick polyline profile
	 * @param pixels
	 * @param width
	 * @param height
	 * @param polyroi
	 * @param connected
	 * @param linewidth
	 * @param or_index
	 * @return
	 */
	public static float[] getProfile(Object pixels,int width,int height,Polygon polyroi,boolean connected,int linewidth,int or_index){
		int length=getPolygonLength(polyroi,connected);
		int[] xvals=polyroi.xpoints;
		int[] yvals=polyroi.ypoints;
		int nlines=polyroi.npoints-1;
		float[] profile=new float[length];
		int counter=0;
		for(int i=0;i<nlines;i++){
			float[] coords={xvals[i],yvals[i],xvals[i+1],yvals[i+1]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				if(templength==1) templength=2;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/linewidth;
				}
				//add on the final point if the profile is not connected
				if(i==(nlines-1)&&!connected && tempfloat.length>1){
					profile[counter+templength-1]+=tempfloat[templength-1]/linewidth;
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines],yvals[nlines],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/linewidth;
				}
			}
		}
		return profile;
	}

	public static float[] get3DProfile(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		int length=get3DPolygonLength(xvals,yvals,zvals,connected);
		int nlines=xvals.length-1;
		float[] profile=new float[length];
		int counter=0;
		for(int i=0;i<nlines;i++){
			float[] coords={xvals[i],yvals[i],xvals[i+1],yvals[i+1]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[i],newcoords[2],newcoords[3],zvals[i+1]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				if(templength==1) templength=2;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/linewidth;
				}
				if(i==(nlines-1)&&!connected && tempfloat.length>1){
					profile[counter+templength-1]+=tempfloat[templength-1]/linewidth;
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines],yvals[nlines],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-0.5f*(float)linewidth;
					distance*=(float)tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)tempoutsign*(float)j;
					}else{
						distance=-(float)tempoutsign*(float)j;
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[nlines-1],newcoords[2],newcoords[3],zvals[0]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/linewidth;
				}
			}
		}
		return profile;
	}
	
	public static float[] get3DThickProfile(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		float zoff=0.5f*linewidth;
		float[] tempzvals=new float[zvals.length];
		for(int i=0;i<zvals.length;i++) tempzvals[i]=zvals[i]-zoff;
		float[] profile=get3DProfile(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
		for(int i=1;i<linewidth;i++){
			for(int j=0;j<zvals.length;j++) tempzvals[j]+=1.0f;
			float[] tempprofile=get3DProfile(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
			for(int j=0;j<profile.length;j++) profile[j]+=tempprofile[j];
		}
		for(int i=0;i<profile.length;i++) profile[i]/=linewidth;
		return profile;
	}
	
	public static float[][] get3DThickStraightened(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		float zoff=0.5f*linewidth;
		float[] tempzvals=new float[zvals.length];
		for(int i=0;i<zvals.length;i++) tempzvals[i]=zvals[i]-zoff;
		float[][] profile=new float[linewidth][];
		profile[0]=get3DStraightened(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
		for(int i=1;i<linewidth;i++){
			for(int j=0;j<zvals.length;j++) tempzvals[j]+=1.0f;
			profile[i]=get3DStraightened(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
		}
		return profile;
	}
	
	/*************
	 * This gives a thick 2D line profile for a single segment
	 * @param coords: {x1,y1,x2,y2}
	 * @param image
	 * @param outsign: -1 for reverse orientation (perpendicular to line)
	 * @param or_index: 0, centered, 1, to the right, 2, to the left
	 * @param linewidth
	 * @param width: image width
	 * @param height: image height
	 * @return
	 */
	public static float[] get2DLineProfile(float[] coords,Object image,int outsign,int or_index,int linewidth,int width,int height){
		float[] profile=null;
		for(int j=0;j<linewidth;j++){
			float distance=0.0f;
			if(or_index==0){
				distance=(float)j-0.5f*(float)linewidth;
				distance*=(float)outsign;
			}else{
				if(or_index==1){
					distance=(float)outsign*(float)j;
				}else{
					distance=-(float)outsign*(float)j;
				}
			}
			float[] newcoords=getParallelLine(coords,distance);
			float[] tempfloat=getLineProfile(newcoords,image,width,height);
			if(profile==null) profile=new float[tempfloat.length];
			int profsize=tempfloat.length;
			if(tempfloat.length>profile.length) profsize=profile.length; //this should never happen
			for(int k=0;k<profsize;k++){
				profile[k]+=tempfloat[k]/linewidth;
			}
		}
		return profile;
	}

	/**********************
	 * This outputs a thin profile for a single segment
	 * @param coords: {x1,y1,x2,y2}
	 * @param image
	 * @param width
	 * @param height
	 * @return
	 */
	public static float[] getLineProfile(float[] coords,Object image,int width,int height){
		int length=(int)get2DLength(coords);
		float[] line=new float[length];
		float xinc,yinc;
		if(coords[2]!=coords[0]&&coords[3]!=coords[1]){
			float slope=(coords[3]-coords[1])/(coords[2]-coords[0]);
			xinc=(float)Math.sqrt(1.0/(1.0+slope*slope));
			if(coords[2]<coords[0]&&xinc>0.0f){
				xinc=-xinc;
			}
			yinc=slope*xinc;
		}else{
			if(coords[2]==coords[0]){
				xinc=0.0f;
				yinc=1.0f;
				if(coords[3]<coords[1]){
					yinc=-1.0f;
				}
			}else{
				xinc=1.0f;
				yinc=0.0f;
				if(coords[2]<coords[0]){
					xinc=-1.0f;
				}
			}
		}
		float x=coords[0];
		float y=coords[1];
		for(int i=0;i<length;i++){
			line[i]=interpolation.interp2D(image,width,height,x,y);
			x+=xinc;
			y+=yinc;
		}
		return line;
	}
	
	public static float[] get3DLineProfile(float[] coords,Object[] image,int width,int height,float zratio){
		//note that the coords are in real z units (xy units), not actual slice units
		float flength=get3DLength(coords);
		int length=(int)flength;
		float[] line=new float[length];
		float xinc=(coords[3]-coords[0])/length;
		float yinc=(coords[4]-coords[1])/length;
		float zinc=(coords[5]-coords[2])/length;
		zinc/=zratio;
		float x=coords[0];
		float y=coords[1];
		float z=coords[2]/zratio;
		for(int i=0;i<length;i++){
			line[i]=interpolation.interp3D(image,width,height,x,y,z);
			x+=xinc;
			y+=yinc;
			z+=zinc;
		}
		return line;
	}

	public static float get2DLength(float[] coords){
		return (float)Math.sqrt((coords[2]-coords[0])*(coords[2]-coords[0])+(coords[3]-coords[1])*(coords[3]-coords[1]));
	}
	
	public static float get3DLength(float[] coords){
		return (float)Math.sqrt((coords[3]-coords[0])*(coords[3]-coords[0])+(coords[4]-coords[1])*(coords[4]-coords[1])+(coords[5]-coords[2])*(coords[5]-coords[2]));
	}
	
	public static int outsign(Polygon polyroi,float[] coords){
		// are we inside or outside the roi?
		// will not work for complex roi shapes
		if(polyroi.npoints>2){
			float[] tempcoords={coords[0],coords[1],0.5f*(coords[0]+coords[2]),0.5f*(coords[1]+coords[3])};
			// want to find the outsign in the middle of the line
			float[] pluscoords=getParallelLine(tempcoords,1.0f);
			if(polyroi.contains(pluscoords[2],pluscoords[3])){
				return -1;
			}
		}
		return 1;
	}
	
	public static int outsign(float[] xvals,float[] yvals,float[] coords){
		// are we inside or outside the roi?
		// will not work for complex roi shapes
		if(xvals.length>2){
			float[] tempcoords={coords[0],coords[1],0.5f*(coords[0]+coords[2]),0.5f*(coords[1]+coords[3])};
			// want to find the outsign in the middle of the line
			float[] pluscoords=getParallelLine(tempcoords,1.0f);
			if(contains(xvals,yvals,pluscoords[2],pluscoords[3])){
				return -1;
			}
		}
		return 1;
	}
	
    /** Returns 'true' if the point (x,y) is inside this polygon. This is a Java
    version of the remarkably small C program by W. Randolph Franklin at
    http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html#The%20C%20Code
    */
    public static boolean contains(float[] xpoints,float[] ypoints,float x, float y) {
        boolean inside = false;
        int npoints=xpoints.length;
        for (int i=0, j=npoints-1; i<npoints; j=i++) {
            if (((ypoints[i]>y)!=(ypoints[j]>y)) &&
            (x<(xpoints[j]-xpoints[i])*(y-ypoints[i])/(ypoints[j]-ypoints[i])+xpoints[i]))
            inside = !inside;
        }
        return inside;
    }
    
    public static int outsign_from_prev(float[] xvals,float[] yvals,int index,int prevoutsign){
    	//here index is the first point of our line segment and the last point of the previous segment
    	//make a small polygon out of the two line segments and check whether the new positive position is on the same or different side as the previous
    	//if same return the previous outsign, if not return its negative
    	int next=index+1;
    	if(next>=xvals.length) next=0;
    	int prev=index-1;
    	float[] tempcoords2={xvals[index],yvals[index],0.5f*(xvals[index]+xvals[next]),0.5f*(yvals[index]+yvals[next])};
    	float[] pluscoords2=getParallelLine(tempcoords2,1.0f);
    	float[] tempcoords1={xvals[prev],yvals[prev],0.5f*(xvals[prev]+xvals[index]),0.5f*(yvals[prev]+yvals[index])};
    	float[] pluscoords1=getParallelLine(tempcoords1,1.0f);
    	float[] xvals2={xvals[prev],xvals[index],xvals[next]};
    	float[] yvals2={yvals[prev],yvals[index],yvals[next]};
    	boolean prevcontains=contains(xvals2,yvals2,pluscoords1[2],pluscoords1[3]);
    	boolean nextcontains=contains(xvals2,yvals2,pluscoords2[2],pluscoords2[3]);
    	if(prevcontains!=nextcontains) return -prevoutsign;
    	else return prevoutsign;
    }

	public static float[] getParallelLine(float[] coords,float distance){
		float xinc,yinc;
		if(coords[2]!=coords[0]&&coords[3]!=coords[1]){
			float slope=(coords[3]-coords[1])/(coords[2]-coords[0]);
			float newslope=-1.0f/slope;
			xinc=(float)Math.sqrt(1.0/(1.0+newslope*newslope));
			yinc=newslope*xinc;
		}else{
			if(coords[2]==coords[0]){
				xinc=1.0f;
				yinc=0.0f;
				if(coords[3]<coords[1]){
					xinc=-1.0f;
				}
			}else{
				xinc=0.0f;
				yinc=1.0f;
				if(coords[2]<coords[0]){
					yinc=-1.0f;
				}
			}
		}
		float[] oc=new float[4];
		oc[0]=distance*xinc+coords[0];
		oc[1]=distance*yinc+coords[1];
		oc[2]=distance*xinc+coords[2];
		oc[3]=distance*yinc+coords[3];
		return oc;
	}
	
	/*********************
	 * this plugin rotates a profile from its current vector to the z axis vector
	 * @param stack: the source image stack
	 * @param width: width of stack
	 * @param height: height of stack
	 * @param center: the point about which rotation occurs
	 * @param zratio: ratio of zres to xyres
	 * @param currvec: the current vector from the center point 
	 * @param newsize: the size (xy) of the rotated image
	 * @param zshift: the shift of the start of the stack "below" the vertex in xy units
	 * @param zsize: the size of the stack in xy units
	 * @return
	 */
	public static float[][] getRotated3DImage(Object[] stack,int width,int height,float[] center,float zratio,float[] currvec,int newsize,float zshift,int zsize){
		float[] zvec= {0.0f,0.0f,1.0f};
		//get the inner angle with the z axis
		float angle=measure_object.get_inner_angle(zvec,currvec);
		//the cross product will give us the axis to rotate about
		float[] crossprod=measure_object.crossProd(zvec,currvec);
		//normalize the cross product
		crossprod=measure_object.norm_vector(crossprod);
		float[][] rotmat=measure_object.getRotationMatrix(crossprod,angle);
		return getRotated3DImage(stack,width,height,center,zratio,rotmat,newsize,zshift,zsize);
	}
	
	/*************
	 * this plugin does the 3D realignment more robustly via a rotation matrix and interpolation
	 * @param stack: the source image stack
	 * @param width: width of stack
	 * @param height: height of stack
	 * @param center: the point about which rotation occurs
	 * @param zratio: ratio of zres to xyres
	 * @param rotmat: the 3 x 3 rotation matrix
	 * @param newsize: the size of the rotated image
	 * @param zshift: the shift of the start of the stack "below" the vertex in xy units
	 * @param zsize: the size of the stack in xy units
	 * @return
	 */
	public static float[][] getRotated3DImage(Object[] stack,int width,int height,float[] center,float zratio,float[][] rotmat,int newsize,float zshift,int zsize){
		//build a rotated image by rotating each voxel to its position in the original image and interpolating
		int rotsize=newsize;
		//int rotheight=4*rotsize;
		int rotheight=zsize;
		float[][] rotated=new float[rotheight][rotsize*rotsize];
		for(int i=0;i<rotheight;i++){
			float zpos=(float)(i-zshift);
			for(int j=0;j<rotsize;j++){
				float ypos=(float)(j-rotsize/2);
				for(int k=0;k<rotsize;k++){
					float xpos=(float)(k-rotsize/2);
					//multiply by the rotation matrix to transform in the old coordinates
					float xoff=rotmat[0][0]*xpos+rotmat[0][1]*ypos+rotmat[0][2]*zpos;
					float yoff=rotmat[1][0]*xpos+rotmat[1][1]*ypos+rotmat[1][2]*zpos;
					float zoff=rotmat[2][0]*xpos+rotmat[2][1]*ypos+rotmat[2][2]*zpos;
					//correct for anisotropic resolution
					zoff/=zratio;
					//add the vertex position back
					xoff+=center[0]; yoff+=center[1]; zoff+=(center[2]/zratio);
					//if(i==0 && j==0 && k==0) IJ.log(""+xoff+" , "+yoff+" , "+zoff);
					//finally interpolate the original image at these points
					rotated[i][k+j*rotsize]=interpolation.interp3D(stack,width,height,xoff,yoff,zoff);
				}
			}
		}
		return rotated;
	}

}
