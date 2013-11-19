/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.awt.*;

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
			length+=(int)Math.sqrt((xvals[nlines-1]-xvals[0])*(xvals[nlines-1]-xvals[0])+(yvals[nlines-1]-yvals[0])*(yvals[nlines-1]-yvals[0]))-1;
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
			length+=(int)get3DLength(new float[]{xvals[nlines-1],yvals[nlines-1],zvals[nlines-1],xvals[0],yvals[0],zvals[0]})-1;
		}else{
			length+=1;
		}
		return length;
	}

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
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
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
			float[] coords={xvals[nlines-1],yvals[nlines-1],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
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
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
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
					profile[(counter+templength-1)*linewidth+j]=tempfloat[templength-1];
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines-1],yvals[nlines-1],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
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
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				if(templength==1) templength=2;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/(float)linewidth;
				}
				if(i==(nlines-1)&&!connected && tempfloat.length>1){
					profile[counter+templength-1]+=tempfloat[templength-1]/(float)linewidth;
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines-1],yvals[nlines-1],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(polyroi,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] tempfloat=getLineProfile(newcoords,pixels,width,height);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/(float)linewidth;
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
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[i],newcoords[2],newcoords[3],zvals[i+1]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				if(templength==1) templength=2;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/(float)linewidth;
				}
				if(i==(nlines-1)&&!connected && tempfloat.length>1){
					profile[counter+templength-1]+=tempfloat[templength-1]/(float)linewidth;
				}
			}
			counter+=(templength-1);
		}
		if(connected){
			float[] coords={xvals[nlines-1],yvals[nlines-1],xvals[0],yvals[0]};
			int templength=0;
			int tempoutsign=outsign(xvals,yvals,coords);
			for(int j=0;j<linewidth;j++){
				float distance=0.0f;
				if(or_index==0){
					distance=(float)j-((float)linewidth)/2.0f;
					distance*=tempoutsign;
				}else{
					if(or_index==1){
						distance=(float)(tempoutsign*j);
					}else{
						distance=(float)(-tempoutsign*j);
					}
				}
				float[] newcoords=getParallelLine(coords,distance);
				float[] newcoords3D={newcoords[0],newcoords[1],zvals[nlines-1],newcoords[2],newcoords[3],zvals[0]};
				float[] tempfloat=get3DLineProfile(newcoords3D,image,width,height,zratio);
				templength=tempfloat.length;
				for(int k=0;k<templength-1;k++){
					profile[counter+k]+=tempfloat[k]/(float)linewidth;
				}
			}
		}
		return profile;
	}
	
	public static float[] get3DThickProfile(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		float zoff=0.5f*(float)linewidth;
		float[] tempzvals=new float[zvals.length];
		for(int i=0;i<zvals.length;i++) tempzvals[i]=zvals[i]-zoff;
		float[] profile=get3DProfile(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
		for(int i=1;i<linewidth;i++){
			for(int j=0;j<zvals.length;j++) tempzvals[j]+=1.0f;
			float[] tempprofile=get3DProfile(image,width,height,xvals,yvals,tempzvals,connected,linewidth,or_index,zratio);
			for(int j=0;j<profile.length;j++) profile[j]+=tempprofile[j];
		}
		for(int i=0;i<profile.length;i++) profile[i]/=(float)linewidth;
		return profile;
	}
	
	public static float[][] get3DThickStraightened(Object[] image,int width,int height,float[] xvals,float[] yvals,float[] zvals,boolean connected,int linewidth,int or_index,float zratio){
		float zoff=0.5f*(float)linewidth;
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

	public static float[] getLineProfile(float[] coords,Object image,int width,int height){
		int length=(int)get2DLength(coords);
		float[] line=new float[length];
		float xinc,yinc;
		if(coords[2]!=coords[0]&&coords[3]!=coords[1]){
			float slope=(coords[3]-coords[1])/(coords[2]-coords[0]);
			xinc=(float)Math.sqrt(1.0/(1.0+(double)(slope*slope)));
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
		float xinc=(coords[3]-coords[0])/(float)length;
		float yinc=(coords[4]-coords[1])/(float)length;
		float zinc=(coords[5]-coords[2])/(float)length;
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
			if(polyroi.contains((double)pluscoords[2],(double)pluscoords[3])){
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
			xinc=(float)Math.sqrt(1.0/(1.0+(double)(newslope*newslope)));
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

}
