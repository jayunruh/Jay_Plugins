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
import jalgs.jstatistics;
import jalgs.profiler;

import java.awt.Polygon;
import java.awt.Rectangle;

public class measure_object{

	static final double HALFPI=1.5707963267949;

	public static float[][] centroids(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[][] outvals=new float[nobjects][];
		for(int i=1;i<=nobjects;i++){
			outvals[i-1]=measure_object.centroid(image,i,width,height);
		}
		return outvals;
	}

	public static float[] areas(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[] outvals=new float[nobjects];
		for(int i=1;i<=nobjects;i++){
			outvals[i-1]=measure_object.area(image,i,width,height);
		}
		return outvals;
	}

	public static float[] widths(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[] outvals=new float[nobjects];
		for(int i=1;i<=nobjects;i++){
			outvals[i-1]=measure_object.objwidth(image,i,width,height);
		}
		return outvals;
	}

	public static float[][] extents(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		int[] wstarts=new int[nobjects];
		int[] wends=new int[nobjects];
		for(int i=0;i<nobjects;i++)
			wstarts[i]=width-1;
		int[] hstarts=new int[nobjects];
		int[] hends=new int[nobjects];
		for(int i=0;i<nobjects;i++)
			hstarts[i]=height-1;
		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
				if(image[i+j*width]>0.0f){
					int id=(int)image[i+j*width]-1;
					if(i<wstarts[id])
						wstarts[id]=i;
					if(i>wends[id])
						wends[id]=i;
					if(j<hstarts[id])
						hstarts[id]=j;
					if(j>hends[id])
						hends[id]=j;
				}
			}
		}
		float[][] outvals=new float[nobjects][2];
		for(int i=0;i<nobjects;i++){
			outvals[i][0]=wends[i]-wstarts[i];
			outvals[i][1]=hends[i]-hstarts[i];
		}
		return outvals;
	}

	public static float[] areas2(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[] outvals=new float[nobjects];
		for(int i=0;i<width*height;i++){
			if(image[i]>0.0f){
				int id=(int)image[i]-1;
				outvals[id]+=1.0f;
			}
		}
		return outvals;
	}

	public static float[] sums(float[] image,Object measurement,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[] outvals=new float[nobjects];
		boolean isshort=(measurement instanceof short[]);
		boolean isbyte=(measurement instanceof byte[]);
		for(int i=0;i<width*height;i++){
			if(image[i]>0.0f){
				int id=(int)image[i]-1;
				if(isbyte)
					outvals[id]+=((byte[])measurement)[i]&0xff;
				else if(isshort)
					outvals[id]+=((short[])measurement)[i]&0xffff;
				else
					outvals[id]+=((float[])measurement)[i];
			}
		}
		return outvals;
	}

	public static float[][] sums(float[] image,Object[] measurement,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[][] outvals=new float[nobjects][measurement.length];
		boolean isshort=(measurement[0] instanceof short[]);
		boolean isbyte=(measurement[0] instanceof byte[]);
		for(int i=0;i<width*height;i++){
			if(image[i]>0.0f){
				int id=(int)image[i]-1;
				for(int j=0;j<measurement.length;j++){
					if(isbyte)
						outvals[id][j]+=((byte[])measurement[j])[i]&0xff;
					else if(isshort)
						outvals[id][j]+=((short[])measurement[j])[i]&0xffff;
					else
						outvals[id][j]+=((float[])measurement[j])[i];
				}
			}
		}
		return outvals;
	}

	public static float[] heights(float[] image,int width,int height){
		int nobjects=(int)measure_object.maxarray(image);
		float[] outvals=new float[nobjects];
		for(int i=1;i<=nobjects;i++){
			outvals[i-1]=measure_object.objheight(image,i,width,height);
		}
		return outvals;
	}

	public static float maxarray(float[] input){
		float max=input[0];
		for(int i=1;i<input.length;i++){
			if(input[i]>max){
				max=input[i];
			}
		}
		return max;
	}

	public static float area(float[] image,float id,int width,int height){
		float count=0.0f;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]==id){
					count+=1.0f;
				}
			}
		}
		return count;
	}
	
	/*************
	 * this method takes a yeast cell mask (object) and uses the transmitted light to find the bud neck and then returns only the bud area portion
	 * @param ellipseparams
	 * @param objects
	 * @param transimg
	 * @param width
	 * @param height
	 * @return
	 */
	public static Object[] getBudObject(float[] ellipseparams,float[] objects,float[] transimg,int width,int height){
		//start by creating a profile along the major axis as identified by the ellipse
		float angle=ellipseparams[2];
		float major=ellipseparams[3];
		if(major<5.0f) return null;
		float x=ellipseparams[0];
		float y=ellipseparams[1];
		float[] vec1=measure_object.get_unit_vector(angle);
		//float[] perpvec={vec1[1],-vec1[0]};
		float[] coords={x-vec1[0]*0.5f*major,y-vec1[1]*0.5f*major,x+vec1[0]*0.5f*major,y+vec1[1]*0.5f*major};
		float[] profile=profiler.get2DLineProfile(coords,transimg,1,0,4,width,height);
		//invert the profile so we can search for maxima
		for(int i=0;i<profile.length;i++) profile[i]=-profile[i];
		//search the first and last 25% for the edges and the central 50% for the bud neck
		int temp=(int)(0.25f*(float)profile.length);
		int temp2=profile.length-1-temp;
		float[] first=(float[])algutils.get_subarray(profile,0,temp);
		float[] last=(float[])algutils.get_subarray(profile,temp2,temp);
		float[] mid=(float[])algutils.get_subarray(profile,temp,temp+temp);
		float firstpos=jstatistics.getstatistic("maxpos",first,null);
		float lastpos=jstatistics.getstatistic("maxpos",last,null);
		lastpos+=(float)temp2;
		float midpos=jstatistics.getstatistic("maxpos",mid,null);
		midpos+=(float)temp;
		//divide the object perpendicular to the bud neck: save only the bud overlay
		boolean lastbud=false;
		if((lastpos-midpos)<(midpos-firstpos)) lastbud=true;
		float[] midcoords={midpos*vec1[0]+coords[0],midpos*vec1[1]+coords[1]};
		float[] newobject=new float[width*height];
		float id=objects[(int)x+width*(int)y];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(objects[j+i*width]==id){
					float[] vec2={(float)j-midcoords[0],(float)i-midcoords[1]};
					float dotprod=vec1[0]*vec2[0]+vec1[1]*vec2[1];
					newobject[j+i*width]=id;
					if(lastbud){if(dotprod<0.0f) newobject[j+i*width]=0.0f;}
					else if(dotprod>0.0f) newobject[j+i*width]=0.0f;
				}
			}
		}
		float[] measurements={firstpos,midpos,lastpos,midcoords[0],midcoords[1],id};
		return new Object[]{newobject,measurements};
	}
	
	/**************
	 * this calculates the tranformation matrix for rotation by angle theta about the axis given by axisvec--must be unit length
	 * @param axisvec: unit length axis of rotation vector
	 * @param theta: angle to rotate
	 * @return
	 */
	public static float[][] getRotationMatrix(float[] axisvec,float theta){
		float[] u=axisvec.clone();
		float cost=(float)Math.cos(theta);
		float sint=(float)Math.sin(theta);
		float[][] r=new float[3][];
		r[0]=new float[]{cost+u[0]*u[0]*(1.0f-cost), u[0]*u[1]*(1.0f-cost)-u[2]*sint, u[0]*u[2]*(1.0f-cost)+u[1]*sint};
		r[1]=new float[]{u[0]*u[1]*(1.0f-cost)+u[2]*sint, cost+u[1]*u[1]*(1.0f-cost), u[1]*u[2]*(1.0f-cost)-u[0]*sint};
		r[2]=new float[]{u[0]*u[2]*(1.0f-cost)-u[1]*sint, u[1]*u[2]*(1.0f-cost)+u[2]*sint, u[2]*u[2]*(1.0f-cost)+cost};
		return r;
	}
	
	/**********************
	 * calculates the cross product of two vectors
	 * @param u: vector 1
	 * @param v: vector 2
	 * @return
	 */
	public static float[] crossProd(float[] u,float[] v){
		float[] cp=new float[3];
		cp[0]=(u[1]*v[2]-u[2]*v[1]);
		cp[1]=(u[2]*v[0]-u[0]*v[2]);
		cp[2]=(u[0]*v[1]-u[1]*v[0]);
		return cp;
	}
	
	/**********************
	 * calculates the cross product of two vectors
	 * @param u: vector 1
	 * @param v: vector 2
	 * @return
	 */
	public static double[] crossProd(double[] u,double[] v){
		double[] cp=new double[3];
		cp[0]=(u[1]*v[2]-u[2]*v[1]);
		cp[1]=(u[2]*v[0]-u[0]*v[2]);
		cp[2]=(u[0]*v[1]-u[1]*v[0]);
		return cp;
	}

	public static float objwidth(float[] image,float id,int width,int height){
		int minpos=width-1;
		int maxpos=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]==id){
					if(j<minpos)
						minpos=j;
					if(j>maxpos)
						maxpos=j;
				}
			}
		}
		return maxpos-minpos;
	}

	public static float objheight(float[] image,float id,int width,int height){
		int minpos=height-1;
		int maxpos=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]==id){
					if(i<minpos)
						minpos=i;
					if(i>maxpos)
						maxpos=i;
				}
			}
		}
		return maxpos-minpos;
	}

	public static float area(Polygon poly){
		double count=0.0;
		Rectangle r=poly.getBounds();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					count+=1.0;
				}
			}
		}
		return (float)count;
	}

	public static float perimeter(Polygon poly){
		int[] xpts=poly.xpoints;
		int[] ypts=poly.ypoints;
		int npts=poly.npoints;
		double dist=Math.sqrt((xpts[0]-xpts[npts-1])*(xpts[0]-xpts[npts-1])+(ypts[0]-ypts[npts-1])*(ypts[0]-ypts[npts-1]));
		for(int i=1;i<poly.npoints;i++){
			dist+=Math.sqrt((xpts[i]-xpts[i-1])*(xpts[i]-xpts[i-1])+(ypts[i]-ypts[i-1])*(ypts[i]-ypts[i-1]));
		}
		return (float)dist;
	}

	public static float circularity(Polygon poly){
		float perim2=perimeter(poly);
		perim2*=perim2;
		float a=area(poly);
		return 4.0f*(float)Math.PI*(a/perim2);
	}

	public static float[] centroid(float[] image,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				count+=image[j+i*width];
				xsum+=image[j+i*width]*(double)j;
				ysum+=image[j+i*width]*(double)i;
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] centroid(float[] image,int width,int height,Rectangle r){
		if(r==null)
			return centroid(image,width,height);
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				count+=image[j+i*width];
				xsum+=image[j+i*width]*(double)j;
				ysum+=image[j+i*width]*(double)i;
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] centroid(float[] image,float id,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]==id){
					count+=1.0;
					xsum+=j;
					ysum+=i;
				}
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] centroid(Polygon poly){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		Rectangle r=poly.getBounds();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					count+=1.0;
					xsum+=j;
					ysum+=i;
				}
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] centroid(byte[] image,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]!=(byte)0){
					count+=1.0;
					xsum+=j;
					ysum+=i;
				}
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] centroid(byte[] image,int width,int height,Rectangle r){
		if(r==null)
			return centroid(image,width,height);
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(image[j+i*width]!=(byte)0){
					count+=1.0;
					xsum+=j;
					ysum+=i;
				}
			}
		}
		float[] centroid=new float[2];
		centroid[0]=(float)(xsum/count+0.5);
		centroid[1]=(float)(ysum/count+0.5);
		return centroid;
	}

	public static float[] radial_profile(float[] image,float id,int width,int height,int nangles){
		float[] cent=centroid(image,id,width,height);
		return radial_profile(image,id,cent[0],cent[1],width,height,nangles);
	}

	public static float[] radial_profile(float[] image,float id,float xcenter,float ycenter,int width,int height,int nangles){
		int maxrad=(int)xcenter-2;
		if((width-3-(int)xcenter)<maxrad){
			maxrad=(width-3-(int)xcenter);
		}
		if(((int)ycenter-2)<maxrad){
			maxrad=(int)ycenter-2;
		}
		if((height-3-(int)ycenter)<maxrad){
			maxrad=(height-3-(int)ycenter);
		}
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[] distances=new float[nangles];
		float[] sel_image=select_object(image,id);
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolate(sel_image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				if(temp>0.5f){
					distances[i]=j;
				}
			}
		}
		return distances;
	}

	public static float[] radial_profile(byte[] image,int width,int height,int nangles){
		float[] cent=centroid(image,width,height);
		return radial_profile(image,cent[0],cent[1],width,height,nangles);
	}

	public static float[] radial_profile(byte[] image,float xcenter,float ycenter,int width,int height,int nangles){
		int maxrad=(int)xcenter-2;
		if((width-3-(int)xcenter)<maxrad){
			maxrad=(width-3-(int)xcenter);
		}
		if(((int)ycenter-2)<maxrad){
			maxrad=(int)ycenter-2;
		}
		if((height-3-(int)ycenter)<maxrad){
			maxrad=(height-3-(int)ycenter);
		}
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[] distances=new float[nangles];
		float[] sel_image=select_object(image);
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolate(sel_image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				if(temp>0.5f){
					distances[i]=j;
				}
			}
		}
		return distances;
	}

	public static float[] radial_profile(byte[] image,float xcenter,float ycenter,int width,int height,float startangle,float dtheta,int maxrad,int nangles){
		float[] distances=new float[nangles];
		float[] sel_image=select_object(image);
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta+startangle;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolate(sel_image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				if(temp>0.5f){
					distances[i]=j;
				}
			}
		}
		return distances;
	}

	public static float[] radial_thickness_profile(byte[] image,int width,int height,int nangles){
		float[] cent=centroid(image,width,height);
		return radial_thickness_profile(image,cent[0],cent[1],width,height,nangles);
	}

	public static float[] radial_thickness_profile(byte[] image,float xcenter,float ycenter,int width,int height,int nangles){
		int maxrad=(int)xcenter-2;
		if((width-3-(int)xcenter)>maxrad){
			maxrad=(width-3-(int)xcenter);
		}
		if(((int)ycenter-2)>maxrad){
			maxrad=(int)ycenter-2;
		}
		if((height-3-(int)ycenter)>maxrad){
			maxrad=(height-3-(int)ycenter);
		}
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[] thickness=new float[nangles];
		float[] sel_image=select_object(image);
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolate(sel_image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				if(temp>0.5f){
					thickness[i]+=1.0;
				}
			}
		}
		return thickness;
	}

	public static float[] radial_sum_profile(float[] image,int width,int height,int nangles){
		float[] cent=centroid(image,width,height);
		return radial_sum_profile(image,cent[0],cent[1],width,height,nangles);
	}

	public static float[] radial_sum_profile(float[] image,float xcenter,float ycenter,int width,int height,int nangles){
		int maxrad=(int)xcenter-2;
		if((width-3-(int)xcenter)>maxrad){
			maxrad=(width-3-(int)xcenter);
		}
		if(((int)ycenter-2)>maxrad){
			maxrad=(int)ycenter-2;
		}
		if((height-3-(int)ycenter)>maxrad){
			maxrad=(height-3-(int)ycenter);
		}
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[] sum=new float[nangles];
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolation.interp2D(image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				sum[i]+=temp;
			}
		}
		return sum;
	}

	public static float[] radial_max_profile(float[] image,int width,int height,int nangles){
		float[] cent=centroid(image,width,height);
		return radial_max_profile(image,cent[0],cent[1],width,height,nangles);
	}

	public static float[] radial_max_profile(float[] image,float xcenter,float ycenter,int width,int height,int nangles){
		int maxrad=(int)xcenter-2;
		if((width-3-(int)xcenter)>maxrad){
			maxrad=(width-3-(int)xcenter);
		}
		if(((int)ycenter-2)>maxrad){
			maxrad=(int)ycenter-2;
		}
		if((height-3-(int)ycenter)>maxrad){
			maxrad=(height-3-(int)ycenter);
		}
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[] max=new float[nangles];
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				float temp=interpolation.interp2D(image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
				if(temp>max[i]){
					max[i]=temp;
				}
			}
		}
		return max;
	}

	public static float[][] radial_straight_profile(float[] image,float xcenter,float ycenter,int width,int height,int nangles,int maxrad){
		float dtheta=2.0f*(float)Math.PI/nangles;
		float[][] profiles=new float[nangles][maxrad];
		for(int i=0;i<nangles;i++){
			float theta=i*dtheta;
			float[] vec=get_unit_vector(theta);
			for(int j=0;j<maxrad;j++){
				profiles[i][j]=interpolation.interp2D(image,width,height,xcenter+vec[0]*j,ycenter-vec[1]*j);
			}
		}
		return profiles;
	}

	public static float[] select_object(float[] image,float id){
		float[] sel_image=new float[image.length];
		for(int i=0;i<image.length;i++){
			if(image[i]==id){
				sel_image[i]=1.0f;
			}
		}
		return sel_image;
	}

	public static float[] select_object(byte[] image){
		float[] sel_image=new float[image.length];
		for(int i=0;i<image.length;i++){
			if(image[i]!=(byte)0){
				sel_image[i]=1.0f;
			}
		}
		return sel_image;
	}

	public static float[] get_unit_vector(float theta1){
		float theta=theta1;
		if(theta>2.0f*(float)Math.PI){
			theta=theta1%(2.0f*(float)Math.PI);
		}
		float[] vector=new float[2];
		if(theta==0.5f*(float)Math.PI){
			vector[1]=1.0f;
			return vector;
		}
		if(theta==1.5f*(float)Math.PI){
			vector[1]=-1.0f;
			return vector;
		}
		if(theta==0.0f){
			vector[0]=1.0f;
			return vector;
		}
		if(theta==(float)Math.PI){
			vector[0]=-1.0f;
			return vector;
		}
		float slope=(float)Math.tan(theta);
		float xinc=(float)Math.abs(1.0/Math.sqrt(1.0+slope*slope));
		float yinc=(float)Math.abs(slope/Math.sqrt(1.0+slope*slope));
		if(theta>(float)Math.PI){
			yinc=-yinc;
		}
		if(theta>(0.5f*(float)Math.PI)&&theta<(1.5f*(float)Math.PI)){
			xinc=-xinc;
		}
		vector[0]=xinc;
		vector[1]=yinc;
		return vector;
	}

	public static float[] get_unit_vector3D(float theta1,float phi1){
		float phi=phi1;
		if(phi<=-0.5f*(float)Math.PI){
			return new float[]{0.0f,0.0f,-1.0f};
		}
		if(phi>=0.5f*(float)Math.PI){
			return new float[]{0.0f,0.0f,1.0f};
		}
		float[] vec2d=get_unit_vector(theta1);
		if(phi==0.0f){
			return new float[]{vec2d[0],vec2d[1],0.0f};
		}
		float cosval=(float)Math.cos(phi);
		float[] vec3d={vec2d[0]*cosval,vec2d[1]*cosval,(float)Math.sin(phi)};
		return vec3d;
	}

	public static float[] rotate_vector3D(float[] input,float[] rotation){
		double tempx=input[0];
		double tempy=input[1];
		double tempz=input[2];
		if(rotation[0]!=0.0){ // rotation about the z axis
			double sinval=Math.sin(rotation[0]);
			double cosval=Math.cos(rotation[0]);
			double tempx1=tempx*cosval-tempy*sinval;
			double tempy1=tempx*sinval+tempy*cosval;
			tempx=tempx1;
			tempy=tempy1;
		}
		if(rotation[1]!=0.0){ // rotation about the y axis
			double sinval=Math.sin(rotation[1]);
			double cosval=Math.cos(rotation[1]);
			double tempx1=tempx*cosval+tempz*sinval;
			double tempz1=-tempx*sinval+tempz*cosval;
			tempx=tempx1;
			tempz=tempz1;
		}
		if(rotation[2]!=0.0){
			double sinval=Math.sin(rotation[2]);
			double cosval=Math.cos(rotation[2]);
			double tempy1=tempy*cosval-tempz*sinval;
			double tempz1=tempy*sinval+tempz*cosval;
			tempy=tempy1;
			tempz=tempz1;
		}
		return new float[]{(float)tempx,(float)tempy,(float)tempz};
	}

	public static float get_angle(float[] vec){
		// the returned angle starts at zero for right pointing vectors
		// and increases as the vector moves counter clockwise up to 2 pi
		double angle=Math.atan2(vec[1],vec[0]);
		if(angle>=0.0){
			return (float)angle;
		}else{
			return (float)(angle+2.0*Math.PI);
		}
	}

	public static float get_angle(float x,float y,float xorigin,float yorigin){
		float[] vec={x-xorigin,y-yorigin};
		return get_angle(vec);
	}

	public static float get_angle2(float[] vec){
		// the returned angle starts at zero for right pointing vectors
		// and increases as the vector moves counter clockwise up to pi
		// the angle goes negative as the vector moves clockwise down to -pi
		double angle=Math.atan2(vec[1],vec[0]);
		return (float)angle;
	}

	public static float get_angle2(float x,float y,float xorigin,float yorigin){
		float[] vec={x-xorigin,y-yorigin};
		return get_angle2(vec);
	}

	public static float get_inner_angle(float[] vec1,float[] vec2){
		float[] tempvec1=norm_vector(vec1);
		float[] tempvec2=norm_vector(vec2);
		float dotprod=0.0f;
		for(int i=0;i<vec1.length;i++) dotprod+=tempvec1[i]*tempvec2[i];
		return (float)Math.acos(dotprod);
	}

	public static float get_inner_angle_points(float[] xpts,float[] ypts){
		float[] vec1={xpts[0]-xpts[1],ypts[0]-ypts[1]};
		float[] vec2={xpts[2]-xpts[1],ypts[2]-ypts[1]};
		return get_inner_angle(vec1,vec2);
	}

	public static float get_inner_angle_points(int[] xpts,int[] ypts){
		float[] vec1={xpts[0]-xpts[1],ypts[0]-ypts[1]};
		float[] vec2={xpts[2]-xpts[1],ypts[2]-ypts[1]};
		return get_inner_angle(vec1,vec2);
	}

	public static float[] norm_vector(float[] vec){
		float length2=0.0f;
		for(int i=0;i<vec.length;i++) length2+=vec[i]*vec[i];
		float length=(float)Math.sqrt(length2);
		float[] retvec=new float[vec.length];
		for(int i=0;i<vec.length;i++) retvec[i]=vec[i]/length;
		return retvec;
	}

	public static float interpolate(float[] image,int width,int height,float x,float y){
		if(x<0.0f){
			return 0.0f;
		}
		if(y<0.0f){
			return 0.0f;
		}
		if(x>width-2){
			return 0.0f;
		}
		if(y>height-2){
			return 0.0f;
		}
		int xprev=(int)x;
		int yprev=(int)y;
		float xrem=x-xprev;
		float yrem=y-yprev;
		float int1=yrem*(image[xprev+1+width*yprev]-image[xprev+width*yprev])+image[xprev+width*yprev];
		float int2=yrem*(image[xprev+1+width*(yprev+1)]-image[xprev+width*(yprev+1)])+image[xprev+width*(yprev+1)];
		return int1+xrem*(int2-int1);
	}

	public static float[] get_ellipse_parameters(Polygon poly){
		// this is a adapted from Bod Rodieck's code in ImageJ (EllipseFitter
		// class) that finds the ellipse
		// parameters that best match the spatial moments of the Roi
		int npixels=0;
		double xsum=0;
		double ysum=0;
		double x2sum=0;
		double y2sum=0;
		double xysum=0;
		Rectangle r=poly.getBounds();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					xsum+=j;
					ysum+=i;
					x2sum+=j*j;
					y2sum+=i*i;
					xysum+=j*i;
					npixels++;
				}
			}
		}
		/*
		 * x2sum+=0.08333333*npixels; y2sum+=0.08333333*npixels;
		 */
		xsum/=npixels;
		ysum/=npixels;
		x2sum/=npixels;
		y2sum/=npixels;
		xysum/=npixels;
		double xvar=x2sum-xsum*xsum;
		double yvar=y2sum-ysum*ysum;
		double xyvar=xysum-xsum*ysum;
		double m4=4.0*Math.abs(yvar*xvar-xyvar*xyvar);
		if(m4<0.000001)
			m4=0.000001;
		double a11=yvar/m4;
		double a12=xyvar/m4;
		double a22=xvar/m4;
		double tmp=a11-a22;
		if(tmp==0.0)
			tmp=0.000001;
		double theta=0.5*Math.atan(2.0*a12/tmp);
		if(theta<0.0)
			theta+=HALFPI;
		if(a12>0.0)
			theta+=HALFPI;
		else if(a12==0.0){
			if(a22>a11){
				theta=0.0;
				tmp=a22;
				a22=a11;
				a11=tmp;
			}else if(a11!=a22)
				theta=HALFPI;
		}
		tmp=Math.sin(theta);
		if(tmp==0.0)
			tmp=0.000001;
		double z=a12*Math.cos(theta)/tmp;
		double major=Math.sqrt(1.0/Math.abs(a22+z));
		double minor=Math.sqrt(1.0/Math.abs(a11-z));
		double scale=Math.sqrt(npixels/(Math.PI*major*minor)); // equalize
		// areas
		major=major*scale*2.0;
		minor=minor*scale*2.0;
		/*
		 * double angle = 180.0 * theta / Math.PI; if (angle == 180.0) angle =
		 * 0.0;
		 */
		double angle=theta;
		if(angle==Math.PI)
			angle=0.0;
		if(major<minor){
			tmp=major;
			major=minor;
			minor=tmp;
		}
		float[] output={(float)xsum,(float)ysum,(float)angle,(float)major,(float)minor};
		return output;
	}

	public static float get_closest_dist(Polygon poly,float x,float y,boolean interp){
		float mindist2=get_dist2(poly.xpoints[0],poly.ypoints[0],x,y);
		int minindex=0;
		for(int i=1;i<poly.npoints;i++){
			float dist2=get_dist2(poly.xpoints[i],poly.ypoints[i],x,y);
			if(dist2<mindist2){
				minindex=i;
				mindist2=dist2;
			}
		}
		if(!interp)
			return (float)Math.sqrt(mindist2);
		int prev=minindex-1;
		if(prev<0)
			prev=poly.npoints-1;
		int next=minindex+1;
		if(next>=poly.npoints)
			next=0;
		float distprev2=get_dist2(x,y,poly.xpoints[prev],poly.ypoints[prev]);
		float distnext2=get_dist2(x,y,poly.xpoints[next],poly.ypoints[next]);
		if(distnext2<distprev2){
			distprev2=distnext2;
			prev=next;
		}
		return get_closest_dist(poly.xpoints[minindex],poly.ypoints[minindex],poly.xpoints[prev],poly.ypoints[prev],x,y);
	}
	
	/***********
	 * is p3 to the left or right of the p1->p2 line when looking down the line from p1 to p2?
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param x3
	 * @param y3
	 * @return
	 */
	public static int whichSide(float x1,float y1,float x2,float y2,float x3,float y3){
		float crossprod=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
		if(crossprod>0.0f) return 1;
		else if(crossprod==0.0f) return 0;
		else return -1;
	}

	public static float get_closest_dist(float x1,float y1,float x2,float y2,float x3,float y3){
		// get the closest distance between p3 and the line p1-p2
		// first get the distance of the tangent intersection from x1,y1
		// this is the dot product of the line and the line start to the point
		// divided by the length of the line
		float length=get_dist(x1,y1,x2,y2);
		float interdist=((x2-x1)*(x3-x1)+(y2-y1)*(y3-y1))/length;
		// don't allow the intersection to move outside the line segment
		if(interdist>length)
			interdist=length;
		if(interdist<0.0f)
			interdist=0.0f;
		// now get the tangent intersection
		float xinc=(x2-x1)/length;
		float yinc=(y2-y1)/length;
		float interx=xinc*interdist+x1;
		float intery=yinc*interdist+y1;
		return get_dist(interx,intery,x3,y3);
	}
	
	public static float get_dist(float x1,float y1,float x2,float y2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}
	
	public static float get_dist2(float x1,float y1,float x2,float y2){
		return ((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}

	public static float[] get_max_coords(float[] image,Polygon poly,int width,int height){
		Rectangle r=new Rectangle(0,0,width,height);
		if(poly!=null){
			r=poly.getBounds();
		}
		float max=image[r.x+r.y*width];
		int maxx=r.x;
		int maxy=r.y;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly==null||poly.contains(j,i)){
					int index=j+i*width;
					if(image[index]>max){
						max=image[index];
						maxx=j;
						maxy=i;
					}
				}
			}
		}
		return new float[]{max,maxx,maxy};
	}

}
