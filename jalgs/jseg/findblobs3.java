/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.jstatistics;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

public class findblobs3 implements Runnable{
	// this class finds contiguous blobs using a flood fill mechanism
	// these blobs can then be sorted using area and circularity
	// object images are typically stored as floating bit images with contiguous
	// blobs having the same index
	int maxStackSize=500; // will be increased as needed
	int[] xstack=new int[maxStackSize];
	int[] ystack=new int[maxStackSize];
	int stackSize;
	public int width,height,nobjects;
	public byte[] runin;
	public float[] runout;

	public findblobs3(int width1,int height1){
		width=width1;
		height=height1;
		nobjects=0;
	}

	public void set_objects(float[] objects){
		nobjects=(int)maxarray(objects);
	}

	public float[] dofindblobs(byte[] data1){
		boolean[] data=new boolean[data1.length];
		for(int i=0;i<data1.length;i++){
			data[i]=(data1[i]!=(byte)0);
		}
		float[] counter=new float[data1.length];
		// raster scan the image to find nonzero values
		int id=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(data[j+i*width]){
					id++;
					fillboolean(data,counter,id,j,i);
				}
			}
		}
		nobjects=id;
		return counter;
	}
	
	public void run(){
		float[] temp=dofindblobs(runin);
		System.arraycopy(temp,0,runout,0,width*height);
	}

	public float[] dofindblobs(Object data1,float thresh){
		byte[] bimage=threshimage(data1,thresh);
		return dofindblobs(bimage);
	}
	
	public float[] coords2objects(float[][] coords,float diameter,boolean renumber){
		//here we use coordinates to create a masked image with indexed objects like one would get from thresholding
		//start with the last coordinate because it will be masked by others
		float[] objects=new float[width*height];
		for(int i=0;i<coords.length;i++){
			int id=coords.length-i;
			findblobs.set_circle_val(id,objects,coords[coords.length-1-i][0],coords[coords.length-1-i][1],diameter,width,height);
		}
		nobjects=coords.length;
		//now set any pixels which are neighbored by different values to zero
		separateobjects(objects);
		//and renumber in standard raster discovery order if asked for
		if(renumber) renumber_objects2(objects);
		return objects;
	}
	
	public float[] outlines2objects(Polygon[] polys){
		float[] objects=new float[width*height];
		for(int i=0;i<polys.length;i++) fillPolygon(objects,polys[i],i+1);
		nobjects=polys.length;
		return objects;
	}
	
	public void fillPolygon(float[] objects,Polygon poly,float val){
		Rectangle r=poly.getBounds();
		Rectangle r2=new Rectangle(0,0,width,height);
		for(int j=r.y;j<(r.y+r.height);j++){
			for(int k=r.x;k<(r.x+r.width);k++){
				if(r2.contains(k,j) && poly.contains(k,j)) objects[k+j*width]=val;
			}
		}
	}
	
	public void fillPolygon(byte[] mask,Polygon poly){
		Rectangle r=poly.getBounds();
		Rectangle r2=new Rectangle(0,0,width,height);
		for(int j=r.y;j<(r.y+r.height);j++){
			for(int k=r.x;k<(r.x+r.width);k++){
				if(r2.contains(k,j) && poly.contains(k,j)) mask[k+j*width]=(byte)255;
			}
		}
	}

	public static byte[] threshimage(Object data,float thresh){
		byte[] bimage=null;
		if(data instanceof float[])
			return threshimage((float[])data,thresh);
		else if(data instanceof short[])
			return threshimage((short[])data,thresh);
		else
			return threshimage((byte[])data,thresh);
	}

	public static byte[] threshimage(float[] data,float thresh){
		byte[] data2=new byte[data.length];
		for(int i=0;i<data.length;i++){
			if(data[i]>=thresh)
				data2[i]=(byte)255;
			else
				data2[i]=(byte)0;
		}
		return data2;
	}

	public static byte[] threshimage(short[] data,float thresh){
		byte[] data2=new byte[data.length];
		for(int i=0;i<data.length;i++){
			if((data[i]&0xffff)>=thresh)
				data2[i]=(byte)255;
			else
				data2[i]=(byte)0;
		}
		return data2;
	}

	public static byte[] threshimage(byte[] data,float thresh){
		byte[] data2=new byte[data.length];
		for(int i=0;i<data.length;i++){
			if((data[i]&0xff)>=thresh)
				data2[i]=(byte)255;
			else
				data2[i]=(byte)0;
		}
		return data2;
	}

	public byte[] tobinary(float[] objects,boolean separate){
		byte[] out=new byte[width*height];
		if(separate){
			float[] objects2=new float[width*height];
			System.arraycopy(objects,0,objects2,0,width*height);
			separateobjects(objects2);
			for(int i=0;i<width*height;i++){
				if(objects2[i]>0.0f){
					out[i]=(byte)255;
				}
			}
		}else{
			for(int i=0;i<width*height;i++){
				if(objects[i]>0.0f){
					out[i]=(byte)255;
				}
			}
		}
		return out;
	}

	public void separateobjects(float[] objects){
		//here we ensure that all objects are separated by black space
		float[] temp=new float[width*height];
		System.arraycopy(objects,0,temp,0,width*height);
		for(int i=1;i<height-1;i++){
			for(int j=1;j<width-1;j++){
				if(objects[j+i*width]>0.0f){
					float value=objects[j+i*width];
					float[] neighbors=getNeighbors(objects,j,i);
					for(int k=0;k<8;k++){
						if(neighbors[k]>0.0f&&neighbors[k]!=value){
							temp[j+i*width]=0.0f;
							break;
						}
					}
				}
			}
		}
		System.arraycopy(temp,0,objects,0,width*height);
	}

	public void separateobjects(float[] objects,Polygon line,boolean close){
		// here we use a polyline to separate two or more objects and then
		// refind the objects
		clear_polyline(objects,line,close);
		byte[] temp=tobinary(objects,true);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	public void separateobjects(float[] objects,Polygon line,boolean close,int id){
		// here we use a polyline to separate two or more objects and then
		// refind the objects
		clear_polyline(objects,line,close,id);
		byte[] temp=tobinary(objects,true);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	public void clear_polyline(float[] objects,Polygon line,boolean close){
		int npts=line.npoints;
		int[] xpts=line.xpoints;
		int[] ypts=line.ypoints;
		for(int i=1;i<npts;i++){
			clear_line(objects,xpts[i-1],ypts[i-1],xpts[i],ypts[i]);
		}
		if(close){
			clear_line(objects,xpts[npts-1],ypts[npts-1],xpts[0],ypts[0]);
		}
	}

	public void clear_polyline(float[] objects,Polygon line,boolean close,int id){
		int npts=line.npoints;
		int[] xpts=line.xpoints;
		int[] ypts=line.ypoints;
		for(int i=1;i<npts;i++){
			clear_line(objects,xpts[i-1],ypts[i-1],xpts[i],ypts[i],id);
		}
		if(close){
			clear_line(objects,xpts[npts-1],ypts[npts-1],xpts[0],ypts[0],id);
		}
	}

	public void clear_line(float[] objects,int x1,int y1,int x2,int y2){
		double length=Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		double xinc=(x2-x1)/length;
		double yinc=(y2-y1)/length;
		double x=x1;
		double y=y1;
		clear_dot(objects,(int)x,(int)y);
		for(int i=1;i<=(int)length;i++){
			x+=xinc;
			y+=yinc;
			clear_dot(objects,(int)x,(int)y);
		}
	}

	public void clear_line(float[] objects,int x1,int y1,int x2,int y2,int id){
		double length=Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		double xinc=(x2-x1)/length;
		double yinc=(y2-y1)/length;
		double x=x1;
		double y=y1;
		clear_dot(objects,(int)x,(int)y,id);
		for(int i=1;i<=(int)length;i++){
			x+=xinc;
			y+=yinc;
			clear_dot(objects,(int)x,(int)y,id);
		}
	}

	public void clear_dot(float[] objects,int x,int y){
		if(x<1||x>=(width-1)||y<1||y>=(height-1)){
			return;
		}
		objects[x+y*width]=0.0f;
		objects[x+1+y*width]=0.0f;
		objects[x+(y+1)*width]=0.0f;
		objects[x+1+(y+1)*width]=0.0f;
	}

	public void clear_dot(float[] objects,int x,int y,int id){
		if(x<1||x>=(width-1)||y<1||y>=(height-1)){
			return;
		}
		if(objects[x+y*width]==id)
			objects[x+y*width]=0.0f;
		if(objects[x+1+y*width]==id)
			objects[x+1+y*width]=0.0f;
		if(objects[x+(y+1)*width]==id)
			objects[x+(y+1)*width]=0.0f;
		if(objects[x+1+(y+1)*width]==id)
			objects[x+1+(y+1)*width]=0.0f;
	}

	public int[] overlap_objects(float[] obj1,float[] obj2){
		int nobj1=(int)maxarray(obj1);
		int nobj2=(int)maxarray(obj2);
		boolean[] stats=new boolean[nobj1];
		for(int i=0;i<width*height;i++){
			if(obj1[i]>0f&&obj2[i]>0f){
				stats[(int)obj1[i]-1]=true;
			}
		}
		int overlaps=0;
		for(int i=0;i<nobj1;i++){
			if(stats[i]){
				overlaps++;
			}
		}
		int[] stats2={overlaps,nobj1-overlaps,nobj2-overlaps};
		return stats2;
	}

	/*********************
	 * dialates objects without letting them merge
	 * @param objects: float indexed objects image
	 */
	public void dilateobjects(float[] objects){
		dilateobjects(objects,true);
	}

	/********************
	 * dialates objects without letting them merge
	 * @param objects: float indexed objects image
	 * @param renumber: whether or not to renumber the objects
	 */
	public void dilateobjects(float[] objects,boolean renumber){
		float[] dilated=new float[objects.length];
		System.arraycopy(objects,0,dilated,0,objects.length);
		// here, if a pixel is neighbored by an object pixel it is added to that
		// object
		// if it is neighbored by two different objects it stays black
		// this means that objects will always stay separate
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				if(objects[j+i*width]==0.0f){
					float[] objarray=getNeighbors(objects,j,i);
					float objid=0.0f;
					for(int k=0;k<8;k++){
						if(objarray[k]>0.0f){
							objid=objarray[k];
							break;
						}
					}
					if(objid>0.0f){
						for(int k=0;k<8;k++){
							if(objarray[k]>0.0f&&objarray[k]!=objid){
								objid=0.0f;
								break;
							}
						}
					}
					dilated[j+i*width]=objid;
				}
			}
		}
		separateobjects(dilated);
		// note we have to refind objects because they may have shifted relative
		// position
		System.arraycopy(dilated,0,objects,0,width*height);
		if(renumber) renumber_objects(objects);
	}

	public void dilateobjects2(float[] objects){
		// here, if a pixel is neighbored by an object pixel it is added to that
		// object
		// objects can blend together
		byte[] temp=tobinary(objects,false);
		(new binary_processing(width,height)).dilate(temp);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	public void erodeobjects(float[] objects){
		float[] eroded=new float[objects.length];
		System.arraycopy(objects,0,eroded,0,objects.length);
		// here, if a pixel is neighbored by a black pixel it is removed from
		// that object
		// if it is neighbored by another object, it is removed as well
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				if(objects[j+i*width]>0.0f){
					float[] objarray=getNeighbors(objects,j,i);
					for(int k=0;k<8;k++){
						if(objarray[k]==0.0f||objarray[k]!=objects[j+i*width]){
							eroded[j+i*width]=0.0f;
						}
					}
				}
			}
		}
		byte[] temp=tobinary(eroded,true);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	public void erodeobjects2(float[] objects){
		float[] eroded=objects.clone();
		float[] skeleton=objects.clone();
		skeletonize(skeleton);
		// here, if a pixel is neighbored by a black pixel it is removed from
		// that object
		// if it is neighbored by another object, it is removed as well
		// this version does not allow an object to be split into two objects
		// this is done by retaining a skeleton of the original image
		// irrespective of the erosion
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				if(objects[j+i*width]>0.0f){
					float[] objarray=getNeighbors(objects,j,i);
					for(int k=0;k<8;k++){
						if(objarray[k]==0.0f||objarray[k]!=objects[j+i*width]){
							eroded[j+i*width]=0.0f;
						}
					}
				}
			}
		}
		for(int i=0;i<width*height;i++){
			if(skeleton[i]>0.0f){
				eroded[i]=skeleton[i];
			}
		}
		// again have to refind objects because of changing positions
		byte[] temp=tobinary(eroded,true);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	public void erodeobjects3(float[] objects){
		float[] eroded=objects.clone();
		// here, if a pixel is neighbored by a black pixel it is removed from
		// that object
		// if it is neighbored by another object, it is removed as well
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				if(objects[j+i*width]>0.0f){
					float[] objarray=getNeighbors(objects,j,i);
					for(int k=0;k<8;k++){
						if(objarray[k]==0.0f||objarray[k]!=objects[j+i*width]){
							eroded[j+i*width]=0.0f;
						}
					}
				}
			}
		}
		// here we don't refind the objects so that numbering is maintained
		System.arraycopy(eroded,0,objects,0,width*height);
	}
	
	public float[] binobjects(float[] objects,int binby,boolean renumber){
		//here we downscale an objects image
		//binned regions containing one object will be assigned to that object
		//binned regions containing two or more objects will be set to zero
		int newwidth=(int)((float)width/(float)binby);
		int newheight=(int)((float)height/(float)binby);
		float[] newobjects=new float[newwidth*newheight];
		for(int i=0;i<newheight;i++){
			for(int j=0;j<newwidth;j++){
				float id=0.0f;
				boolean found1=false;
				boolean found2=false;
				for(int k=0;k<binby;k++){
					for(int l=0;l<binby;l++){
						float val=objects[j*binby+l+(i*binby+k)*width];
						if(val>0.0f){
							if(!found1){
								id=val;
								found1=true;
							} else {
								if(val!=id) found2=true;
							}
						}
					}
				}
				if(found1 && !found2) newobjects[j+i*newwidth]=id;
			}
		}
		width=newwidth;
		height=newheight;
		if(renumber) renumber_objects2(newobjects);
		return newobjects;
	}
	
	/****************
	 * enlarges the image maintaining object boundaries
	 * @param objects: the objects images
	 * @param enlargeby: the scaling factor
	 */
	public float[] enlarge(float[] objects,int enlargeby){
		int newwidth=width*enlargeby;
		int newheight=height*enlargeby;
		float[] newobjects=new float[newwidth*newheight];
		for(int i=0;i<newheight;i++){
			int oldy=(int)((float)i/(float)enlargeby);
			for(int j=0;j<newwidth;j++){
				int oldx=(int)((float)j/(float)enlargeby);
				newobjects[j+i*newwidth]=objects[oldx+width*oldy];
			}
		}
		width=newwidth;
		height=newheight;
		return newobjects;
	}

	public void closeobjects(float[] objects){
		// here we dilate, then rerecognize objects, then erode retaining a
		// skeleton
		// this will combine objects that dilate into each other
		dilateobjects2(objects);
		erodeobjects2(objects);
	}

	public void clear_edges(float[] objects){
		clear_edges(objects,false);
	}

	public void clear_edges(float[] objects,boolean renumber){
		for(int i=0;i<width;i++){
			if(objects[i]>0.0f){
				delete_object(objects,objects[i]);
			}
			if(objects[i+(height-1)*width]>0.0f){
				delete_object(objects,objects[i+(height-1)*width]);
			}
		}
		for(int i=0;i<height;i++){
			if(objects[i*width]>0.0f){
				delete_object(objects,objects[i*width]);
			}
			if(objects[i*width+width-1]>0.0f){
				delete_object(objects,objects[i*width+width-1]);
			}
		}
		if(renumber)
			renumber_objects(objects);
	}

	public void clear_borders(float[] objects,int border){
		for(int i=0;i<width;i++){
			for(int j=0;j<border;j++){
				if(objects[j*width+i]>0.0f){
					delete_object(objects,objects[j*width+i]);
				}
			}
			for(int j=(height-border);j<height;j++){
				if(objects[j*width+i]>0.0f){
					delete_object(objects,objects[j*width+i]);
				}
			}
		}
		for(int i=border;i<(height-border);i++){
			for(int j=0;j<border;j++){
				if(objects[j+i*width]>0.0f){
					delete_object(objects,objects[j+i*width]);
				}
			}
			for(int j=(width-border);j<width;j++){
				if(objects[j+i*width]>0.0f){
					delete_object(objects,objects[j+i*width]);
				}
			}
		}
	}

	public void delete_object(float[] objects,float id){
		nobjects--;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(objects[j+i*width]==id){
					objects[j+i*width]=0.0f;
				}
			}
		}
	}

	public void delete_object(float[] objects,float id,Rectangle r){
		nobjects--;
		for(int j=r.y;j<(r.y+r.height);j++){
			for(int k=r.x;k<(r.x+r.width);k++){
				if(objects[k+j*width]==id){
					objects[k+j*width]=0.0f;
				}
			}
		}
	}

	public void copy_object(float[] source,float srcid,float[] newobjects,float newid){
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(source[j+i*width]==srcid){
					newobjects[j+i*width]=newid;
				}
			}
		}
	}

	public void add_object(float[] objects,Polygon poly){
		nobjects++;
		float maxid=maxarray(objects);
		fill_poly_id(objects,maxid+1.0f,poly);
		renumber_objects2(objects);
	}

	public void fill_poly_id(float[] objects,float id,Polygon poly){
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(poly.contains(j,i)){
					objects[j+i*width]=id;
				}
			}
		}
	}

	public void change_id(float[] objects,float oldid,float newid){
		for(int i=0;i<objects.length;i++){
			if(objects[i]==oldid){
				objects[i]=newid;
			}
		}
	}

	public void combine_objects(float[] objects,Polygon poly){
		// first find out which two objects are in the polygon
		Rectangle bounds=poly.getBounds();
		int objcounter=0;
		float[] ids=new float[nobjects];
		for(int i=bounds.y;i<=(bounds.y+bounds.height);i++){
			for(int j=bounds.x;j<=(bounds.x+bounds.width);j++){
				float tempid=objects[j+i*width];
				if(poly.contains(j,i)&&tempid>0.0f){
					boolean found=false;
					for(int k=0;k<objcounter;k++){
						if(tempid==ids[k]){
							found=true;
							break;
						}
					}
					if(!found){
						ids[objcounter]=tempid;
						objcounter++;
					}
				}
			}
		}
		if(objcounter>1){
			for(int i=1;i<objcounter;i++){
				change_id(objects,ids[i],ids[0]);
			}
			fill_poly_id(objects,ids[0],poly);
			renumber_objects(objects);
		}
	}

	/*************
	 * adds the polygon region to the object
	 * @param objects
	 * @param poly
	 */
	public void expand_object(float[] objects,Polygon poly){
		// first find out which object is in the polygon
		Rectangle bounds=poly.getBounds();
		float obj1=0.0f;
		findobj1: for(int i=bounds.y;i<=(bounds.y+bounds.height);i++){
			for(int j=bounds.x;j<=(bounds.x+bounds.width);j++){
				if(poly.contains(j,i)&&objects[j+i*width]>0.0f){
					obj1=objects[j+i*width];
					break findobj1;
				}
			}
		}
		if(obj1==0.0f){
			return;
		}
		fill_poly_id(objects,obj1,poly);
	}

	/**********
	 * here we replace the object partially contained in the polygon with the polygon itself
	 * @param objects
	 * @param poly
	 */
	public void edit_object(float[] objects,Polygon poly){
		// first find out which object is in the polygon
		Rectangle bounds=poly.getBounds();
		float obj1=0.0f;
		findobj1: for(int i=bounds.y;i<=(bounds.y+bounds.height);i++){
			for(int j=bounds.x;j<=(bounds.x+bounds.width);j++){
				if(poly.contains(j,i)&&objects[j+i*width]>0.0f){
					obj1=objects[j+i*width];
					break findobj1;
				}
			}
		}
		if(obj1==0.0f){
			return;
		}
		delete_object(objects,obj1);
		add_object(objects,poly);
	}

	/***********
	 * this implementation doesn't necessarily preserve discovery order but is more efficient
	 * @param objects
	 */
	public void renumber_objects(float[] objects){
		//build a table with the new numbers
		int maxid=(int)maxarray(objects);
		int[] hist=new int[maxid+1];
		for(int i=0;i<objects.length;i++){
			if(objects[i]>0.0f) hist[(int)objects[i]]++;
		}
		int[] newnumbers=new int[maxid+1];
		int index=1;
		for(int i=1;i<hist.length;i++){
			if(hist[i]>0){
				newnumbers[i]=index;
				index++;
			}
		}
		nobjects=index-1;
		for(int i=0;i<objects.length;i++){
			if(objects[i]>0.0f){
				int temp=(int)objects[i];
				objects[i]=newnumbers[temp];
			}
		}
	}
	
	/****************
	 * this version recalculates discovery order but is less efficient
	 * it also separates objects in the process
	 * @param objects
	 */
	public void renumber_objects2(float[] objects){
		//this is a less efficient implementation
		byte[] temp=tobinary(objects,true);
		float[] temp2=dofindblobs(temp);
		System.arraycopy(temp2,0,objects,0,width*height);
	}

	private float maxarray(float[] input){
		float max=input[0];
		for(int i=1;i<input.length;i++){
			if(input[i]>max){
				max=input[i];
			}
		}
		return max;
	}

	public void dilate_all_lines(int length,byte[] image){
		byte[] dilated=new byte[image.length];
		System.arraycopy(image,0,dilated,0,image.length);
		for(int i=0;i<8;i++){
			dilate_lines(i,length,dilated);
		}
		System.arraycopy(dilated,0,image,0,dilated.length);
	}

	public void dilate_lines(int directionindex,int length,byte[] image){
		// here we dilate if length pixels along direction index are selected
		// direction indices are e,ne,n,nw,w,sw,s,se
		byte[] dilated=new byte[image.length];
		System.arraycopy(image,0,dilated,0,image.length);
		for(int i=length;i<(height-length-1);i++){
			for(int j=length;j<(width-length-1);j++){
				if(image[j+i*width]==(byte)0){
					boolean test=true;
					for(int k=1;k<=length;k++){
						switch(directionindex){
						case 0:
							test=image[j+i*width+k]!=(byte)0;
							break;
						case 1:
							test=image[j+(i-k)*width+k]!=(byte)0;
							break;
						case 2:
							test=image[j+(i-k)*width]!=(byte)0;
							break;
						case 3:
							test=image[j+(i-k)*width-k]!=(byte)0;
							break;
						case 4:
							test=image[j+i*width-k]!=(byte)0;
							break;
						case 5:
							test=image[j+(i+k)*width-k]!=(byte)0;
							break;
						case 6:
							test=image[j+(i+k)*width]!=(byte)0;
							break;
						case 7:
							test=image[j+(i+k)*width+k]!=(byte)0;
							break;
						}
						if(!test){
							break;
						}
					}
					if(test){
						dilated[j+i*width]=(byte)255;
					}
				}else{
					dilated[j+i*width]=(byte)255;
				}
			}
		}
		System.arraycopy(dilated,0,image,0,dilated.length);
	}

	// Does a 4-connected flood fill copied from ImageJ--modified to fill true
	// with false
	// the object gets deleted from data
	public boolean fillboolean(boolean[] data,float[] counter,int id,int x,int y){
		stackSize=0;
		push(x,y);
		while(true){
			x=popx();
			if(x==-1)
				return true;
			y=popy();
			if(!data[x+y*width])
				continue;
			int x1=x;
			int x2=x;
			while(x1>=0&&data[x1+y*width])
				x1--; // find start of scan-line
			x1++;
			while(x2<width&&data[x2+y*width])
				x2++; // find end of scan-line
			x2--;
			fillLine(data,counter,id,x1,x2,y); // fill scan-line
			boolean inScanLine=false;
			if(y>0){
				for(int i=x1;i<=x2;i++){ // find contiguous scan-lines
					// above this one
					if(!inScanLine&&data[i+(y-1)*width]){
						push(i,y-1);
						inScanLine=true;
					}else if(inScanLine&&(!data[i+(y-1)*width]))
						inScanLine=false;
				}
			}
			inScanLine=false;
			if(y<(height-1)){
				for(int i=x1;i<=x2;i++){ // find contiguous scan-lines
					// below this one
					if(!inScanLine&&data[i+(y+1)*width]){
						push(i,y+1);
						inScanLine=true;
					}else if(inScanLine&&(!data[i+(y+1)*width]))
						inScanLine=false;
				}
			}
		}
	}
	
	/*****************
	 * gets all of the fill limits for the objects (xmin,xmax,ymin,ymax)
	 * @param objects
	 * @return
	 */
	public int[][] getallfilllimits(float[] objects){
		//int nobjects=(int)maxarray(objects);
		int[][] lims=new int[nobjects][4];
		for(int i=0;i<nobjects;i++){
			lims[i][0]=width-1; lims[i][1]=0; lims[i][2]=height-1; lims[i][3]=0;
		}
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(objects[j+i*width]>0.0f){
					int id=(int)objects[j+i*width]-1;
					if(j<lims[id][0]) lims[id][0]=j;
					if(j>lims[id][1]) lims[id][1]=j;
					if(i<lims[id][2]) lims[id][2]=i;
					if(i>lims[id][3]) lims[id][3]=i;
				}
			}
		}
		for(int i=0;i<nobjects;i++){
			if(lims[i][1]<lims[i][0]){
				//here we never found an object, so reset full limits
				lims[i][0]=0; lims[i][1]=width-1; lims[i][2]=0; lims[i][3]=height-1;
			}
		}
		return lims;
	}

	/************
	 * here we flood fill an object just to find its extent
	 * @param objects
	 * @param id
	 * @param x
	 * @param y
	 * @return
	 */
	public int[] getfilllimits(float[] objects,int id,int x,int y){
		float fid=id;
		float[] data=objects.clone();
		stackSize=0;
		push(x,y);
		int[] limits={x,x,y,y};
		while(true){
			x=popx();
			if(x==-1)
				return limits;
			y=popy();
			if(y<limits[2])
				limits[2]=y;
			if(y>limits[3])
				limits[3]=y;
			int x1=x;
			int x2=x;
			while(x1>=0&&data[x1+y*width]==fid)
				x1--;
			x1++;
			if(x1<limits[0])
				limits[0]=x1;
			while(x2<width&&data[x2+y*width]==fid)
				x2++;
			x2--;
			if(x2>limits[1])
				limits[1]=x2;
			fillLine(data,x1,x2,y);
			boolean inScanLine=false;
			if(y>0){
				for(int i=x1;i<=x2;i++){ // find contiguous scan-lines
					// above this one
					if(!inScanLine&&data[i+(y-1)*width]==fid){
						push(i,y-1);
						inScanLine=true;
					}else if(inScanLine&&(data[i+(y-1)*width]!=fid))
						inScanLine=false;
				}
			}
			inScanLine=false;
			if(y<(height-1)){
				for(int i=x1;i<=x2;i++){ // find contiguous scan-lines
					// below this one
					if(!inScanLine&&data[i+(y+1)*width]==fid){
						push(i,y+1);
						inScanLine=true;
					}else if(inScanLine&&(data[i+(y+1)*width]!=fid))
						inScanLine=false;
				}
			}
		}
	}

	final void push(int x,int y){
		stackSize++;
		if(stackSize==maxStackSize){
			int[] newXStack=new int[maxStackSize*2];
			int[] newYStack=new int[maxStackSize*2];
			System.arraycopy(xstack,0,newXStack,0,maxStackSize);
			System.arraycopy(ystack,0,newYStack,0,maxStackSize);
			xstack=newXStack;
			ystack=newYStack;
			maxStackSize*=2;
		}
		xstack[stackSize-1]=x;
		ystack[stackSize-1]=y;
	}

	final int popx(){
		if(stackSize==0)
			return -1;
		else
			return xstack[stackSize-1];
	}

	final int popy(){
		int value=ystack[stackSize-1];
		stackSize--;
		return value;
	}

	final void fillLine(boolean[] data,float[] counter,int id,int x1,int x2,int y){
		if(x1>x2){
			int t=x1;
			x1=x2;
			x2=t;
		}
		for(int x=x1;x<=x2;x++){
			data[x+y*width]=false;
			counter[x+y*width]=id;
		}
	}

	final void fillLine(float[] data,int x1,int x2,int y){
		if(x1>x2){
			int t=x1;
			x1=x2;
			x2=t;
		}
		for(int x=x1;x<=x2;x++){
			data[x+y*width]=0.0f;
		}
	}

	public void fill_holes(float[] objects){
		fill_holes(objects,get_object_outlines(objects));
	}

	public void fill_holes(float[] objects,Polygon[] outlines){
		for(int i=0;i<outlines.length;i++){
			float id=i+1;
			Rectangle r=outlines[i].getBounds();
			for(int j=r.y;j<(r.y+r.height);j++){
				for(int k=r.x;k<(r.x+r.width);k++){
					if(outlines[i].contains(k,j)){
						if(objects[k+j*width]==0.0f){
							objects[k+j*width]=id;
						}
					}
				}
			}
		}
	}

	public double get_perimeter(Polygon poly){
		int polypts=poly.npoints;
		int[] xpts=poly.xpoints;
		int[] ypts=poly.ypoints;
		double distance=0.0;
		for(int i=1;i<polypts;i++){
			distance+=Math.sqrt((xpts[i]-xpts[i-1])*(xpts[i]-xpts[i-1])+(ypts[i]-ypts[i-1])*(ypts[i]-ypts[i-1]));
		}
		distance+=Math.sqrt((xpts[0]-xpts[polypts-1])*(xpts[0]-xpts[polypts-1])+(ypts[0]-ypts[polypts-1])*(ypts[0]-ypts[polypts-1]));
		return distance;
	}
	
	public double get_smoothed_perimeter(Polygon poly){
		//adapted from the ImageJ PolygonRoi class
		int nPoints=poly.npoints;
		int[] xp=poly.xpoints;
		int[] yp=poly.ypoints;
		double length=0.0;
		double dx = (xp[0]+xp[1]+xp[2])/3.0-xp[0];
        double dy = (yp[0]+yp[1]+yp[2])/3.0-yp[0];
        length += Math.sqrt(dx*dx+dy*dy);
        for (int i=1; i<nPoints-2; i++) {
            dx = (xp[i+2]-xp[i-1])/3.0; // = (x[i]+x[i+1]+x[i+2])/3-(x[i-1]+x[i]+x[i+1])/3
            dy = (yp[i+2]-yp[i-1])/3.0; // = (y[i]+y[i+1]+y[i+2])/3-(y[i-1]+y[i]+y[i+1])/3
            length += Math.sqrt(dx*dx+dy*dy);
        }
        dx = xp[nPoints-1]-(xp[nPoints-3]+xp[nPoints-2]+xp[nPoints-1])/3.0;
        dy = yp[nPoints-1]-(yp[nPoints-3]+yp[nPoints-2]+yp[nPoints-1])/3.0;
        length += Math.sqrt(dx*dx+dy*dy);
        dx = xp[nPoints-1]-xp[0];
        dy = yp[nPoints-1]-yp[0];
        length += Math.sqrt(dx*dx+dy*dy);
		return length;
	}

	public void filter_area(float[] objects,int[] arealims){
		filter_area(objects,arealims,false);
	}

	/****************************
	 * filers out objects that are outside the arealims
	 * @param objects: the float indexed objects image
	 * @param arealims: an int array with lower then upper area limits
	 * @param renumber: whether or not to renumber after filtering
	 */
	public void filter_area(float[] objects,int[] arealims,boolean renumber){
		int[] hist=get_areas(objects);
		boolean[] delete=new boolean[hist.length];
		for(int i=0;i<hist.length;i++){
			if(hist[i]<arealims[0]||hist[i]>arealims[1]){
				delete[i]=true;
				//delete_object(objects,i+1);
			}
		}
		for(int i=0;i<width*height;i++){
			if(objects[i]>0.0f){
				int index=(int)objects[i];
				if(delete[index-1]) objects[i]=0.0f;
			}
		}
		if(renumber)
			renumber_objects(objects);
	}

	public void filter_clusters(float[] objects,int minsize,int maxsize,int dilations){
		//clusters are defined as contiguous regions after dilations number of dilations
		//here we delete clusters with less objects than minsize and more objects than maxsize
		float[] newobj=objects.clone();
		for(int i=0;i<dilations;i++){
			dilateobjects2(newobj);
		}
		int[] cluster_count=new int[nobjects];
		float curroldobj=0.0f;
		for(int i=0;i<newobj.length;i++){
			if(newobj[i]>0.0f){
				// we are in a cluster
				if(objects[i]>curroldobj){
					cluster_count[(int)newobj[i]-1]++;
					curroldobj=objects[i];
				}
			}
		}
		// filter the clusters
		for(int i=0;i<newobj.length;i++){
			if(newobj[i]>0.0f){
				int nclusters=cluster_count[(int)newobj[i]-1];
				if(nclusters<minsize||nclusters>maxsize){
					objects[i]=0.0f;
				}
			}
		}
	}

	/*****************
	 * here we filter for area and circularity
	 * @param objects: the indexed objects image
	 * @param outlines: the object outlines
	 * @param limits: minarea,maxarea,mincirc,maxcirc
	 */
	public void filter_area_circ(float[] objects,Polygon[] outlines,float[] limits){
		int minarea=(int)limits[0];
		int maxarea=(int)limits[1];
		if(maxarea==0){
			maxarea=width*height;
		}
		float mincirc=limits[2];
		float maxcirc=limits[3];
		for(int i=0;i<outlines.length;i++){
			float id=i+1;
			int area=0;
			Rectangle r=outlines[i].getBounds();
			for(int j=r.y;j<(r.y+r.height);j++){
				for(int k=r.x;k<(r.x+r.width);k++){
					if(objects[k+j*width]==id){
						area++;
					}
				}
			}
			float perimeter=(float)get_smoothed_perimeter(outlines[i]);
			float circ=4.0f*(float)Math.PI*(area/(perimeter*perimeter));
			if(area<minarea||area>maxarea||circ<mincirc||circ>maxcirc){
				delete_object(objects,id,r);
			}
		}
	}

	public void filter_area_circ(float[] objects,float[] limits){
		filter_area_circ(objects,get_object_outlines(objects),limits);
	}

	/***************
	 * this plugin uses dilation to generate a border around each object
	 * @param objects: the objects
	 * @param circrad: the width of the border
	 * @return
	 */
	public float[] get_circ(float[] objects,int circrad){
		float[] temp=objects.clone();
		//start by dilating the objects (without touching or renumbering)
		for(int i=0;i<circrad;i++)
			dilateobjects(temp,false);
		//now delete the original object
		for(int i=0;i<temp.length;i++){
			if(objects[i]>0.0f)
				temp[i]=0.0f;
		}
		return temp;
	}
	
	/**********************
	 * this plugin uses dilation to generate a border around each object with a gap in between
	 * @param objects: the objects
	 * @param circrad: the border width
	 * @param circgap: the gap width
	 * @return
	 */
	public float[] get_circ(float[] objects,int circrad,int circgap){
		float[] temp=objects.clone();
		for(int i=0;i<circgap;i++)
			dilateobjects(temp,false);
		return get_circ(temp,circrad);
	}

	/************
	 * gets the areas (in pixels) of all of the objects
	 * @param objects
	 * @return
	 */
	public int[] get_areas(float[] objects){
		int nblobs=get_nblobs(objects);
		int[] hist=new int[nblobs];
		for(int i=0;i<width*height;i++){
			if(objects[i]>0.0f){
				hist[(int)objects[i]-1]++;
			}
		}
		return hist;
	}
	
	/*****************
	 * gets the areas, perimeters, and circularities of all of the objects
	 * @param objects
	 * @return
	 */
	public float[][] get_area_perim_circ(float[] objects){
		//get the outlines and then call the more specific method
		return get_area_perim_circ(objects,get_object_outlines(objects));
	}
	
	/************
	 * gets the areas, perimeters, and circularities (1=round, 0=linear) of all of the objects
	 * @param objects
	 * @param outlines: polygon outlines of the objects
	 * @return
	 */
	public float[][] get_area_perim_circ(float[] objects,Polygon[] outlines){
		//returns three arrays, the first with areas and the second with perimeters, then circularities
		float[][] retvals=new float[3][outlines.length];
		for(int i=0;i<outlines.length;i++){
			float id=i+1;
			int area=0;
			Rectangle r=outlines[i].getBounds();
			for(int j=r.y;j<(r.y+r.height);j++){
				for(int k=r.x;k<(r.x+r.width);k++){
					if(objects[k+j*width]==id){
						area++;
					}
				}
			}
			float perimeter=(float)get_smoothed_perimeter(outlines[i]);
			float circ=4.0f*(float)Math.PI*(area/(perimeter*perimeter));
			retvals[0][i]=(float)area;
			retvals[1][i]=perimeter;
			retvals[2][i]=circ;
		}
		return retvals;
	}

	/*************
	 * gets the "stat" measurement for object with id value
	 * @param objects: the objects
	 * @param id: the id of the selected object
	 * @param measurement: the measurement image
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float get_object_stats(float[] objects,int id,Object measurement,String stat){
		return jstatistics.getstatistic(stat,measurement,width,height,get_object_mask(objects,id),null);
	}

	/*************
	 * gets the "stat" measurement for object with id value over an array of measurement images
	 * @param objects: the objects
	 * @param id: the id of the selected object
	 * @param measurement: the measurement image array
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float[] get_object_stats(float[] objects,int id,Object[] measurement,String stat){
		boolean[] mask=get_object_mask(objects,id);
		float[] stats=new float[measurement.length];
		for(int i=0;i<measurement.length;i++){
			stats[i]=jstatistics.getstatistic(stat,measurement[i],width,height,mask,null);
		}
		return stats;
	}
	
	/**************
	 * gets the "stat" measurement for object with id value over a measurement image
	 * @param objects: the objects
	 * @param id: the id of the selected object
	 * @param measurement: the measurement image
	 * @param lims: an array of integers with the boundaries for the selected object
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float get_object_stats(float[] objects,int id,Object measurement,int[] lims,String stat){
		return jstatistics.getstatistic(stat,measurement,width,height,get_object_mask(objects,id,lims),lims,null);
	}

	/**************
	 * gets the "stat" measurement for object with id value over an array of measurement images
	 * @param objects: the objects
	 * @param id: the id of the selected object
	 * @param measurement: the measurement image
	 * @param lims: an array of integers with the boundaries for the selected object
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float[] get_object_stats(float[] objects,int id,Object[] measurement,int[] lims,String stat){
		boolean[] mask=get_object_mask(objects,id,lims);
		float[] stats=new float[measurement.length];
		for(int i=0;i<measurement.length;i++){
			stats[i]=jstatistics.getstatistic(stat,measurement[i],width,height,mask,lims,null);
		}
		return stats;
	}

	/**************
	 * gets all of the object stats for an array of measurement images
	 * @param objects: the objects
	 * @param measurement: the measurement image array
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float[][] get_all_object_stats(float[] objects,Object[] measurement,String stat){
		float[][] allstats=new float[nobjects][];
		for(int i=0;i<nobjects;i++)
			allstats[i]=get_object_stats(objects,i+1,measurement,stat);
		return allstats;
	}

	/**************
	 * gets all of the object stats for a measurement image
	 * @param objects: the objects
	 * @param measurement: the measurement image
	 * @param stat: the statistic to measure
	 * @return
	 */
	public float[] get_all_object_stats(float[] objects,Object measurement,String stat){
		float[] allstats=new float[nobjects];
		for(int i=0;i<nobjects;i++)
			allstats[i]=get_object_stats(objects,i+1,measurement,stat);
		return allstats;
	}
	
	/**************
	 * gets all of the object stats for an array of measurement images
	 * @param objects: the objects
	 * @param measurement: the measurement image array
	 * @param lims: the boundaries of the objects
	 * @param stat: the measurement statistic
	 * @return
	 */
	public float[][] get_all_object_stats(float[] objects,Object[] measurement,int[][] lims,String stat){
		float[][] allstats=new float[nobjects][];
		for(int i=0;i<nobjects;i++)
			allstats[i]=get_object_stats(objects,i+1,measurement,lims[i],stat);
		return allstats;
	}
	
	/***********
	 * gets all of the object stats for a measurement image
	 * @param objects: the objects
	 * @param measurement: the measurement image
	 * @param lims: the boundaries of the objects
	 * @param stat: the measurement statistic
	 * @return
	 */
	public float[] get_all_object_stats(float[] objects,Object measurement,int[][] lims,String stat){
		float[] allstats=new float[nobjects];
		for(int i=0;i<nobjects;i++)
			allstats[i]=get_object_stats(objects,i+1,measurement,lims[i],stat);
		return allstats;
	}

	public int[] get_cluster_ids(float[] clusters,float[] objects,int clusterid){
		//here clusters is an image of objects containing "objects"
		//we find the cluster corresponding to clusterid and return its component objects
		List<Integer> list=new ArrayList<Integer>();
		float oldobj=0.0f;
		for(int i=0;i<clusters.length;i++){
			if(clusters[i]==clusterid){
				if(objects[i]>oldobj){
					list.add((int)objects[i]);
					oldobj=objects[i];
				}
			}
		}
		int[] output=new int[list.size()];
		for(int i=0;i<list.size();i++)
			output[i]=list.get(i).intValue();
		return output;
	}
	
	public boolean[] find_clusters(float[] objects){
		//here we find which objects are "clusters"
		//these are non-contiguous objects
		int nobj=(int)maxarray(objects);
		int oldnblobs=nobjects;
		boolean[] iscluster=new boolean[nobj];
		float[] newid=new float[nobj+1];
		float[] renumbered=dofindblobs(objects,0.5f);
		for(int i=0;i<objects.length;i++){
			if(objects[i]>0.0f){
				int id=(int)objects[i];
				if(newid[id]==0.0f) newid[id]=renumbered[i];
				else{
					if(newid[id]!=renumbered[i]) iscluster[id-1]=true;
				}
			}
		}
		nobjects=oldnblobs;
		return iscluster;
	}

	public float[] get_cluster_stats(float[] clusters,float[] objects,int clusterid,Object measurement,String stat){
		int[] clusterids=get_cluster_ids(clusters,objects,clusterid);
		return get_cluster_stats(objects,clusterids,measurement,stat);
	}

	public float[] get_cluster_stats(float[] objects,int[] clusterids,Object measurement,String stat){
		float[] output=new float[clusterids.length];
		for(int i=0;i<clusterids.length;i++){
			float temp=jstatistics.getstatistic(stat,measurement,width,height,get_object_mask(objects,clusterids[i]),null);
			output[i]=temp;
		}
		return output;
	}

	public float[] get_object_spectrum(float[] objects,int id,Object[] measurement,String stat){
		float[] spectrum=new float[measurement.length];
		boolean[] mask=get_object_mask(objects,id);
		for(int i=0;i<measurement.length;i++){
			spectrum[i]=jstatistics.getstatistic(stat,measurement[i],width,height,mask,null);
		}
		return spectrum;
	}

	public boolean[] get_object_mask(float[] objects,int id){
		boolean[] mask=new boolean[width*height];
		for(int i=0;i<width*height;i++){
			if(objects[i]==id){
				mask[i]=true;
			}
		}
		return mask;
	}
	
	public boolean[] get_object_mask(float[] objects,int id,int[] lims){
		boolean[] mask=new boolean[width*height];
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
				if(objects[j+i*width]==id) mask[j+i*width]=true;
			}
		}
		return mask;
	}

	public byte[] get_object_mask2(float[] objects,int id){
		byte[] mask=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(objects[i]==id){
				mask[i]=(byte)255;
			}
		}
		return mask;
	}
	
	public byte[] get_object_mask2(float[] objects,int id,int[] lims){
		byte[] mask=new byte[width*height];
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
				if(objects[j+i*width]==id) mask[j+i*width]=(byte)255;
			}
		}
		return mask;
	}

	public int get_nblobs(float[] objects){
		return (int)maxarray(objects);
	}

	public Polygon[] get_object_outlines(float[] objects){
		int[][] coords=get_objects_coords(objects); //note this renumbers the objects
		return get_object_outlines(objects,coords);
	}

	public Polygon[] get_object_outlines(float[] objects,int[][] coords){
		Polygon[] outlines=new Polygon[coords.length];
		for(int i=0;i<coords.length;i++){
			if(coords[i]!=null) {
				int[][] outline=traceEdge(coords[i][0],coords[i][1],objects);
				outlines[i]=new Polygon(outline[0],outline[1],outline[0].length);
			}
		}
		return outlines;
	}

	public Polygon get_object_outline(float[] objects,int id){
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(objects[j+i*width]==id){
					return get_object_outline(objects,j,i);
				}
			}
		}
		return null;
	}
	
	//returns a list of all of the points in this object
	public List<float[]> getObjectPoints(float[] objects,int[] bounds,int id){
		List<float[]> coords=new ArrayList<float[]>();
		for(int i=bounds[2];i<=bounds[3];i++){
			for(int j=bounds[0];j<=bounds[1];j++){
				if(objects[j+i*width]==(float)id) coords.add(new float[]{j,i});
			}
		}
		return coords;
	}
	
	public List<float[]> getObjectEdgePoints(float[] objects,int[] bounds,int id){
		List<float[]> coords=new ArrayList<float[]>();
		for(int i=bounds[2];i<=bounds[3];i++){
			for(int j=bounds[0];j<=bounds[1];j++){
				float[] neighbors=getNeighbors(objects,j,i);
				int temp=0;
				for(int k=0;k<neighbors.length;k++) if(neighbors[k]!=(float)id) temp++;
				if(temp>0) coords.add(new float[]{j,i});
			}
		}
		return coords;
	}
	
	public Polygon get_cluster_outline(float[] objects,float[] clusters,int clusterid){
		//this combines fragmented clusters by adding an isthmus between them
		int[] ids=get_cluster_ids(clusters,objects,clusterid);
		Polygon[] polys=new Polygon[ids.length];
		for(int i=0;i<ids.length;i++){
			polys[i]=get_object_outline(objects,ids[i]);
		}
		//now we need to combine all of the polygons by adding an isthmus between them
		return combine_polygons(polys);
	}
	
	public float[] link_objects(float[] objects,int[] linkids){
		//here we set all linkids objects equal to the lowest one's id
		//then need to renumber everything carefully afterwards
		int tempnobjects=(int)maxarray(objects);
		int[] hist=new int[tempnobjects+1];
		int minid=(int)jstatistics.getstatistic("Min",linkids,null);
		for(int i=0;i<objects.length;i++){
			if(objects[i]>0.0f){
    			for(int j=0;j<linkids.length;j++){
    				if((int)objects[i]==linkids[j]){
    					objects[i]=(float)minid;
    					break;
    				}
    			}
    			hist[(int)objects[i]]++;
			}
		}
		//now the renumbering
		//build a table with the new numbers
		int[] newnumbers=new int[tempnobjects+1];
		int index=1;
		for(int i=1;i<hist.length;i++){
			if(hist[i]>0){
				newnumbers[i]=index;
				index++;
			}
		}
		nobjects=index-1;
		for(int i=0;i<objects.length;i++){
			if(objects[i]>0.0f){
				int temp=(int)objects[i];
				objects[i]=newnumbers[temp];
			}
		}
		return objects;
	}
	
	public Polygon combine_polygons(Polygon[] polys){
		Polygon poly=new Polygon(polys[0].xpoints.clone(),polys[0].ypoints.clone(),polys[0].npoints);
		for(int i=1;i<polys.length;i++){
			poly=combine_polygons(poly,polys[i]);
		}
		return poly;
	}
	
	public Polygon combine_polygons(Polygon poly1,Polygon poly2){
		//start by finding the closest pair of points
		int index1=0;
		int index2=0;
		float closest2=(poly1.xpoints[0]-poly2.xpoints[0])*(poly1.xpoints[0]-poly2.xpoints[0])+(poly1.ypoints[0]-poly2.ypoints[0])*(poly1.ypoints[0]-poly2.ypoints[0]);
		for(int i=0;i<poly1.npoints;i++){
			for(int j=0;j<poly2.npoints;j++){
				float temp=(poly1.xpoints[i]-poly2.xpoints[j])*(poly1.xpoints[i]-poly2.xpoints[j])+(poly1.ypoints[i]-poly2.ypoints[j])*(poly1.ypoints[i]-poly2.ypoints[j]);
				if(temp<closest2){
					closest2=temp;
					index1=i;
					index2=j;
				}
			}
		}
		int[] newxpts=new int[poly1.npoints+poly2.npoints+2];
		int[] newypts=new int[poly1.npoints+poly2.npoints+2];
		int counter=0;
		for(int i=0;i<=index1;i++){
			newxpts[counter]=poly1.xpoints[i];
			newypts[counter]=poly1.ypoints[i];
			counter++;
		}
		for(int i=index2;i<poly2.npoints;i++){
			newxpts[counter]=poly2.xpoints[i];
			newypts[counter]=poly2.ypoints[i];
			counter++;
		}
		for(int i=0;i<=index2;i++){
			newxpts[counter]=poly2.xpoints[i];
			newypts[counter]=poly2.ypoints[i];
			counter++;
		}
		for(int i=index1;i<poly1.npoints;i++){
			newxpts[counter]=poly1.xpoints[i];
			newypts[counter]=poly1.ypoints[i];
			counter++;
		}
		return new Polygon(newxpts,newypts,newxpts.length);
	}

	public Polygon get_object_outline(float[] objects,int x,int y){
		int[][] outline=traceEdge(x,y,objects);
		return new Polygon(outline[0],outline[1],outline[0].length);
	}

	public int[][] get_objects_coords(float[] objects){
		renumber_objects(objects);
		int[][] coords=new int[nobjects][];
		int counter=0;
		float currobject=0.0f;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				// outline the raster first point of every object
				// raster should find objects in increasing id order
				if(objects[j+i*width]>currobject){
					currobject=objects[j+i*width];
					int[] temp={j,i};
					coords[counter]=temp;
					counter++;
					// if(counter>=nobjects){break searchloop;}
				}
			}
		}
		return coords;
	}

	public int[] get_object_coords(float[] objects,int id){
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(objects[j+i*width]==id){
					int[] temp={j,i};
					return temp;
				}
			}
		}
		return null;
	}

	static int[] table={0,0,0,1,0,0,1,3,0,0,3,1,1,0,1,3,0,0,0,0,0,0,0,0,2,0,2,0,3,0,3,3,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,3,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,0,0,0,0,3,0,0,0,3,0,2,0,0,0,3,1,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,2,3,1,3,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,0,1,0,0,0,1,0,0,0,0,0,0,0,0,3,3,0,1,0,0,0,0,2,2,0,0,2,0,0,0};

	public void skeletonize(float[] image){
		// this algorithm and its functions is modified from the ImageJ version
		// in which in turn was based on
		// Zhang and Suen (CACM, March 1984, 236-239)
		// the algorithm came from the orignal BinaryProcessor ImageJ class
		int pass=0;
		int pixelsRemoved;
		//need to eliminate edge pixels, clearing edge objects can be too dramatic
		//clear_edges(image);
		//try just clearing the edge pixels
		for(int i=0;i<width;i++){image[i]=0.0f; image[i+(height-1)*width]=0.0f;}
		for(int i=0;i<height;i++){image[width*i]=0.0f; image[width*i+width-1]=0.0f;}
		do{
			pixelsRemoved=thin(pass++,table,image);
			pixelsRemoved=thin(pass++,table,image);
		}while(pixelsRemoved>0);
	}

	public int thin(int pass,int[] table,float[] image){
		float p1,p2,p3,p4,p5,p6,p7,p8,p9;
		int inc=height/25;
		if(inc<1)
			inc=1;
		float bgColor=0.0f;
		float[] pixels2=image.clone();
		int index,code;
		float v;
		int offset,rowOffset=width;
		int pixelsRemoved=0;
		for(int y=1;y<(height-1);y++){
			offset=1+y*width;
			for(int x=1;x<(width-1);x++){
				p5=pixels2[offset];
				v=p5;
				if(v!=bgColor){
					p1=pixels2[offset-rowOffset-1];
					p2=pixels2[offset-rowOffset];
					p3=pixels2[offset-rowOffset+1];
					p4=pixels2[offset-1];
					p6=pixels2[offset+1];
					p7=pixels2[offset+rowOffset-1];
					p8=pixels2[offset+rowOffset];
					p9=pixels2[offset+rowOffset+1];
					index=0;
					if(p1==p5)
						index|=1;
					if(p2==p5)
						index|=2;
					if(p3==p5)
						index|=4;
					if(p6==p5)
						index|=8;
					if(p9==p5)
						index|=16;
					if(p8==p5)
						index|=32;
					if(p7==p5)
						index|=64;
					if(p4==p5)
						index|=128;
					code=table[index];
					if((pass&1)==1){ // odd pass
						if(code==2||code==3){
							v=bgColor;
							pixelsRemoved++;
						}
					}else{ // even pass
						if(code==1||code==3){
							v=bgColor;
							pixelsRemoved++;
						}
					}
				}
				image[offset++]=v;
			}
		}
		return pixelsRemoved;
	}

	public static Polygon dilate_polygon(Polygon poly){
		Rectangle r=poly.getBounds();
		float[] temp=new float[(r.width+4)*(r.height+4)];
		int xoff=r.x-2;
		int yoff=r.y-2;
		for(int i=0;i<r.height;i++){
			for(int j=0;j<r.width;j++){
				if(poly.contains(j+xoff,i+yoff)){
					temp[j+2+(i+2)*(r.width+4)]=1.0f;
				}
			}
		}
		findblobs3 fb=new findblobs3(r.width+4,r.height+4);
		fb.dilateobjects(temp);
		Polygon poly2=fb.get_object_outline(temp,1);
		poly2.translate(xoff,yoff);
		return poly2;
	}

	/***************
	 * this version doesn't include the center pixel
	 * @param objects
	 * @param x
	 * @param y
	 * @return
	 */
	public float[] getNeighbors(float[] objects,int x,int y){
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		float[] temp=new float[8];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2+=2;
		temp[4]=objects[temp2];
		temp2+=(width-2);
		temp[5]=objects[temp2];
		temp2++;
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		return temp;
	}

	/************
	 * this version includes the center pixel
	 * @param objects
	 * @param x
	 * @param y
	 * @return
	 */
	public float[] getNeighbors2(float[] objects,int x,int y){
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		float[] temp=new float[9];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[temp2];
		temp2++;
		temp[1]=objects[temp2];
		temp2++;
		temp[2]=objects[temp2];
		temp2+=(width-2);
		temp[3]=objects[temp2];
		temp2++;
		temp[4]=objects[temp2];
		temp2++;
		temp[5]=objects[temp2];
		temp2+=(width-2);
		temp[6]=objects[temp2];
		temp2++;
		temp[7]=objects[temp2];
		temp2++;
		temp[8]=objects[temp2];
		return temp;
	}

	// trace edge algorithm and its intro and functions adapted from the ImageJ
	// Wand class code
	/*
	 * Trace the outline, starting at a point (startX, startY). Pixel (startX-1,
	 * startY) must be outside, (startX, startY) must be inside, or reverse.
	 * Otherwise an endless loop will occur (and eat up all memory). Traces
	 * 8-connected inside pixels unless fourConnected is true. Returns whether
	 * the selection created encloses an 'inside' area and not an inner hole.
	 */
	public int[][] traceEdge(int startX,int startY,float[] data){
		// Let us name the crossings between 4 pixels vertices, then the
		// vertex (x,y) marked with '+', is between pixels (x-1, y-1) and (x,y):
		//
		// pixel x-1 x
		// y-1 |
		// ----+----
		// y |
		//
		// The four principal directions are numbered such that the direction
		// number * 90 degrees gives the angle in the mathematical sense; and
		// the directions to the adjacent pixels (for inside(x,y,direction) are
		// at (number * 90 - 45) degrees:
		// walking pixel
		// directions: 1 directions: 2 | 1
		// 2 + 0 ----+----
		// 3 3 | 0
		//
		// Directions, like angles, are cyclic; direction -1 = direction 3, etc.
		//
		// The algorithm: We walk along the border, from one vertex to the next,
		// with the outside pixels always being at the left-hand side.
		// For 8-connected tracing, we always trying to turn left as much as
		// possible, to encompass an area as large as possible.
		// Thus, when walking in direction 1 (up, -y), we start looking
		// at the pixel in direction 2; if it is inside, we proceed in this
		// direction (left); otherwise we try with direction 1 (up); if pixel 1
		// is not inside, we must proceed in direction 0 (right).
		//
		// 2 | 1 (i=inside, o=outside)
		// direction 2 < ---+---- > direction 0
		// o | i
		// ^ direction 1 = up = starting direction
		//
		// For 4-connected pixels, we try to go right as much as possible:
		// First try with pixel 1; if it is outside we go in direction 0
		// (right).
		// Otherwise, we examine pixel 2; if it is outside, we go in
		// direction 1 (up); otherwise in direction 2 (left).
		//
		// When moving a closed loop, 'direction' gets incremented or
		// decremented
		// by a total of 360 degrees (i.e., 4) for counterclockwise and
		// clockwise
		// loops respectively. As the inside pixels are at the right side, we
		// have
		// got an outline of inner pixels after a cw loop (direction decremented
		// by 4).
		//
		int[][] outarray=new int[2][1000];
		int npoints=0;
		final int startDirection;
		if(inside(startX,startY,data)) // inside at left, outside right
			startDirection=1; // starting in direction 1 = up
		else{
			startDirection=3; // starting in direction 3 = down
			startY++; // continue after the boundary that has direction 3
		}
		int x=startX;
		int y=startY;
		int direction=startDirection;
		do{
			int newDirection;
			newDirection=direction;
			do{
				if(!inside(x,y,newDirection,data))
					break;
				newDirection++;
			}while(newDirection<direction+2);
			newDirection--;
			outarray[0][npoints]=x;
			outarray[1][npoints]=y;
			npoints++;
			if(npoints>=1000){
				outarray[0]=expand_vector_1000(outarray[0]);
				outarray[1]=expand_vector_1000(outarray[1]);
			}
			switch(newDirection&3){ // '& 3' is remainder modulo 4
			case 0:
				x++;
				break;
			case 1:
				y--;
				break;
			case 2:
				x--;
				break;
			case 3:
				y++;
				break;
			}
			direction=newDirection;
		}while(x!=startX||y!=startY||(direction&3)!=startDirection);
		int[][] new_out_array=new int[2][npoints];
		System.arraycopy(outarray[0],0,new_out_array[0],0,npoints);
		System.arraycopy(outarray[1],0,new_out_array[1],0,npoints);
		return new_out_array;
	}

	private int[] expand_vector_1000(int[] vector){
		int[] temp=new int[vector.length+1000];
		System.arraycopy(vector,0,temp,0,vector.length);
		return temp;
	}

	// check pixel at (x,y), whether it is inside traced area
	private boolean inside(int x,int y,float[] data){
		if(x<0||x>=width||y<0||y>=height)
			return false;
		return data[x+width*y]>0.0f;
	}

	// check pixel in a given direction from vertex (x,y)
	private boolean inside(int x,int y,int direction,float[] data){
		switch(direction&3){ // '& 3' is remainder modulo 4
		case 0:
			return inside(x,y,data);
		case 1:
			return inside(x,y-1,data);
		case 2:
			return inside(x-1,y-1,data);
		case 3:
			return inside(x-1,y,data);
		}
		return false; // will never occur, needed for the compiler
	}

}
