/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.*;

import java.util.*;
import java.awt.Polygon;
import java.awt.Rectangle;

public class findblobs3D{
	// this class finds contiguous blobs using a flood fill mechanism
	// these blobs can then be sorted using area and circularity
	// object images are typically stored as floating bit images with contiguous
	// blobs having the same index
	public int width,height,depth,nobjects;
	public findblobs3 fb;
	int maxStackSize=500; // will be increased as needed
	int[][] stack=new int[maxStackSize][3];
	int stackSize;

	public findblobs3D(int width1,int height1,int depth1){
		width=width1;
		height=height1;
		depth=depth1;
		fb=new findblobs3(width,height);
		nobjects=0;
	}

	public float[][] dofindblobs(byte[][] data1){
		float[][] temp=new float[depth][];
		int[] sliceblobs=new int[depth];
		for(int i=0;i<depth;i++){
			temp[i]=fb.dofindblobs(data1[i]);
			sliceblobs[i]=fb.nobjects;
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,id,k,j,i);
					}
				}
			}
		}
		nobjects=id;
		return objects;
	}

	public float[][] dofindblobs(float[][] data1,float thresh){
		float[][] temp=new float[data1.length][];
		int[] sliceblobs=new int[depth];
		for(int i=0;i<data1.length;i++){
			temp[i]=fb.dofindblobs(data1[i],thresh);
			sliceblobs[i]=fb.nobjects;
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,id,k,j,i);
					}
				}
			}
		}
		nobjects=id;
		return objects;
	}

	public byte[][] tobinary(float[][] objects,boolean separate){
		byte[][] out=new byte[depth][width*height];
		for(int i=0;i<depth;i++){
			for(int j=0;j<width*height;j++)
				if(objects[i][j]>0.0f)
					out[i][j]=(byte)255;
		}
		return out;
	}

	public void clear_edges(float[][] objects,boolean renumber){
		for(int i=0;i<width*height;i++){
			// get the top and bottom
			if(objects[0][i]>0.0f){
				delete_object(objects,objects[0][i]);
			}
			if(objects[depth-1][i]>0.0f){
				delete_object(objects,objects[depth-1][i]);
			}
		}
		for(int k=1;k<(depth-1);k++){
			// now the sides
			for(int i=0;i<width;i++){
				if(objects[k][i]>0.0f){
					delete_object(objects,objects[k][i]);
				}
				if(objects[k][i+(height-1)*width]>0.0f){
					delete_object(objects,objects[k][i+(height-1)*width]);
				}
			}
			for(int i=0;i<height;i++){
				if(objects[k][i*width]>0.0f){
					delete_object(objects,objects[k][i*width]);
				}
				if(objects[k][i*width+width-1]>0.0f){
					delete_object(objects,objects[k][i*width+width-1]);
				}
			}
		}
		if(renumber) renumber_objects(objects);
	}
	
	public void erodeobjects(float[][] objects,int threshold){
		//pixels neighbored in 3D by space or another object are deleted
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					float val=cloned[i][k+j*width];
					float[] neighbors=getNeighbors(cloned,k,j,i);
					int count=0;
					for(int l=0;l<neighbors.length;l++) if(neighbors[l]!=val) count++;
					if(count>threshold) objects[i][j*width+k]=0.0f;
				}
			}
		}
		byte[][] temp=tobinary(objects,false);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);
	}
	
	private void copy_objects(float[][] source,float[][] dest){
		for(int i=0;i<source.length;i++){
			System.arraycopy(source[i],0,dest[i],0,source[i].length);
		}
	}

	public void delete_object(float[][] objects,float id){
		nobjects--;
		for(int k=0;k<depth;k++){
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(objects[k][j+i*width]==id){
						objects[k][j+i*width]=0.0f;
					}
				}
			}
		}
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

	private float maxarray(float[][] input){
		float max=maxarray(input[0]);
		for(int i=1;i<input.length;i++){
			float tm=maxarray(input[i]);
			if(tm>max)
				max=tm;
		}
		return max;
	}

	public int get_nblobs(float[][] objects){
		return (int)maxarray(objects);
	}

	public void renumber_objects(float[][] objects){
		int maxid=(int)maxarray(objects);
		boolean[] occupancy=new boolean[maxid+1];
		for(int j=0;j<depth;j++){
			for(int i=0;i<width*height;i++){
				if(objects[j][i]>0.0f){
					occupancy[(int)objects[j][i]]=true;
				}
			}
		}
		int counter=0;
		float[] newids=new float[maxid+1];
		for(int i=1;i<=maxid;i++){
			if(occupancy[i]){
				counter++;
				newids[i]=(float)counter;
			}
		}
		nobjects=counter;
		for(int j=0;j<depth;j++){
			for(int i=0;i<width*height;i++){
				if(objects[j][i]>0.0f){
					objects[j][i]=newids[(int)objects[j][i]];
				}
			}
		}
	}

	public float[] getNeighbors(float[][] objects,int x,int y,int z){
		//this excludes the current pixel
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		if(z==0||z>=(depth-1)){
			return null;
		}
		float[] temp=new float[26];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[z-1][temp2];
		temp2++;
		temp[1]=objects[z-1][temp2];
		temp2++;
		temp[2]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[3]=objects[z-1][temp2];
		temp2++;
		temp[4]=objects[z-1][temp2];
		temp2++;
		temp[5]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[6]=objects[z-1][temp2];
		temp2++;
		temp[7]=objects[z-1][temp2];
		temp2++;
		temp[8]=objects[z-1][temp2];
		temp2=x-1+(y-1)*width;
		temp[9]=objects[z][temp2];
		temp2++;
		temp[10]=objects[z][temp2];
		temp2++;
		temp[11]=objects[z][temp2];
		temp2+=(width-2);
		temp[12]=objects[z][temp2];
		temp2+=2;
		temp[13]=objects[z][temp2];
		temp2+=(width-2);
		temp[14]=objects[z][temp2];
		temp2++;
		temp[15]=objects[z][temp2];
		temp2++;
		temp[16]=objects[z][temp2];
		temp2=x-1+(y-1)*width;
		temp[17]=objects[z+1][temp2];
		temp2++;
		temp[18]=objects[z+1][temp2];
		temp2++;
		temp[19]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[20]=objects[z+1][temp2];
		temp2++;
		temp[21]=objects[z+1][temp2];
		temp2++;
		temp[22]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[23]=objects[z+1][temp2];
		temp2++;
		temp[24]=objects[z+1][temp2];
		temp2++;
		temp[25]=objects[z+1][temp2];
		return temp;
	}

	public float[] getNeighbors2(float[][] objects,int x,int y,int z){
		//here we include the center pixel
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		if(z==0||z>=(depth-1)){
			return null;
		}
		float[] temp=new float[27];
		int temp2=x-1+(y-1)*width;
		temp[0]=objects[z-1][temp2];
		temp2++;
		temp[1]=objects[z-1][temp2];
		temp2++;
		temp[2]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[3]=objects[z-1][temp2];
		temp2++;
		temp[4]=objects[z-1][temp2];
		temp2++;
		temp[5]=objects[z-1][temp2];
		temp2+=(width-2);
		temp[6]=objects[z-1][temp2];
		temp2++;
		temp[7]=objects[z-1][temp2];
		temp2++;
		temp[8]=objects[z-1][temp2];
		temp2=x-1+(y-1)*width;
		temp[9]=objects[z][temp2];
		temp2++;
		temp[10]=objects[z][temp2];
		temp2++;
		temp[11]=objects[z][temp2];
		temp2+=(width-2);
		temp[12]=objects[z][temp2];
		temp2++;
		temp[13]=objects[z][temp2];
		temp2++;
		temp[14]=objects[z][temp2];
		temp2+=(width-2);
		temp[15]=objects[z][temp2];
		temp2++;
		temp[16]=objects[z][temp2];
		temp2++;
		temp[17]=objects[z][temp2];
		temp2=x-1+(y-1)*width;
		temp[18]=objects[z+1][temp2];
		temp2++;
		temp[19]=objects[z+1][temp2];
		temp2++;
		temp[20]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[21]=objects[z+1][temp2];
		temp2++;
		temp[22]=objects[z+1][temp2];
		temp2++;
		temp[23]=objects[z+1][temp2];
		temp2+=(width-2);
		temp[24]=objects[z+1][temp2];
		temp2++;
		temp[25]=objects[z+1][temp2];
		temp2++;
		temp[26]=objects[z+1][temp2];
		return temp;
	}

	// Does a 4-connected 3D flood fill
	// the data input stack has been flood filled already with the 2D algorithm
	// the object gets deleted from data
	public boolean fill3D(float[][] data,float[][] counter,int[] sliceblobs,int id,int x,int y,int z){
		stackSize=0;
		push(x,y,z);
		while(true){
			int[] cds=pop();
			x=cds[0];
			if(x<0)
				return true;
			y=cds[1];
			z=cds[2];
			if(data[z][x+y*width]==0.0f)
				continue;
			int[] lims2d=fb.getfilllimits(data[z],(int)data[z][x+y*width],x,y);
			//if(lims2d[0]<0) lims2d[0]=0;  if(lims2d[1]>=width) lims2d[1]=(width-1);
			//if(lims2d[2]<0) lims2d[2]=0;  if(lims2d[3]>=height) lims2d[3]=(height-1);
			fillRegion(data,counter,id,lims2d,z,(int)data[z][x+y*width]); // fill scan-region
			if(z>0){
				boolean[] prevselected=new boolean[sliceblobs[z-1]];
				for(int i=lims2d[2];i<=lims2d[3];i++){ // find unique scan-regions above this one
					for(int j=lims2d[0];j<=lims2d[1];j++){
						if(counter[z][j+i*width]==(float)id){
							int tempid=(int)data[z-1][j+i*width];
							if(tempid>0 && !prevselected[tempid-1]){
								prevselected[tempid-1]=true;
								push(j,i,z-1);
							}
						}
					}
				}
			}
			if(z<(depth-1)){
				boolean[] prevselected=new boolean[sliceblobs[z+1]];
				for(int i=lims2d[2];i<=lims2d[3];i++){ // find unique scan-regions below this one
					for(int j=lims2d[0];j<=lims2d[1];j++){
						if(counter[z][j+i*width]==(float)id){
							int tempid=(int)data[z+1][j+i*width];
							if(tempid>0&&!prevselected[tempid-1]){
								prevselected[tempid-1]=true;
								push(j,i,z+1);
							}
						}
					}
				}
			}
		}
	}

	final void fillRegion(float[][] data,float[][] counter,int id,int[] lims2d,int z,int oldid){
		for(int i=lims2d[2];i<=lims2d[3];i++){
			for(int j=lims2d[0];j<=lims2d[1];j++){
				if(data[z][j+i*width]==oldid){
					data[z][j+i*width]=0.0f;
					counter[z][j+i*width]=id;
				}
			}
		}
	}

	final void push(int x,int y,int z){
		stackSize++;
		if(stackSize==maxStackSize){
			int[][] newStack=new int[maxStackSize*2][];
			System.arraycopy(stack,0,newStack,0,maxStackSize);
			stack=newStack;
			maxStackSize*=2;
		}
		stack[stackSize-1]=new int[]{x,y,z};
	}

	final int[] pop(){
		if(stackSize<1)
			return new int[]{-1,-1,-1};
		else{
			int[] vals=stack[stackSize-1];
			stackSize--;
			return vals;
		}
	}
	
	public int[] get_areas(float[][] objects){
		int nblobs=get_nblobs(objects);
		int[] hist=new int[nblobs];
		for(int j=0;j<depth;j++){
    		for(int i=0;i<width*height;i++){
    			if(objects[j][i]>0.0f){
    				hist[(int)objects[j][i]-1]++;
    			}
    		}
		}
		return hist;
	}
	
	public void filter_area(float[][] objects,int[] arealims,boolean renumber){
		int[] hist=get_areas(objects);
		for(int i=0;i<hist.length;i++){
			if(hist[i]<arealims[0]||hist[i]>arealims[1]){
				delete_object(objects,(float)(i+1));
			}
		}
		if(renumber)
			renumber_objects(objects);
	}
	
	public Object[] get_object_outline(float[][] objects,int id){
		//this is incorrect--need to look for multiple polygons per slice
		Polygon[] outlines=new Polygon[objects.length];
		int[] pos=new int[objects.length];
		int counter=0;
		for(int i=0;i<objects.length;i++){
			Polygon temp=fb.get_object_outline(objects[i],id);
			if(temp!=null){
				outlines[counter]=temp;
				pos[counter]=i;
				counter++;
			}
		}
		Polygon[] outlines2=new Polygon[counter];
		System.arraycopy(outlines,0,outlines2,0,counter);
		int[] pos2=new int[counter];
		System.arraycopy(pos,0,pos2,0,counter);
		return new Object[]{outlines2,pos2};
	}

}
