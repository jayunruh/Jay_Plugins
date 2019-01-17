/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.algutils;
import jalgs.gui_interface;
import jalgs.jstatistics;

import java.awt.Polygon;
import java.util.ArrayList;
import java.util.List;

import quickhull3d.Point3d;
import quickhull3d.QuickHull3D;

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
	public gui_interface gui;

	public findblobs3D(int width1,int height1,int depth1){
		width=width1;
		height=height1;
		depth=depth1;
		fb=new findblobs3(width,height);
		nobjects=0;
		gui=null;
	}
	
	public findblobs3D(int width1,int height1,int depth1,gui_interface gui){
		width=width1;
		height=height1;
		depth=depth1;
		fb=new findblobs3(width,height);
		nobjects=0;
		this.gui=gui;
	}
	
	public void set_objects(float[][] objects){
		nobjects=(int)maxarray(objects);
	}
	
	public float[][] dofindblobs(Object[] data1){
		if(data1[0] instanceof byte[]){
			byte[][] temp=new byte[data1.length][];
			for(int i=0;i<data1.length;i++) temp[i]=(byte[])data1[i];
			return dofindblobs(temp);
		} else {
			float[][] temp=new float[data1.length][];
			for(int i=0;i<temp.length;i++) temp[i]=(float[])data1[i];
			return dofindblobs(temp,0.5f);
		}
	}

	public float[][] dofindblobs(byte[][] data1){
		float[][] temp=new float[depth][];
		int[] sliceblobs=new int[depth];
		int[][][] filllims=new int[depth][][];
		for(int i=0;i<depth;i++){
			temp[i]=fb.dofindblobs(data1[i]);
			sliceblobs[i]=fb.nobjects;
			filllims[i]=fb.getallfilllimits(temp[i]);
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		// objects are removed from the temp array as we fill them
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,filllims,id,k,j,i);
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		nobjects=id;
		return objects;
	}

	public float[][] dofindblobs(float[][] data1,float thresh){
		float[][] temp=new float[data1.length][];
		int[] sliceblobs=new int[depth];
		int[][][] filllims=new int[depth][][];
		for(int i=0;i<data1.length;i++){
			temp[i]=fb.dofindblobs(data1[i],thresh);
			sliceblobs[i]=fb.nobjects;
			filllims[i]=fb.getallfilllimits(temp[i]);
		}
		float[][] objects=new float[depth][width*height];
		// now go through and assemble the 2D objects into 3D objects
		// objects are removed from the temp array as we fill them
		int id=0;
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					if(temp[i][k+j*width]>0.0f){
						id++;
						fill3D(temp,objects,sliceblobs,filllims,id,k,j,i);
					}
				}
			}
		}
		nobjects=id;
		return objects;
	}
	
	public byte[][] tobinary(float[][] objects){
		return tobinary(objects,false);
	}

	public byte[][] tobinary(float[][] objects,boolean separate){
		byte[][] out=new byte[depth][width*height];
		if(separate) separateobjects(objects);
		for(int i=0;i<depth;i++){
			for(int j=0;j<width*height;j++)
				if(objects[i][j]>0.0f)
					out[i][j]=(byte)255;
		}
		return out;
	}

	public void clear_edges(float[][] objects,boolean renumber){
		int totlength=width*height+width*depth+height*depth;
		int counter=0;
		for(int i=0;i<width*height;i++){
			// get the top and bottom
			if(objects[0][i]>0.0f){
				delete_object(objects,objects[0][i]);
			}
			if(objects[depth-1][i]>0.0f){
				delete_object(objects,objects[depth-1][i]);
			}
			counter++;
			if(gui!=null) gui.showProgress(counter,totlength);
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
				counter++;
				if(gui!=null) gui.showProgress(counter,totlength);
			}
			for(int i=0;i<height;i++){
				if(objects[k][i*width]>0.0f){
					delete_object(objects,objects[k][i*width]);
				}
				if(objects[k][i*width+width-1]>0.0f){
					delete_object(objects,objects[k][i*width+width-1]);
				}
				counter++;
				if(gui!=null) gui.showProgress(counter,totlength);
			}
		}
		if(renumber) renumber_objects(objects);
	}
	
	public void separateobjects(float[][] objects){
		//here non-zero pixels neighboring another object are turned black
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					float val=cloned[i][k+j*width];
					if(val>0.0f){
						float[] neighbors=getNeighbors(cloned,k,j,i);
						for(int l=0;l<neighbors.length;l++){
							if(neighbors[l]>0.0f && neighbors[l]!=val){
								objects[i][j*width+k]=0.0f;
								break;
							}
						}
					}
				}
			}
		}
		/*byte[][] temp=tobinary(objects,false);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);*/
	}
	
	public void erodeobjects(float[][] objects,int threshold){
		//pixels neighbored in 3D by space or another object are deleted
		//objects are allowed to separate (become multiple new objects)
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
	
	public void erodeobjects2(float[][] objects,int threshold){
		//pixels neighbored in 3D by space or another object are deleted
		//objects are not allowed to separate
		//they are not renumbered either
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
		//byte[][] temp=tobinary(objects,false);
		//float[][] temp2=dofindblobs(temp);
		//copy_objects(temp2,objects);
	}
	
	public void dilateobjects(float[][] objects){
		//pixels neighbored in 3D by one object are added to that object
		//objects are not allowed to merge;
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					if(cloned[i][k+j*width]==0.0f){
    					float[] neighbors=getNeighbors(cloned,k,j,i);
    					float val=0.0f; boolean found=false; boolean found2=false;
    					for(int l=0;l<neighbors.length;l++){
    						if(neighbors[l]>0.0f){
    							if(!found){
    								found=true;
    								val=neighbors[l];
    							} else {
    								if(neighbors[l]!=val){
    									found2=true; break;
    								}
    							}
    						}
    					}
    					if(found && !found2) objects[i][j*width+k]=val;
					}
				}
			}
		}
		byte[][] temp=tobinary(objects);
		float[][] temp2=dofindblobs(temp);
		copy_objects(temp2,objects);
	}
	
	public void dilateobjects2(float[][] objects){
		//pixels neighbored in 3D by one object are added to that object
		//here objects are allowed to merge
		float[][] cloned=algutils.clone_multidim_array(objects);
		for(int i=1;i<(depth-1);i++){
			for(int j=1;j<(height-1);j++){
				for(int k=1;k<(width-1);k++){
					if(cloned[i][k+j*width]==0.0f){
    					float[] neighbors=getNeighbors(cloned,k,j,i);
    					for(int l=0;l<neighbors.length;l++){
    						if(neighbors[l]>0.0f){
    							objects[i][j*width+k]=1.0f;
    							break;
    						}
    					}
					}
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
	
	public void delete_object(float[][] objects,float id,int[] filllims){
		nobjects--;
		for(int k=filllims[4];k<=filllims[5];k++){
			for(int i=filllims[2];i<=filllims[3];i++){
				for(int j=filllims[0];j<=filllims[1];j++){
					if(objects[k][j+i*width]==id){
						objects[k][j+i*width]=0.0f;
					}
				}
			}
		}
	}
	
	public void delete_objects(float[][] objects,int[] ids){
		nobjects-=ids.length;
		for(int k=0;k<depth;k++){
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					int id=(int)objects[k][j+i*width];
					if(id>0){
						for(int l=0;l<ids.length;l++){
							if(id==ids[l]){objects[k][j+i*width]=0.0f; break;}
						}
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
				newids[i]=counter;
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
	public boolean fill3D(float[][] data,float[][] counter,int[] sliceblobs,int[][][] allfilllims,int id,int x,int y,int z){
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
			//int[] lims2d=fb.getfilllimits(data[z],(int)data[z][x+y*width],x,y);
			int[] lims2d=allfilllims[z][(int)data[z][x+y*width]-1];
			//if(lims2d[0]<0) lims2d[0]=0;  if(lims2d[1]>=width) lims2d[1]=(width-1);
			//if(lims2d[2]<0) lims2d[2]=0;  if(lims2d[3]>=height) lims2d[3]=(height-1);
			fillRegion(data,counter,id,lims2d,z,(int)data[z][x+y*width]); // fill scan-region
			if(z>0){
				boolean[] prevselected=new boolean[sliceblobs[z-1]];
				for(int i=lims2d[2];i<=lims2d[3];i++){ // find unique scan-regions above this one
					for(int j=lims2d[0];j<=lims2d[1];j++){
						if(counter[z][j+i*width]==id){
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
						if(counter[z][j+i*width]==id){
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
	
	public int[][] getallfilllimits(float[][] objects){
		//int nobjects=(int)maxarray(objects);
		//the limits are lowerx,upperx,lowery,uppery,lowerz,upperz
		int[][] lims=new int[nobjects][6];
		for(int i=0;i<nobjects;i++){
			lims[i][0]=width-1; lims[i][1]=0;
			lims[i][2]=height-1; lims[i][3]=0;
			lims[i][4]=depth-1; lims[i][5]=0;
		}
		for(int k=0;k<depth;k++){
    		for(int i=0;i<height;i++){
    			for(int j=0;j<width;j++){
    				if(objects[k][j+i*width]>0.0f){
    					int id=(int)objects[k][j+i*width]-1;
    					if(j<lims[id][0]) lims[id][0]=j;
    					if(j>lims[id][1]) lims[id][1]=j;
    					if(i<lims[id][2]) lims[id][2]=i;
    					if(i>lims[id][3]) lims[id][3]=i;
    					if(k<lims[id][4]) lims[id][4]=k;
    					if(k>lims[id][5]) lims[id][5]=k;
    				}
    			}
    		}
		}
		for(int i=0;i<nobjects;i++){
			if(lims[i][1]<lims[i][0]){
				//here we never found an object, so reset full limits
				lims[i][0]=0; lims[i][1]=width-1;
				lims[i][2]=0; lims[i][3]=height-1;
				lims[i][4]=depth-1; lims[i][5]=0;
			}
		}
		return lims;
	}
	
	public Object[] projectObjects(float[][] objects,int method,boolean keep3D){
		//here we find the z projection for each object (0 is largest, 1 is sum) and place it in a 3D space
		//if keep3D is false, we place it in the first slice
		//need to get a z area profile for each object
		int[][] zprof=new int[nobjects][depth];
		for(int i=0;i<depth;i++) {
			for(int j=0;j<height;j++) {
				for(int k=0;k<width;k++) {
					int id=(int)objects[i][k+j*width];
					if(id>0 && id<=nobjects) {
						zprof[id-1][i]++;
					}
				}
			}
		}
		//now find either the average z position or the maximum z position
		float[] zposstats=new float[nobjects];
		float[] areas=new float[nobjects];
		for(int i=0;i<nobjects;i++) {
			areas[i]=jstatistics.getstatistic("Sum",zprof[i],null);
			if(method==0) zposstats[i]=jstatistics.getstatistic("MaxPos",zprof[i],null);
			if(method==1) zposstats[i]=jstatistics.getstatistic("AvgPos",zprof[i],null);
		}
		float[][] newobj=null;
		if(keep3D) newobj=new float[depth][width*height];
		else newobj=new float[1][width*height];
		int[][] filllims=getallfilllimits(objects);
		for(int i=0;i<nobjects;i++) {
			if(areas[i]>0) {
				if(method==0) {
					//copy the largest slice to the new image
					int bestslice=(int)zposstats[i];
					if(bestslice<0) bestslice=0; if(bestslice>=depth) bestslice=(depth-1);
					for(int j=filllims[i][2];j<=filllims[i][3];j++) {
						for(int k=filllims[i][0];k<=filllims[i][1];k++) {
							int id=(int)objects[bestslice][k+j*width];
							if(id==(i+1)) {
								if(keep3D) newobj[bestslice][k+j*width]=(float)id;
								else newobj[0][k+j*width]=(float)id;
							}
						}
					}
				} else {
					//project the object onto its best slice
					int bestslice=Math.round(zposstats[i]);
					if(bestslice<0) bestslice=0; if(bestslice>=depth) bestslice=(depth-1);
					for(int l=filllims[i][4];l<=filllims[i][5];l++) {
    					for(int j=filllims[i][2];j<=filllims[i][3];j++) {
    						for(int k=filllims[i][0];k<=filllims[i][1];k++) {
    							int id=(int)objects[l][k+j*width];
    							if(id==(i+1)) {
    								if(keep3D) newobj[bestslice][k+j*width]=(float)id;
    								else newobj[0][k+j*width]=(float)id;
    							}
    						}
    					}
					}
				}
			}
		}
		return new Object[] {newobj,zposstats};
	}
	
	/********************************************************
	 * gets all of of the object stats
	 * @param objects
	 * @param measurement: the intensity image to measure
	 * @param lims
	 * @param stat
	 * @return
	 */
	public float[] get_all_object_stats(float[][] objects,Object[] measurement,int[][] lims,String stat){
		float[] stats=new float[nobjects];
		int[] areas=get_areas(objects);
		for(int i=0;i<nobjects;i++){
			float[] temp=new float[areas[i]];
			int counter=0;
			if(measurement[0] instanceof float[]){
    			for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((float[])measurement[j])[xyindex];
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0] instanceof short[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((short[])measurement[j])[xyindex]&0xffff;
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0] instanceof byte[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							temp[counter]=((byte[])measurement[j])[xyindex]&0xff;
    							counter++;
    						}
    					}
    				}
    			}
			}
			stats[i]=jstatistics.getstatistic(stat,temp,null);
			if(gui!=null) gui.showProgress(i,nobjects);
		}
		return stats;
	}
	
	/********************************************************
	 * gets all of of the object stats from multiple channels
	 * @param objects
	 * @param measurement: the intensity image to measure
	 * @param lims
	 * @param stat
	 * @return
	 */
	public float[][] get_all_object_stats(float[][] objects,Object[][] measurement,int[][] lims,String stat){
		float[][] stats=new float[measurement.length][nobjects];
		int[] areas=get_areas(objects);
		for(int i=0;i<nobjects;i++){
			float[][] temp=new float[measurement.length][areas[i]];
			int counter=0;
			if(measurement[0][0] instanceof float[]){
    			for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							for(int m=0;m<measurement.length;m++) temp[m][counter]=((float[])measurement[m][j])[xyindex];
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0][0] instanceof short[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							for(int m=0;m<measurement.length;m++) temp[m][counter]=((short[])measurement[m][j])[xyindex]&0xffff;
    							counter++;
    						}
    					}
    				}
    			}
			} else if(measurement[0][0] instanceof byte[]){
				for(int j=lims[i][4];j<=lims[i][5];j++){
    				for(int k=lims[i][2];k<=lims[i][3];k++){
    					for(int l=lims[i][0];l<=lims[i][1];l++){
    						int xyindex=l+k*width;
    						if((int)objects[j][xyindex]==(i+1)){
    							for(int m=0;m<measurement.length;m++) temp[m][counter]=((byte[])measurement[m][j])[xyindex]&0xff;
    							counter++;
    						}
    					}
    				}
    			}
			}
			for(int j=0;j<temp.length;j++) stats[j][i]=jstatistics.getstatistic(stat,temp[j],null);
			if(gui!=null) gui.showProgress(i,nobjects);
		}
		return stats;
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
	
	public int[] get_areas(float[][] objects,float[][] lims){
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
		//this version uses fill limits
		int[] hist=get_areas(objects);
		int[][] filllims=getallfilllimits(objects);
		for(int i=0;i<hist.length;i++){
			if(hist[i]<arealims[0]||hist[i]>arealims[1]){
				delete_object(objects,i+1,filllims[i]);
			}
			if(gui!=null) gui.showProgress(i,hist.length);
		}
		if(renumber)
			renumber_objects(objects);
	}
	
	public void filter_area_surface(float[][] objects,int[] arealims,int threshcount,boolean renumber){
		float[][] stats=getCentroidsAreasSurface(objects,threshcount);
		int[][] filllims=getallfilllimits(objects);
		for(int i=0;i<stats.length;i++){
			if((int)stats[3][i]<arealims[0]||(int)stats[3][i]>arealims[1] || (int)stats[4][i]<arealims[2]||(int)stats[4][i]>arealims[3]){
				delete_object(objects,i+1,filllims[i]);
			}
			if(gui!=null) gui.showProgress(i,stats.length);
		}
		if(renumber)
			renumber_objects(objects);
	}
	
	public float[][] getCentroidsAreas(float[][] objects){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][4];
		for(int i=0;i<objects.length;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
						centroids[index-1][0]+=k;
						centroids[index-1][1]+=j;
						centroids[index-1][2]+=i;
						centroids[index-1][3]+=1.0f;
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
		}
		return centroids;
	}
	
	public float[][] getCentroidsAreasAvgs(float[][] objects,Object[] measurement){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][5];
		for(int i=0;i<objects.length;i++){
			float[] tempmeas=algutils.convert_arr_float2(measurement[i]);
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
						centroids[index-1][0]+=k;
						centroids[index-1][1]+=j;
						centroids[index-1][2]+=i;
						centroids[index-1][3]+=1.0f;
						centroids[index-1][4]+=tempmeas[k+j*width];
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
			centroids[i][4]/=centroids[i][3];
		}
		return centroids;
	}
	
	/********************
	 * here we have multiple measurement channels
	 * @param objects
	 * @param measurement
	 * @return
	 */
	public float[][] getCentroidsAreasAvgs(float[][] objects,Object[][] measurement){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][4+measurement.length];
		for(int i=0;i<objects.length;i++){ //loop over z
			float[][] tempmeas=new float[measurement.length][];
			for(int j=0;j<measurement.length;j++) tempmeas[j]=algutils.convert_arr_float2(measurement[j][i]);
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
						centroids[index-1][0]+=k;
						centroids[index-1][1]+=j;
						centroids[index-1][2]+=i;
						centroids[index-1][3]+=1.0f;
						for(int l=0;l<measurement.length;l++) centroids[index-1][4+l]+=tempmeas[l][k+j*width];
					}
				}
				if(gui!=null) gui.showProgress(j+i*height,depth*height);
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
			for(int l=0;l<measurement.length;l++) centroids[i][4+l]/=centroids[i][3];
		}
		return centroids;
	}
	
	public float[][] getCentroidsAreasSurface(float[][] objects,int threshcount){
		int nblobs=get_nblobs(objects);
		float[][] centroids=new float[nblobs][5];
		//assume that edge voxels are surface
		for(int i=0;i<depth;i++){
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					int index=(int)objects[i][k+j*width];
					if(index>0){
    					float[] neighbors=getNeighbors(objects,k,j,i);
    					if(neighbors==null) centroids[index-1][4]+=1.0f;
    					else {
        					int count=0;
        					for(int l=0;l<neighbors.length;l++) if((int)neighbors[l]!=index) count++;
        					if(count>=threshcount) centroids[index-1][4]+=1.0f;
    					}
    					centroids[index-1][3]+=1.0f;
    					centroids[index-1][0]+=k;
    					centroids[index-1][1]+=j;
    					centroids[index-1][2]+=i;
					}
				}
			}
		}
		for(int i=0;i<nblobs;i++){
			centroids[i][0]/=centroids[i][3];
			centroids[i][1]/=centroids[i][3];
			centroids[i][2]/=centroids[i][3];
		}
		return centroids;
	}
	
	public float[] getLongestDimensions(float[][] objects){
		int nblobs=get_nblobs(objects);
		int[][] bounds=getallfilllimits(objects);
		float[] outarr=new float[nblobs];
		for(int i=0;i<nblobs;i++){
			outarr[i]=getLongestDimension(objects,bounds[i],i+1);
		}
		return outarr;
	}
	
	/**************
	 * gets and arraylist of all of the coordinates in an object
	 * @param objects
	 * @param bounds
	 * @param id
	 * @return
	 */
	public List<float[]> getObjectCoords(float[][] objects,int[] bounds,int id){
		List<float[]> coords=new ArrayList<float[]>();
		for(int i=bounds[4];i<=bounds[5];i++){
			for(int j=bounds[2];j<=bounds[3];j++){
				for(int k=bounds[0];k<=bounds[1];k++){
					if((int)objects[i][k+j*width]==id){
						coords.add(new float[]{k,j,i});
					}
				}
			}
		}
		return coords;
	}
	
	public List<float[]> getObjectSurfCoords(float[][] objects,int[] bounds,int id,int edgethresh){
		List<float[]> coords=new ArrayList<float[]>();
		for(int i=bounds[4];i<=bounds[5];i++){
			for(int j=bounds[2];j<=bounds[3];j++){
				for(int k=bounds[0];k<=bounds[1];k++){
					if((int)objects[i][k+j*width]==id){
						float[] neighbors=getNeighbors(objects,k,j,i);
						if(neighbors==null) coords.add(new float[]{k,j,i}); //on the image edge
						else{
    						int temp=0;
    						for(int l=0;l<neighbors.length;l++) if(neighbors[l]!=(float)id) temp++;
    						if(temp>=edgethresh) coords.add(new float[]{k,j,i});
						}
					}
				}
			}
		}
		return coords;
	}
	
	public List<float[]> getConvexHull(float[][] objects,int[] bounds,int id){
		List<float[]> surf=getObjectSurfCoords(objects,bounds,id,1);
		Point3d[] points=new Point3d[surf.size()];
		for(int i=0;i<surf.size();i++) points[i]=new Point3d(surf.get(i)[0],surf.get(i)[1],surf.get(i)[2]);
		QuickHull3D hull=new QuickHull3D(points);
		Point3d[] hullpoints=hull.getVertices();
		List<float[]> hullpoints2=new ArrayList<float[]>();
		for(int i=0;i<hullpoints.length;i++){
			hullpoints2.add(new float[]{(float)hullpoints[i].x,(float)hullpoints[i].y,(float)hullpoints[i].z});
		}
		return hullpoints2;
	}
	
	/****************
	 * returns the max distance between two object points
	 * @param objects
	 * @param bounds
	 * @param id
	 * @return
	 */
	public float getLongestDimension(float[][] objects,int[] bounds,int id){
		List<float[]> coords=getObjectCoords(objects,bounds,id);
		float maxdist2=0.0f;
		for(int i=0;i<coords.size();i++){
			float[] query=coords.get(i);
			for(int j=i+1;j<coords.size();j++){
				float[] bait=coords.get(j);
				float tempdist2=(query[0]-bait[0])*(query[0]-bait[0])+(query[1]-bait[1])*(query[1]-bait[1])+(query[2]-bait[2])*(query[2]-bait[2]);
				if(tempdist2>maxdist2) maxdist2=tempdist2;
			}
		}
		return (float)Math.sqrt(maxdist2);
	}
	
	public Object[] get_object_outline(float[][] objects,int id){
		//this is incorrect--need to look for multiple polygons per slice
		List<Polygon> polylist=new ArrayList<Polygon>();
		List<Integer> poslist=new ArrayList<Integer>();
		for(int i=0;i<objects.length;i++){
			byte[] temp=new byte[objects[i].length];
			boolean found=false;
			for(int j=0;j<objects[i].length;j++) if(objects[i][j]==(float)id){temp[j]=(byte)255; found=true;}
			if(found){
				float[] idobj=fb.dofindblobs(temp);
				Polygon[] polys=fb.get_object_outlines(idobj);
				if(polys==null) continue;
				for(int j=0;j<polys.length;j++){
					polylist.add(polys[j]);
					poslist.add(i);
				}
			}
		}
		Polygon[] outlines=new Polygon[polylist.size()];
		int[] pos=new int[polylist.size()];
		for(int i=0;i<outlines.length;i++){
			outlines[i]=polylist.get(i);
			pos[i]=(int)poslist.get(i);
		}
		return new Object[]{outlines,pos};
	}

}