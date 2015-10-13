/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.jsort;

public class findblobs{
	//this is a class that finds objects by a max not mask approach
	//objects are never allowed to be less than searchd distance apart 
	public int minarea,maxarea,searchd,maxblobs,width,height;
	public float thresh,minsep,edgebuf;
	public boolean usemaxpt;

	public findblobs(int width1,int height1,float[] criteria){
		width=width1;
		height=height1;
		usemaxpt=false;
		setcriteria(criteria);
	}

	public void setcriteria(float[] criteria){
		minarea=(int)criteria[0]; // min blob area
		maxarea=(int)criteria[1]; // max blob area
		searchd=(int)criteria[2]; // search diameter, not used
		maxblobs=(int)criteria[3]; // max number of blobs
		thresh=criteria[4]; // threshold
		minsep=criteria[5]; // min blob separation in pixels
		edgebuf=criteria[6]; // edge buffer in pixels
	}

	public float[][] dofindblobs2(float[] data,float[] mask){
		// here we find blobs a slightly different way
		// searchd is not used here, rather we use minsep
		boolean[] binmask=new boolean[width*height];
		int maxpt=0;
		int maxx=0;
		int maxy=0;
		float maxval=0.0f;
		boolean ptfound=false;
		boolean atedge=false;
		float[][] tempstats=new float[maxblobs][5];
		int blobcounter=0;
		int validblobs=0;
		int searchr=(int)(0.5f*minsep+0.5f);
		int searchr2=searchr*searchr;
		do{
			// start by finding the maximum point
			ptfound=false;
			maxpt=maxnotmask(data,binmask);
			maxval=data[maxpt];
			maxy=maxpt/width;
			maxx=maxpt-maxy*width;
			if(maxval>=thresh){
				ptfound=true;
				atedge=false;
				int upperx=(maxx+searchr);
				if(upperx>=width){
					upperx=width-1;
					atedge=true;
				}
				int lowerx=(maxx-searchr);
				if(lowerx<0){
					lowerx=0;
					atedge=true;
				}
				int uppery=(maxy+searchr);
				if(uppery>=height){
					uppery=height-1;
					atedge=true;
				}
				int lowery=(maxy-searchr);
				if(lowery<0){
					lowery=0;
					atedge=true;
				}
				if(maxx<edgebuf){
					atedge=true;
				}
				if(maxx>(width-1-edgebuf)){
					atedge=true;
				}
				if(maxy<edgebuf){
					atedge=true;
				}
				if(maxy>(height-1-edgebuf)){
					atedge=true;
				}
				if(!atedge){
					blobcounter++;
					float intval=0.0f;
					tempstats[blobcounter-1][2]=maxval;
					for(int i=lowery;i<=uppery;i++){
						for(int j=lowerx;j<=upperx;j++){
							if(((i-maxy)*(i-maxy)+(j-maxx)*(j-maxx))<=searchr2){
								int index=j+i*width;
								binmask[index]=true;
								if(!Float.isNaN(data[index])){
    								if(data[index]>=thresh){
    									mask[index]=blobcounter;
    									tempstats[blobcounter-1][0]+=data[index]*j;
    									tempstats[blobcounter-1][1]+=data[index]*i;
    									intval+=data[index];
    									tempstats[blobcounter-1][3]+=1.0f;
    								}
								}
							}
						}
					}
					tempstats[blobcounter-1][0]/=intval;
					tempstats[blobcounter-1][1]/=intval;
					if(usemaxpt){
						tempstats[blobcounter-1][0]=maxx;
						tempstats[blobcounter-1][1]=maxy;
					}
					if(tempstats[blobcounter-1][3]>=minarea&&tempstats[blobcounter-1][3]<=maxarea){
						tempstats[blobcounter-1][4]=1.0f;
						validblobs++;
					}
				}else{
					for(int i=lowery;i<=uppery;i++){
						for(int j=lowerx;j<=upperx;j++){
							if(((i-maxy)*(i-maxy)+(j-maxx)*(j-maxx))<=searchr2){
								binmask[j+i*width]=true;
							}
						}
					}
				}
			}
		}while(ptfound&&blobcounter<maxblobs);
		float[][] stats=new float[validblobs][4];
		int counter=0;
		for(int i=0;i<blobcounter;i++){
			if(tempstats[i][4]==1.0f){
				stats[counter][0]=tempstats[i][0];
				stats[counter][1]=tempstats[i][1];
				stats[counter][2]=tempstats[i][2];
				stats[counter][3]=tempstats[i][3];
				for(int j=0;j<width*height;j++){
					if(mask[j]==i+1){
						mask[j]=counter+1;
					}
				}
				counter++;
			}else{
				for(int j=0;j<width*height;j++){
					if(mask[j]==i+1){
						mask[j]=0.0f;
					}
				}
			}
		}
		// stats are centerx,centery,integral,area
		return stats;
	}
	
	public float[][] dofindblobs3(float[] data,float[] mask){
		//this is just like dofindblobs2 but we sort the intensities and keep the sort index
		//that speeds up the max search each time dramatically
		boolean[] binmask=new boolean[width*height];
		int[] order1=jsort.get_javasort_order(data);
		int[] order=new int[width*height];
		int temp=width*height;
		for(int i=0;i<temp;i++){
			order[i]=order1[temp-1-i];
		}
		order1=null;
		int maxpt=0;
		int maxx=0;
		int maxy=0;
		float maxval=0.0f;
		boolean ptfound=false;
		boolean atedge=false;
		float[][] tempstats=new float[maxblobs][5];
		int blobcounter=0;
		int validblobs=0;
		int searchr=(int)(0.5f*minsep+0.5f);
		int searchr2=searchr*searchr;
		int previndex=0;
		do{
			// start by finding the maximum point
			ptfound=false;
			int[] temp1=maxnotmask(data,binmask,order,previndex);
			if(temp1==null) break;
			maxpt=temp1[0];
			previndex=temp1[1];
			maxval=data[maxpt];
			maxy=maxpt/width;
			maxx=maxpt-maxy*width;
			if(maxval>=thresh){
				ptfound=true;
				atedge=false;
				int upperx=(maxx+searchr);
				if(upperx>=width){
					upperx=width-1;
					atedge=true;
				}
				int lowerx=(maxx-searchr);
				if(lowerx<0){
					lowerx=0;
					atedge=true;
				}
				int uppery=(maxy+searchr);
				if(uppery>=height){
					uppery=height-1;
					atedge=true;
				}
				int lowery=(maxy-searchr);
				if(lowery<0){
					lowery=0;
					atedge=true;
				}
				if(maxx<edgebuf){
					atedge=true;
				}
				if(maxx>(width-1-edgebuf)){
					atedge=true;
				}
				if(maxy<edgebuf){
					atedge=true;
				}
				if(maxy>(height-1-edgebuf)){
					atedge=true;
				}
				if(!atedge){
					blobcounter++;
					float intval=0.0f;
					tempstats[blobcounter-1][2]=maxval;
					for(int i=lowery;i<=uppery;i++){
						for(int j=lowerx;j<=upperx;j++){
							if(((i-maxy)*(i-maxy)+(j-maxx)*(j-maxx))<=searchr2){
								int index=j+i*width;
								binmask[index]=true;
								if(!Float.isNaN(data[index])){
    								if(data[index]>=thresh){
    									mask[index]=blobcounter;
    									tempstats[blobcounter-1][0]+=data[index]*j;
    									tempstats[blobcounter-1][1]+=data[index]*i;
    									intval+=data[index];
    									tempstats[blobcounter-1][3]+=1.0f;
    								}
								}
							}
						}
					}
					tempstats[blobcounter-1][0]/=intval;
					tempstats[blobcounter-1][1]/=intval;
					if(usemaxpt){
						tempstats[blobcounter-1][0]=maxx;
						tempstats[blobcounter-1][1]=maxy;
					}
					if(tempstats[blobcounter-1][3]>=minarea&&tempstats[blobcounter-1][3]<=maxarea){
						tempstats[blobcounter-1][4]=1.0f;
						validblobs++;
					}
				}else{
					for(int i=lowery;i<=uppery;i++){
						for(int j=lowerx;j<=upperx;j++){
							if(((i-maxy)*(i-maxy)+(j-maxx)*(j-maxx))<=searchr2){
								binmask[j+i*width]=true;
							}
						}
					}
				}
			}
		}while(ptfound&&blobcounter<maxblobs);
		float[][] stats=new float[validblobs][4];
		int counter=0;
		for(int i=0;i<blobcounter;i++){
			if(tempstats[i][4]==1.0f){
				stats[counter][0]=tempstats[i][0];
				stats[counter][1]=tempstats[i][1];
				stats[counter][2]=tempstats[i][2];
				stats[counter][3]=tempstats[i][3];
				for(int j=0;j<width*height;j++){
					if(mask[j]==i+1){
						mask[j]=counter+1;
					}
				}
				counter++;
			}else{
				for(int j=0;j<width*height;j++){
					if(mask[j]==i+1){
						mask[j]=0.0f;
					}
				}
			}
		}
		// stats are centerx,centery,integral,area
		return stats;
	}
	
	public float[] get_indexed_objects(float[][] coords,boolean renumber){
		//here we use the coordinates to create a masked image with indexed objects like one would get from thresholding
		//this works best if usemaxpt is selected
		//start with the last coordinate because it will be masked by others
		float[] objects=new float[width*height];
		for(int i=0;i<coords.length;i++){
			int id=coords.length-i;
			set_circle_val(id,objects,coords[coords.length-1-i][0],coords[coords.length-1-i][1],minsep,width,height);
		}
		//now set any pixels which are neighbored by different values to zero
		findblobs3 fb=new findblobs3(width,height);
		fb.separateobjects(objects);
		//and renumber in standard raster discovery order if asked for
		if(renumber) fb.renumber_objects(objects);
		return objects;
	}
	
	public static void set_circle_val(float val,float[] image,float x,float y,float diameter,int width,int height){
		float rad=0.5f*diameter;
		float rad2=rad*rad;
		int xstart=(int)(x-rad); if(xstart<0) xstart=0;
		int xend=1+(int)(x+rad); if(xend>=width) xend=width-1;
		int ystart=(int)(y-rad); if(ystart<0) ystart=0;
		int yend=1+(int)(y+rad); if(yend>=height) yend=height-1;
		for(int i=ystart;i<=yend;i++){
			for(int j=xstart;j<=xend;j++){
				float d2=(j-x)*(j-x)+(i-y)*(i-y);
				if(d2<=rad2) image[j+i*width]=val;
			}
		}
	}

	private int maxnotmask(float[] data,boolean[] binmask){
		float max=-1000000.0f;
		int imax=0;
		for(int i=0;i<data.length;i++){
			if(Float.isNaN(data[i])) continue;
			if(data[i]>max&&!binmask[i]){
				max=data[i];
				imax=i;
			}
		}
		return imax;
	}
	
	private int[] maxnotmask(float[] data,boolean[] binmask,int[] order,int start){
		for(int i=start;i<data.length;i++){
			int index=order[i];
			if(Float.isNaN(data[index])) continue;
			if(!binmask[index]){
				return new int[]{index,i};
			}
		}
		return null;
	}
}
