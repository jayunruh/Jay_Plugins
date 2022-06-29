package jalgs.jseg;

import java.util.ArrayList;
import java.util.List;

import jalgs.algutils;

public class FindBlobsFast{
	//this class attempts to subdivide an image into multiple parts and find each set of blobs independently
	//then combine the sets together by reconciling the edges
	//this should save significant time for large image analysis
	public int minarea,maxarea,searchd,maxblobs,width,height,subsize;
	public float thresh,minsep,edgebuf,zedgebuf;
	public boolean usemaxpt;
	public findblobs fb;
	
	public FindBlobsFast(int width,int height,int subsize,float[] criteria) {
		this.width=width;
		this.height=height;
		this.subsize=subsize;
		fb=new findblobs(subsize,subsize,criteria);
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
		float[] tempcriteria=criteria.clone();
		//tempcriteria[6]=1+(int)(0.5f*(float)minsep);
		tempcriteria[6]=0;
		fb.setcriteria(tempcriteria);
	}
	
	public float[][] dofindblobs3(float[] data) {
		//the mask contains the indexing on output
		//work our way through the image creating subsize x subsize thumbnails, searching them, and then merging them
		//will enforce the edge buffer retroactively after searching
		int overlap=1+(int)minsep;
		int nwide=(int)((float)width/(float)(subsize-overlap)); //this is the number of full images (minus overlap) in x
		int xextra=width-nwide*(subsize-overlap);
		int nhigh=(int)((float)height/(float)(subsize-overlap)); //this is the number of full images (minus overlap) in y
		int yextra=height-nhigh*(subsize-overlap);
		if(thresh!=fb.thresh) fb.thresh=thresh;
		List<float[]> points=new ArrayList<float[]>();
		for(int i=0;i<nhigh;i++) {
			int ypos=i*(subsize-overlap);
			for(int j=0;j<nwide;j++) {
				int xpos=j*(subsize-overlap);
				//get the subregion
				//float[] subdata=algutils.get_region2(data,xpos,ypos,subsize,subsize,width,height);
				float[] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//find the maxima in that subregion
				float[][] temppts=fb.dofindblobs3(subdata,new float[subsize*subsize]);
				//now corect those values based on our position
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][4]+=(float)xpos;
					temppts[k][5]+=(float)ypos;
				}
				//merge the points into the master list
				mergeLists(points,temppts);
			}
			//at the end of each row, need to add the extra partial image, pad it to keep the size the same
			if(xextra>minsep) {
				int xpos=nwide*(subsize-overlap);
				//float[] subdata=algutils.get_region2(data,xpos,ypos,xextra,subsize,width,height);
				//subdata=padImage(subdata,xextra,subsize,subsize,subsize);
				float[] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				float[][] temppts=fb.dofindblobs3(subdata,new float[subsize*subsize]);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][4]+=(float)xpos;
					temppts[k][5]+=(float)ypos;
				}
				mergeLists(points,temppts);
			}
		}
		//in the last row, we need to pad all of the partial images
		if(yextra>minsep) {
			int ypos=nhigh*(subsize-overlap);
			for(int j=0;j<nwide;j++) {
				int xpos=j*(subsize-overlap);
				//float[] subdata=algutils.get_region2(data,xpos,ypos,subsize,yextra,width,height);
				//subdata=padImage(subdata,subsize,yextra,subsize,subsize);
				float[] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				float[][] temppts=fb.dofindblobs3(subdata,new float[subsize*subsize]);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][4]+=(float)xpos;
					temppts[k][5]+=(float)ypos;
				}
				mergeLists(points,temppts);
			}
			//finally the corner right image
			if(xextra>minsep) {
				int xpos=nwide*(subsize-overlap);
				float[] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//float[] subdata=algutils.get_region2(data,xpos,ypos,xextra,yextra,width,height);
				//subdata=padImage(subdata,xextra,yextra,subsize,subsize);
				float[][] temppts=fb.dofindblobs3(subdata,new float[subsize*subsize]);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][4]+=(float)xpos;
					temppts[k][5]+=(float)ypos;
				}
				mergeLists(points,temppts);
			}
		}
		//finally we enforce the edge buffer and output the points
		return enforceEdges(points);
	}
	
	public void mergeLists(List<float[]> master,float[][] source) {
		//here we add the points in the source list to the master list
		//overlaps must be reconciled by max intensity
		// stats are centerx,centery,integral,area,maxx,maxy,maxval
		//note that the centerx, centery values might be center of mass
		//assume that the sourcelist is sorted from high intensity to low
		if(master.size()==0) {
			//if the master list is empty
			for(int i=0;i<source.length;i++) master.add(source[i]);
			return;
		}
		//float mindist2=minsep*minsep; 
		int searchr=(int)(0.5f*minsep+0.5f);
		int searchr2=searchr*searchr; //calculate the squared value to avoid square roots later
		int oldsize=master.size();
		for(int i=0;i<source.length;i++) {
			float x=source[i][4];
			float y=source[i][5];
			//for each point start by finding the closest master point
			float[] temp0=master.get(0);
			float closestdist2=(x-temp0[4])*(x-temp0[4])+(y-temp0[5])*(y-temp0[5]);
			int closestpoint=0;
			for(int j=1;j<oldsize;j++) {
				float[] temp=master.get(j);
				float dist2=(x-temp[4])*(x-temp[4])+(y-temp[5])*(y-temp[5]);
				if(dist2<=closestdist2) {
					closestdist2=dist2;
					closestpoint=j;
				}
			}
			//now that we found the closest point, see if it is too close
			if(closestdist2<searchr2) {
				//this point is too close, keep the one with higher intensity
				if(source[i][6]>master.get(closestpoint)[6]) {
					master.set(closestpoint,source[i]);
				}
			} else {
				//this point doesn't conflict with the original list, add it to the list at the end
				master.add(source[i]);
			}
		}
		return;
	}
	
	public float[][] dofindblobs3D(float[][] data,float searchrz,float zedgebuf) {
		//work our way through the image creating subsize x subsize thumbnails, searching them, and then merging them
		//will enforce the edge buffer retroactively after searching
		//for now, do not use subregions in z
		int searchr=(int)(0.5f*minsep+0.5f);
		float zratio=searchrz/(float)searchr;
		int overlap=1+(int)minsep;
		//int zoverlap=1+(int)searchrz;
		int nwide=(int)((float)width/(float)(subsize-overlap)); //this is the number of full images (minus overlap) in x
		int xextra=width-nwide*(subsize-overlap);
		int nhigh=(int)((float)height/(float)(subsize-overlap)); //this is the number of full images (minus overlap) in y
		int yextra=height-nhigh*(subsize-overlap);
		if(thresh!=fb.thresh) fb.thresh=thresh;
		List<float[]> points=new ArrayList<float[]>();
		for(int i=0;i<nhigh;i++) {
			int ypos=i*(subsize-overlap);
			for(int j=0;j<nwide;j++) {
				int xpos=j*(subsize-overlap);
				//get the subregion
				float[][] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//float[][] subdata=algutils.get_region2(data,xpos,ypos,subsize,subsize,width,height);
				//find the maxima in that subregion
				float[][] temppts=fb.dofindblobs3D(subdata,null,searchrz,0);
				//now corect those values based on our position
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][5]+=(float)xpos;
					temppts[k][6]+=(float)ypos;
				}
				//merge the points into the master list
				mergeLists3D(points,temppts,zratio);
			}
			//at the end of each row, need to add the extra partial image, pad it to keep the size the same
			if(xextra>minsep) {
				int xpos=nwide*(subsize-overlap);
				float[][] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//float[][] subdata=algutils.get_region2(data,xpos,ypos,xextra,subsize,width,height);
				//subdata=padImage(subdata,xextra,subsize,subsize,subsize);
				float[][] temppts=fb.dofindblobs3D(subdata,null,searchrz,0);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][5]+=(float)xpos;
					temppts[k][6]+=(float)ypos;
				}
				mergeLists3D(points,temppts,zratio);
			}
		}
		//in the last row, we need to pad all of the partial images
		if(yextra>minsep) {
			int ypos=nhigh*(subsize-overlap);
			for(int j=0;j<nwide;j++) {
				int xpos=j*(subsize-overlap);
				float[][] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//float[][] subdata=algutils.get_region2(data,xpos,ypos,subsize,yextra,width,height);
				//subdata=padImage(subdata,subsize,yextra,subsize,subsize);
				float[][] temppts=fb.dofindblobs3D(subdata,null,searchrz,0);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][5]+=(float)xpos;
					temppts[k][6]+=(float)ypos;
				}
				mergeLists3D(points,temppts,zratio);
			}
			//finally the corner right image
			if(xextra>minsep) {
				int xpos=nwide*(subsize-overlap);
				float[][] subdata=getRegionPad(data,xpos,ypos,subsize,subsize,width,height);
				//float[][] subdata=algutils.get_region2(data,xpos,ypos,xextra,yextra,width,height);
				//subdata=padImage(subdata,xextra,yextra,subsize,subsize);
				float[][] temppts=fb.dofindblobs3D(subdata,null,searchrz,0);
				for(int k=0;k<temppts.length;k++) {
					temppts[k][0]+=(float)xpos;
					temppts[k][1]+=(float)ypos;
					temppts[k][5]+=(float)xpos;
					temppts[k][6]+=(float)ypos;
				}
				mergeLists3D(points,temppts,zratio);
			}
		}
		//finally we enforce the edge buffer and output the points
		return enforce3DEdges(points,zedgebuf,data.length);
	}
	
	public void mergeLists3D(List<float[]> master,float[][] source,float zratio) {
		//here we add the points in the source list to the master list
		//overlaps must be reconciled by max intensity
		// stats are centerx,centery,integral,area,maxx,maxy,maxval
		//note that the centerx, centery values might be center of mass
		//assume that the sourcelist is sorted from high intensity to low
		if(master.size()==0) {
			//if the master list is empty
			for(int i=0;i<source.length;i++) master.add(source[i]);
			return;
		}
		//float mindist2=minsep*minsep; //calculate the squared value to avoid square roots later
		int searchr=(int)(0.5f*minsep+0.5f);
		int searchr2=searchr*searchr; //calculate the squared value to avoid square roots later
		int oldsize=master.size();
		for(int i=0;i<source.length;i++) {
			float x=source[i][5];
			float y=source[i][6];
			float z=source[i][7];
			//for each point start by finding the closest master point
			float[] temp0=master.get(0);
			float closestdist2=(x-temp0[5])*(x-temp0[5])+(y-temp0[6])*(y-temp0[6])+(z-temp0[7])*(z-temp0[7])/(zratio*zratio);
			int closestpoint=0;
			for(int j=1;j<oldsize;j++) {
				float[] temp=master.get(j);
				float dist2=(x-temp[5])*(x-temp[5])+(y-temp[6])*(y-temp[6])+(z-temp[7])*(z-temp[7])/(zratio*zratio);
				if(dist2<=closestdist2) {
					closestdist2=dist2;
					closestpoint=j;
				}
			}
			//now that we found the closest point, see if it is too close
			if(closestdist2<searchr2) {
				//this point is too close, keep the one with higher intensity
				if(source[i][8]>master.get(closestpoint)[8]) {
					master.set(closestpoint,source[i]);
				}
			} else {
				//this point doesn't conflict with the original list, add it to the list at the end
				master.add(source[i]);
			}
		}
		return;
	}
	
	public float[] getRegionPad(Object image,int xpos,int ypos,int rwidth,int rheight,int width,int height){
		//this method gets a region of a image and pads it if it is clipped
		if(xpos>=width) return null;
		if(ypos>=height) return null;
		if(xpos<-rwidth) return null;
		if(ypos<-rheight) return null;
		int xstart=xpos;
		if(xstart<0) xstart=0;
		int xend=xpos+rwidth-1;
		if(xend>=width) xend=width-1;
		int trwidth=xend-xstart+1;
		int ystart=ypos;
		if(ystart<0) ystart=0;
		int yend=ypos+rheight-1;
		if(yend>=height) yend=height-1;
		int trheight=yend-ystart+1;
		float[] temp=algutils.get_region2(image,xstart,ystart,trwidth,trheight,width,height);
		//this is technically not quite right because we pad down and to the right but it works for this code
		return padImage(temp,trwidth,trheight,rwidth,rheight);
	}
	
	public float[][] getRegionPad(Object[] image,int xpos,int ypos,int rwidth,int rheight,int width,int height){
		//this method gets a region of a image and pads it if it is clipped
		if(xpos>=width) return null;
		if(ypos>=height) return null;
		if(xpos<-rwidth) return null;
		if(ypos<-rheight) return null;
		int xstart=xpos;
		if(xstart<0) xstart=0;
		int xend=xpos+rwidth-1;
		if(xend>=width) xend=width-1;
		int trwidth=xend-xstart+1;
		int ystart=ypos;
		if(ystart<0) ystart=0;
		int yend=ypos+rheight-1;
		if(yend>=height) yend=height-1;
		int trheight=yend-ystart+1;
		float[][] padded=new float[image.length][];
		for(int i=0;i<image.length;i++){
			float[] temp=algutils.get_region2(image[i],xstart,ystart,trwidth,trheight,width,height);
			//this is technically not quite right because we pad down and to the right but it works for this code
			padded[i]=padImage(temp,trwidth,trheight,rwidth,rheight);
		}
		return padded;
	}
	
	public float[] padImage(float[] image,int width,int height,int newwidth,int newheight) {
		//this method pads simply by copying the image into the upper left hand corner of a blank image
		//make sure newwidth and newheight are greater than the old ones
		float[] newimage=new float[newwidth*newheight];
		for(int i=0;i<height;i++) {
			System.arraycopy(image,i*width,newimage,i*newwidth,width);
		}
		return newimage;
	}
	
	public float[][] padImage(float[][] image,int width,int height,int newwidth,int newheight) {
		//this method pads simply by copying the image into the upper left hand corner of a blank image
		//make sure newwidth and newheight are greater than the old ones
		float[][] newimage=new float[image.length][];
		for(int i=0;i<image.length;i++) newimage[i]=padImage(image[i],width,height,newwidth,newheight);
		return newimage;
	}
	
	public float[][] enforce3DEdges(List<float[]> points,float zedgebuf,int slices){
		//here we enforce the edge buffer and output the points
		List<float[]> filtered=new ArrayList<float[]>();
		for(int i=0;i<points.size();i++) {
			float[] temp=points.get(i);
			if(temp[5]<edgebuf || temp[6]<edgebuf || temp[5]>(width-edgebuf-1) || temp[6]>(height-edgebuf-1) || temp[7]<zedgebuf || temp[7]>(slices-zedgebuf-1)) {}
			//if(temp[5]<edgebuf || temp[6]<edgebuf || temp[5]>(width-edgebuf-1) || temp[6]>(height-edgebuf-1)) {}
			else{filtered.add(temp);}
		}
		float[][] outpoints=new float[filtered.size()][];
		for(int i=0;i<filtered.size();i++) outpoints[i]=filtered.get(i);
		return outpoints;
	}
	
	public float[][] enforceEdges(List<float[]> points){
		//here we enforce the edge buffer and output the points
		List<float[]> filtered=new ArrayList<float[]>();
		for(int i=0;i<points.size();i++) {
			float[] temp=points.get(i);
			if(temp[4]<edgebuf || temp[5]<edgebuf || temp[4]>(width-edgebuf-1) || temp[5]>(height-edgebuf-1)) {}
			else{filtered.add(temp);}
		}
		float[][] outpoints=new float[filtered.size()][];
		for(int i=0;i<filtered.size();i++) outpoints[i]=filtered.get(i);
		return outpoints;
	}

}
