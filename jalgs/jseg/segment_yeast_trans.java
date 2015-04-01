/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import jalgs.jstatistics;

public class segment_yeast_trans{
	public static byte[] segment_image(float[] image,int width,int height,float sobel_thresh,int minarea,int maxarea){
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		return segment_sobel(smoothed,width,height,sobel_thresh,minarea,maxarea);
	}

	public static byte[] segment_image2(float[] image,int width,int height,float sobel_thresh,int minarea,int maxarea){
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		float[] percent={sobel_thresh};
		float thresh=jstatistics.getstatistic("Percentile",smoothed,percent);
		return segment_sobel(smoothed,width,height,thresh,minarea,maxarea);
	}

	public static byte[] segment_sobel(float[] sobel,int width,int height,float thresh,int minarea,int maxarea){
		findblobs3 fb=new findblobs3(width,height);
		binary_processing bp=new binary_processing(width,height);
		float[] skel=fb.dofindblobs(sobel,thresh);
		fb.clear_edges(skel);
		fb.closeobjects(skel);
		fb.closeobjects(skel);
		fb.skeletonize(skel);
		fb.dilateobjects(skel);
		byte[] objects=fb.tobinary(skel,false);
		bp.fillholes(objects);
		for(int i=0;i<width*height;i++){
			if(skel[i]>0.0f){
				objects[i]=(byte)0;
			}
		}
		bp.fillholes(objects);
		skel=fb.dofindblobs(objects);
		int[] arealims={minarea,maxarea};
		fb.filter_area(skel,arealims);
		fb.dilateobjects(skel);
		objects=fb.tobinary(skel,true);
		return objects;
	}
}
