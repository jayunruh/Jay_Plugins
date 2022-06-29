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
	
	/*************
	 * this is the default implementation
	 * @param image
	 * @param width
	 * @param height
	 * @param sobel_thresh: a raw sobel (slope) edge threshold
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return : a binary thresholded image
	 */
	public static byte[] segment_image(float[] image,int width,int height,float sobel_thresh,int minarea,int maxarea){
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		return segment_sobel(smoothed,width,height,sobel_thresh,minarea,maxarea);
	}

	/*************
	 * this version segments at a percentile of the sobel intensity distribution
	 * @param image
	 * @param width
	 * @param height
	 * @param sobel_thresh: the percentile to threshold at (0 to 100)
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return : a binary thresholded image
	 */
	public static byte[] segment_image2(float[] image,int width,int height,float sobel_thresh,int minarea,int maxarea){
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		float[] percent={sobel_thresh};
		float thresh=jstatistics.getstatistic("Percentile",smoothed,percent);
		return segment_sobel(smoothed,width,height,thresh,minarea,maxarea);
	}
	
	/***************
	 * this thresholds at a multiple of a statistic (see jalgs.jstatistics) specified by stat
	 * @param image
	 * @param width
	 * @param height
	 * @param sobel_thresh_mult: the statistic multiplier
	 * @param stat: the statistic for the threshold
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return : a binary thresholded image
	 */
	public static byte[] segment_image3(float[] image,int width,int height,float sobel_thresh_mult,String stat,int minarea,int maxarea){
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		float thresh=sobel_thresh_mult*jstatistics.getstatistic(stat,smoothed,null);
		return segment_sobel(smoothed,width,height,thresh,minarea,maxarea);
	}
	
	/********************
	 * this version outputs floating point numbered objects instead of a binary image
	 * @param image
	 * @param width
	 * @param height
	 * @param sobel_thresh_mult: the statistic multiplier
	 * @param stat: the statistic for the threshold
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return : a binary thresholded image
	 */
	public static float[] segment_image4(float[] image,int width,int height,float sobel_thresh_mult,String stat,int minarea,int maxarea){
		//this version outputs the indexed objects
		float[] smoothed=jsmooth.smooth2D(image,width,height);
		smoothed=((new jsobel(width,height)).do_sobel(smoothed))[0];
		float thresh=sobel_thresh_mult*jstatistics.getstatistic(stat,smoothed,null);
		return segment_sobel2(smoothed,width,height,thresh,minarea,maxarea);
	}

	/********************
	 * this segments a sobel image and returns a binary image
	 * @param sobel
	 * @param width
	 * @param height
	 * @param thresh: the raw threshold
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return
	 */
	public static byte[] segment_sobel(float[] sobel,int width,int height,float thresh,int minarea,int maxarea){
		float[] skel=segment_sobel2(sobel,width,height,thresh,minarea,maxarea);
		return (new findblobs3(width,height)).tobinary(skel,true);
	}
	
	/*****************
	 * this segments a sobel image and returns floating point indexed image
	 * @param sobel
	 * @param width
	 * @param height
	 * @param thresh: the raw threshold
	 * @param minarea: min area of segmented objects
	 * @param maxarea: max area of segmented objects
	 * @return
	 */
	public static float[] segment_sobel2(float[] sobel,int width,int height,float thresh,int minarea,int maxarea){
		//this version outputs an indexed objects image
		findblobs3 fb=new findblobs3(width,height);
		binary_processing bp=new binary_processing(width,height);
		//first simply threshold the image and get the indexed objects
		float[] skel=fb.dofindblobs(sobel,thresh);
		fb.clear_edges(skel);
		fb.closeobjects(skel); //binary close
		fb.closeobjects(skel);
		fb.skeletonize(skel);
		fb.dilateobjects(skel); //prevents objects from merging--good idea?
		byte[] objects=fb.tobinary(skel,false);
		bp.fillholes(objects);
		//now remove the skeletons from original image to separate the objects
		for(int i=0;i<width*height;i++){
			if(skel[i]>0.0f){
				objects[i]=(byte)0;
			}
		}
		bp.fillholes(objects);
		//refind the objects
		skel=fb.dofindblobs(objects);
		//filter on size
		int[] arealims={minarea,maxarea};
		fb.filter_area(skel,arealims);
		fb.dilateobjects(skel);
		return skel;
		//objects=fb.tobinary(skel,true);
		//return objects;
	}
}
