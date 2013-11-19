/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class jhough_transform_21{
	// this version of the hough transform uses Canny edge detection and the
	// sobel gradient
	// to draw lines normal to each edge point from minradius to maxradius and
	// -minradius to -maxradius
	// these lines are accumulated and threshholded to find the circle centers
	// the circle radii are then found by finding the radius with the most edge
	// points on it
	// I have yet to find an application where it works well
	public int minrad,maxrad,width,height,accumbin;
	public float min_dist,canny_thresh,circle_thresh;
	int[][][] pixelshifts;
	int[] nangles;

	public jhough_transform_21(int minrad1,int maxrad1,float canny_thresh1,float circle_thresh1,float min_dist1,int width1,int height1,int accumbin1){
		minrad=minrad1;
		maxrad=maxrad1;
		canny_thresh=canny_thresh1;
		circle_thresh=circle_thresh1;
		min_dist=min_dist1;
		width=width1;
		height=height1;
		accumbin=accumbin1;
		setup_pixel_shifts();
	}

	private void setup_pixel_shifts(){
		pixelshifts=new int[maxrad-minrad+1][2][];
		nangles=new int[maxrad-minrad+1];
		for(int radius=minrad;radius<=maxrad;radius++){
			int rindex=radius-minrad;
			double dangle=1.0/(double)radius;
			nangles[rindex]=(int)(2.0*Math.PI*(double)radius);
			pixelshifts[rindex][0]=new int[nangles[rindex]];
			pixelshifts[rindex][1]=new int[nangles[rindex]];
			for(int i=0;i<nangles[rindex];i++){
				double x=(double)radius*Math.cos(dangle*(double)i);
				double y=(double)radius*Math.sin(dangle*(double)i);
				pixelshifts[rindex][0][i]=(int)x;
				pixelshifts[rindex][1][i]=(int)y;
			}
		}
	}

	public int[] get_accum(float[] image){
		// start by doing the canny edge detection
		float[][] dxdy=(new jsobel(width,height)).do_sobel2(image);
		float[][] mag_angle=(new jsobel(width,height)).get_mag_angle(dxdy[0],dxdy[1]);
		byte[] edges=(new canny(canny_thresh,canny_thresh/2.0f,width,height)).find_edges(image,mag_angle);
		int[] accum=new int[width*height];
		// draw lines from each edge point along the edge gradient
		for(int i=0;i<height;i++){
			for(int k=0;k<width;k++){
				int index=k+i*width;
				if(edges[index]==(byte)255){
					float xinc=dxdy[0][index]/mag_angle[0][index];
					float yinc=dxdy[1][index]/mag_angle[0][index];
					for(int j=minrad;j<=maxrad;j++){
						int x=k+(int)(j*xinc);
						int y=i+(int)(j*yinc);
						if(x>=0&&x<width&&y>=0&&y<height){
							accum[x+y*width]++;
						}else{
							break;
						}
					}
					for(int j=-minrad;j>=-maxrad;j--){
						int x=k+(int)(j*xinc);
						int y=i+(int)(j*yinc);
						if(x>=0&&x<width&&y>=0&&y<height){
							accum[x+y*width]++;
						}else{
							break;
						}
					}
				}
			}
		}
		if(accumbin>1){
			bin_image(accum);
		}
		return accum;
	}

	public Object[] detect_circles(float[] image,boolean getradii,boolean drawcircles){
		// start by doing the canny edge detection
		float[][] dxdy=(new jsobel(width,height)).do_sobel2(image);
		float[][] mag_angle=(new jsobel(width,height)).get_mag_angle(dxdy[0],dxdy[1]);
		byte[] edges=(new canny(canny_thresh,canny_thresh/2.0f,width,height)).find_edges(image,mag_angle);
		int[] accum=new int[width*height];
		// draw lines from each edge point along the edge gradient
		for(int i=0;i<height;i++){
			for(int k=0;k<width;k++){
				int index=k+i*width;
				if(edges[index]==(byte)255){
					float xinc=dxdy[0][index]/mag_angle[0][index];
					float yinc=dxdy[1][index]/mag_angle[0][index];
					for(int j=minrad;j<=maxrad;j++){
						int x=k+(int)(j*xinc);
						int y=i+(int)(j*yinc);
						if(x>=0&&x<width&&y>=0&&y<height){
							accum[x+y*width]++;
						}else{
							break;
						}
					}
					for(int j=-minrad;j>=-maxrad;j--){
						int x=k+(int)(j*xinc);
						int y=i+(int)(j*yinc);
						if(x>=0&&x<width&&y>=0&&y<height){
							accum[x+y*width]++;
						}else{
							break;
						}
					}
				}
			}
		}
		if(accumbin>1){
			bin_image(accum);
		}
		int[] accum1=accum.clone();
		// now threshhold the accumulated image to get the centers and do max
		// not mask to get unique centers
		int[][] centers=new int[10000][2];
		int max=0;
		int counter=0;
		do{
			int[] temp=maxarray(accum);
			max=accum[temp[0]+width*temp[1]];
			if(max>(int)circle_thresh){
				centers[counter][0]=temp[0];
				centers[counter][1]=temp[1];
				clear_circle(accum,centers[counter]);
				counter++;
			}
		}while(max>(int)circle_thresh&&counter<10000);
		int[][] centers1=new int[counter][2];
		for(int i=0;i<counter;i++){
			centers1[i][0]=centers[i][0];
			centers1[i][1]=centers[i][1];
		}
		centers=null;
		// find the radii of the circles if desired
		int[] radii=null;
		byte[] retimage=null;
		if(getradii){
			radii=new int[counter];
			for(int i=0;i<counter;i++){
				int bestrad=minrad;
				float bestint=circle_integral(edges,centers1[i],minrad);
				for(int j=minrad+1;j<=maxrad;j++){
					float integral=circle_integral(edges,centers1[i],j);
					if(integral>bestint){
						bestint=integral;
						bestrad=j;
					}
				}
				radii[i]=bestrad;
			}
			if(drawcircles){
				retimage=new byte[width*height];
				for(int i=0;i<counter;i++){
					retimage[centers1[i][0]+width*centers1[i][1]]=(byte)255;
					// draw the circles on the output image
					draw_circle(retimage,centers1[i],radii[i]);
				}
			}
		}

		// dxdy is a 2 by width*height array giving the gradient vectors from
		// the sobel operator
		// edges is canny detected width*height byte array image
		// accum is a width*height integer array with the hough lines
		// centers1 is a n by 2 array with the detected centers (padded with
		// zeros)
		// radii is an int array with the detected radii (or null)
		// retimage is a byte array width*height image with circles drawn on it
		// (or null)
		Object[] retvals=new Object[]{dxdy,edges,accum1,centers1,radii,retimage};
		return retvals;
	}

	private void bin_image(int[] arr){
		for(int i=0;i<height;i+=accumbin){
			for(int j=0;j<width;j+=accumbin){
				int temp=0;
				for(int k=0;k<accumbin;k++){
					for(int l=0;l<accumbin;l++){
						temp+=arr[j+l+(k+i)*width];
					}
				}
				for(int k=0;k<accumbin;k++){
					for(int l=0;l<accumbin;l++){
						arr[j+l+(k+i)*width]=temp;
					}
				}
			}
		}
	}

	private float circle_integral(byte[] image,int[] center,int radius){
		int radindex=radius-minrad;
		float integral=0.0f;
		for(int i=0;i<nangles[radindex];i++){
			int x=center[0]+pixelshifts[radindex][0][i];
			int y=center[1]+pixelshifts[radindex][1][i];
			if(x>=0&&x<width&&y>=0&&y<height){
				if(image[x+y*width]==(byte)255){
					integral+=1.0f/(float)nangles[radindex];
				}
			}
		}
		return integral;
	}

	private void clear_circle(int[] image,int[] center){
		int min_dist_sq=(int)(min_dist*min_dist);
		int imin_dist=(int)min_dist;
		for(int i=-imin_dist;i<=imin_dist;i++){
			for(int j=-imin_dist;j<=imin_dist;j++){
				if((j*j+i*i)<=min_dist_sq){
					int x=j+center[0];
					int y=i+center[1];
					if(y>=0&&y<height&&x>=0&&x<width){
						image[x+y*width]=0;
					}
				}
			}
		}
	}

	private void draw_circle(byte[] image,int[] center,int radius){
		int radindex=radius-minrad;
		for(int i=0;i<nangles[radindex];i++){
			int x=center[0]+pixelshifts[radindex][0][i];
			int y=center[1]+pixelshifts[radindex][1][i];
			if(x>=0&&x<width&&y>=0&&y<height){
				image[x+y*width]=(byte)255;
			}
		}
	}

	private int[] maxarray(int[] arr){
		int max=arr[0];
		int xmax=0;
		int ymax=0;
		for(int i=1;i<height-1;i++){
			for(int j=1;j<width-1;j++){
				int temp=j+i*width;
				if(arr[temp]>max){
					max=arr[temp];
					xmax=j;
					ymax=i;
				}
			}
		}
		int[] temp2={xmax,ymax};
		return temp2;
	}

}
