/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class jhough_transform{
	int[][][] pixelshifts;
	int[] nangles;
	int startradius,endradius,width,height;

	public jhough_transform(int startrad1,int endrad1,int width1,int height1){
		startradius=startrad1;
		endradius=endrad1;
		width=width1;
		height=height1;
		setup_pixel_shifts();
	}

	private void setup_pixel_shifts(){
		pixelshifts=new int[endradius-startradius+1][2][];
		nangles=new int[endradius-startradius+1];
		for(int radius=startradius;radius<=endradius;radius++){
			int rindex=radius-startradius;
			double dangle=1.0/radius;
			nangles[rindex]=(int)(2.0*Math.PI*radius);
			pixelshifts[rindex][0]=new int[nangles[rindex]];
			pixelshifts[rindex][1]=new int[nangles[rindex]];
			for(int i=0;i<nangles[rindex];i++){
				double x=radius*Math.cos(dangle*i);
				double y=radius*Math.sin(dangle*i);
				pixelshifts[rindex][0][i]=(int)x;
				pixelshifts[rindex][1][i]=(int)y;
			}
		}
	}

	public float[] hough_transform(float[] image,int radius){
		int rindex=radius-startradius;
		float[] newimage=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float integral=0.0f;
				for(int k=0;k<nangles[rindex];k++){
					int x=j+pixelshifts[rindex][0][k];
					if(x<0){
						x+=width;
					}
					if(x>=width){
						x-=width;
					}
					int y=i+pixelshifts[rindex][1][k];
					if(y<0){
						y+=height;
					}
					if(y>=height){
						y-=height;
					}
					integral+=image[x+y*width];
				}
				newimage[j+i*width]=integral/nangles[rindex];
			}
		}
		return newimage;
	}

	public float[] hough_transform(byte[] image,int radius){
		int rindex=radius-startradius;
		float[] newimage=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float integral=0.0f;
				for(int k=0;k<nangles[rindex];k++){
					int x=j+pixelshifts[rindex][0][k];
					if(x<0){
						x+=width;
					}
					if(x>=width){
						x-=width;
					}
					int y=i+pixelshifts[rindex][1][k];
					if(y<0){
						y+=height;
					}
					if(y>=height){
						y-=height;
					}
					if(image[x+y*width]==(byte)255){
						integral+=1.0f;
					}
				}
				newimage[j+i*width]=integral/nangles[rindex];
			}
		}
		return newimage;
	}

	public byte[] inverse_hough_transform(byte[] image,int radius){
		// here we basically create a circle for every image pixel above
		// threshhold
		int rindex=radius-startradius;
		byte[] newimage=new byte[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]>(byte)0){
					for(int k=0;k<nangles[rindex];k++){
						int x=j+pixelshifts[rindex][0][k];
						if(x<0){
							x+=width;
						}
						if(x>=width){
							x-=width;
						}
						int y=i+pixelshifts[rindex][1][k];
						if(y<0){
							y+=height;
						}
						if(y>=height){
							y-=height;
						}
						newimage[x+y*width]=(byte)255;
					}
				}
			}
		}
		return newimage;
	}

	public byte[] inverse_hough_transform(float[] image,float threshhold,int radius){
		// here we basically create a circle for every image pixel above
		// threshhold
		int rindex=radius-startradius;
		byte[] newimage=new byte[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]>threshhold){
					for(int k=0;k<nangles[rindex];k++){
						int x=j+pixelshifts[rindex][0][k];
						if(x<0){
							x+=width;
						}
						if(x>=width){
							x-=width;
						}
						int y=i+pixelshifts[rindex][1][k];
						if(y<0){
							y+=height;
						}
						if(y>=height){
							y-=height;
						}
						newimage[x+y*width]=(byte)255;
					}
				}
			}
		}
		return newimage;
	}

	public void inverse_hough_transform(float[] image,byte[] newimage,float threshhold,int radius){
		// here we basically create a circle for every image pixel above
		// threshhold
		int rindex=radius-startradius;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]>threshhold){
					for(int k=0;k<nangles[rindex];k++){
						int x=j+pixelshifts[rindex][0][k];
						if(x<0){
							x+=width;
						}
						if(x>=width){
							x-=width;
						}
						int y=i+pixelshifts[rindex][1][k];
						if(y<0){
							y+=height;
						}
						if(y>=height){
							y-=height;
						}
						newimage[x+y*width]=(byte)255;
					}
				}
			}
		}
	}

	public float[] do_hough_transform_slice(byte[] image,float threshhold,int startradius1,int endradius1){
		// start by looping through the radii and getting the maximum projection
		// hough transform along with a maximum radius image
		float[] maxproj=new float[width*height];
		int[] radproj=new int[width*height];
		float[] outobjects=new float[width*height];
		for(int i=startradius1;i<=endradius1;i++){
			// get the hough transform
			float[] tpixels=hough_transform(image,i);
			// update the maximum projection
			for(int j=0;j<width*height;j++){
				if(tpixels[j]>maxproj[j]){
					maxproj[j]=tpixels[j];
					radproj[j]=i;
				}
			}
		}
		// now create the object image by using a max>threshhold not mask
		// approach
		float objid=0.0f;
		float max=0.0f;
		do{
			float[] temp=getmax(maxproj);
			max=temp[0];
			int maxx=(int)temp[1];
			int maxy=(int)temp[2];
			if(max>=threshhold){
				objid+=1.0;
				int maxradius=radproj[maxx+width*maxy];
				int startx=maxx-maxradius;
				if(startx<0){
					startx=0;
				}
				int endx=maxx+maxradius;
				if(endx>=width){
					endx=width-1;
				}
				int starty=maxy-maxradius;
				if(starty<0){
					starty=0;
				}
				int endy=maxy+maxradius;
				if(endy>=height){
					endy=height-1;
				}
				for(int i=starty;i<=endy;i++){
					int ydiff=i-maxy;
					for(int j=startx;j<=endx;j++){
						int xdiff=j-maxx;
						double r=Math.sqrt(ydiff*ydiff+xdiff*xdiff);
						if(r<=maxradius){
							maxproj[j+i*width]=0.0f;
							if(outobjects[j+i*width]==0.0f){
								outobjects[j+i*width]=objid;
							}
						}
					}
				}
			}else{
				break;
			}
		}while(max>=threshhold);
		return outobjects;
	}

	public float[] getmax(float[] data){
		float max=data[0];
		int maxx=0;
		int maxy=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(data[j+i*width]>max){
					max=data[j+i*width];
					maxx=j;
					maxy=i;
				}
			}
		}
		float[] ret={max,maxx,maxy};
		return ret;
	}
}
