/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.process.FloatProcessor;
import jalgs.interpolation;
import jalgs.jfft.crosscorr2D;
import jalgs.jseg.jsmooth;

import java.awt.Polygon;

public class STICS_map{
	public crosscorr2D cc2D;
	public int subsize,stepsize,xregions,yregions,shift,linewidth;
	public float[][] velocities;

	public STICS_map(int subsize,int stepsize){
		cc2D=new crosscorr2D(subsize,subsize);
		this.subsize=subsize;
		this.stepsize=stepsize;
		linewidth=2;
	}

	public void update_STICS_map(float[][] svel,int xregions,int yregions,float frametime,float scaling){
		// here we reconstruct this object from a scaled velocity map
		this.xregions=xregions;
		this.yregions=yregions;
		this.shift=1;
		// have to unscale the velocities
		velocities=new float[3][xregions*yregions];
		for(int i=0;i<xregions*yregions;i++){
			velocities[0][i]=svel[0][i]*frametime*shift/scaling;
			velocities[1][i]=svel[1][i]*frametime*shift/scaling;
			velocities[2][i]=(float)Math.sqrt(velocities[0][i]*velocities[0][i]+velocities[1][i]*velocities[1][i]);
		}
	}

	public void update_STICS_map(Object[] image,int width,int height,int slices,Polygon roi,int shift){
		update_STICS_map(image,width,height,0,slices,roi,shift);
	}

	public void update_STICS_map(Object[] image,int width,int height,int startslice,int slices1,Polygon roi,int shift){
		boolean[][] mask=new boolean[1][width*height];
		if(roi==null){
			for(int i=0;i<width*height;i++)
				mask[0][i]=true;
		}else{
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(roi.contains(j,i))
						mask[0][j+i*width]=true;
				}
			}
		}
		update_STICS_map(image,width,height,startslice,slices1,mask,shift);
	}

	public void update_STICS_map(Object[] image,int width,int height,int startslice,int slices1,boolean[][] mask,int shift){
		this.shift=shift;
		int slices=slices1;
		if((startslice+slices)>image.length)
			slices=image.length-startslice;
		xregions=1+(int)(((float)width-(float)subsize)/(float)stepsize);
		int newwidth=xregions*subsize;
		yregions=1+(int)(((float)height-(float)subsize)/(float)stepsize);
		int newheight=yregions*subsize;
		int counter=0;
		int minsize=(int)(0.25f*subsize*subsize);
		float[] STICS=new float[newwidth*newheight];
		int mi1=0;
		int mi2=0;
		for(int m=0;m<yregions;m++){
			int ypos=m*stepsize;
			for(int l=0;l<xregions;l++){
				int xpos=l*stepsize;
				for(int k=startslice;k<(startslice+slices-shift);k++){
					// copy the subarray for correlation
					float[] subarray1=new float[subsize*subsize];
					float[] subarray2=new float[subsize*subsize];
					if(mask.length>1){
						mi1=k;
						mi2=k+shift;
					}
					float avg1=0.0f;
					float avg2=0.0f;
					int count=0;
					for(int i=ypos;i<(ypos+subsize);i++){
						for(int j=xpos;j<(xpos+subsize);j++){
							if(mask[mi1][j+i*width]&&mask[mi2][j+i*width]){
								subarray1[(i-ypos)*subsize+(j-xpos)]=((float[])image[k])[i*width+j];
								subarray2[(i-ypos)*subsize+(j-xpos)]=((float[])image[k+shift])[i*width+j];
								avg1+=subarray1[(i-ypos)*subsize+(j-xpos)];
								avg2+=subarray2[(i-ypos)*subsize+(j-xpos)];
								count++;
							}
						}
					}
					if(count>=minsize&&!is_flat(subarray1)){
						// if pixels are outside the roi, set them to the
						// average values
						if(count<subsize*subsize){
							avg1/=(float)count;
							avg2/=(float)count;
							for(int i=ypos;i<(ypos+subsize);i++){
								for(int j=xpos;j<(xpos+subsize);j++){
									if(!(mask[mi1][j+i*width]&&mask[mi2][j+i*width])){
										subarray1[(i-ypos)*subsize+(j-xpos)]=avg1;
										subarray2[(i-ypos)*subsize+(j-xpos)]=avg2;
									}
								}
							}
						}
						subarray1=cc2D.docrosscorr2D(subarray2,subarray1,true,true,false,false);
						// copy the cc to the correlation array
						for(int i=0;i<subsize;i++){
							for(int j=0;j<subsize;j++){
								STICS[(i+m*subsize)*newwidth+j+l*subsize]+=subarray1[i*subsize+j]/(float)(slices-shift);
							}
						}
					}
				}
				counter++;
				IJ.showProgress(counter,xregions*yregions);
			}
		}
		// now find the velocities
		STICS=jsmooth.smooth2D(STICS,newwidth,newheight);
		velocities=new float[3][xregions*yregions]; // record x vel,y vel,
		// and mag
		for(int m=0;m<yregions;m++){
			for(int l=0;l<xregions;l++){
				// first copy the appropriate subarray
				float[] subarray=new float[subsize*subsize];
				for(int i=0;i<subsize;i++){
					for(int j=0;j<subsize;j++){
						subarray[j+i*subsize]=STICS[(i+m*subsize)*newwidth+j+l*subsize];
					}
				}
				// now find the max pixel
				int maxj=subsize/2;
				int maxi=subsize/2;
				float max=subarray[subsize/2+subsize*subsize/2];
				for(int i=0;i<subsize;i++){
					for(int j=0;j<subsize;j++){
						float dist=(float)Math.sqrt((i-subsize/2)*(i-subsize/2)+(j-subsize/2)*(j-subsize/2));
						if(dist<((float)subsize/4.0f)){
							if(subarray[j+i*subsize]>max){
								max=subarray[j+i*subsize];
								maxj=j;
								maxi=i;
							}
						}
					}
				}
				if(maxi<=0){
					maxi=1;
				}
				if(maxi>=(subsize-1)){
					maxi=subsize-2;
				}
				if(maxj<=0){
					maxj=1;
				}
				if(maxj>=(subsize-1)){
					maxj=subsize-2;
				}
				// get the interpolated maximum
				float[] tempmax=interpolation.get_local_max(subarray,maxj,maxi,subsize,subsize);
				float xsum=tempmax[0];
				float ysum=tempmax[1];
				xsum-=0.5f*(float)subsize;
				ysum-=0.5f*(float)subsize;
				float mag=(float)Math.sqrt(xsum*xsum+ysum*ysum);
				velocities[0][m*xregions+l]=xsum;
				velocities[1][m*xregions+l]=ysum;
				velocities[2][m*xregions+l]=mag;
			}
		}
	}

	public FloatProcessor get_map(float scaling,float frametime,int spacing,boolean center,boolean norm,float multiplier,int length,float magthresh){
		FloatProcessor fp=new FloatProcessor(spacing*(xregions+1),spacing*(yregions+1));
		fp.setLineWidth(linewidth);
		for(int m=0;m<yregions;m++){
			int ypos=m*spacing+spacing/2;
			for(int l=0;l<xregions;l++){
				int xpos=l*spacing+spacing/2;
				float mag=velocities[2][l+m*xregions];
				float xvel=velocities[0][l+m*xregions];
				float yvel=velocities[1][l+m*xregions];
				float scaled_mag=(mag/(float)shift)*scaling/frametime;
				if(mag>magthresh){
					fp.setValue(scaled_mag);
					float xnorm=multiplier*xvel;
					float ynorm=multiplier*yvel;
					if(norm){
						xnorm=(float)length*xvel/mag;
						ynorm=(float)length*yvel/mag;
					}
					if(center){
						jutils.draw_arrow(fp,xpos+spacing/2-(int)(0.5f*xnorm),ypos+spacing/2-(int)(0.5f*ynorm),xpos+spacing/2+(int)(0.5f*xnorm),ypos+spacing/2+(int)(0.5f*ynorm));
					}else{
						jutils.draw_arrow(fp,xpos+spacing/2,ypos+spacing/2,xpos+spacing/2+(int)xnorm,ypos+spacing/2+(int)ynorm);
					}
				}
			}
		}
		return fp;
	}

	public float[][] get_scaled_velocities(float scaling,float frametime,int spacing){
		float[][] output=new float[5][xregions*yregions];
		for(int m=0;m<yregions;m++){
			int ypos=m*spacing+spacing/2;
			for(int l=0;l<xregions;l++){
				int xpos=l*spacing+spacing/2;
				float mag=velocities[2][l+m*xregions];
				float xvel=velocities[0][l+m*xregions];
				float yvel=velocities[1][l+m*xregions];
				output[0][l+m*xregions]=(xvel/(float)shift)*scaling/frametime;
				output[1][l+m*xregions]=(yvel/(float)shift)*scaling/frametime;
				output[2][l+m*xregions]=(mag/(float)shift)*scaling/frametime;
				output[3][l+m*xregions]=(float)xpos+0.5f*(float)spacing;
				output[3][l+m*xregions]*=scaling;
				output[4][l+m*xregions]=(float)ypos+0.5f*(float)spacing;
				output[4][l+m*xregions]*=scaling;
			}
		}
		return output;
	}

	public boolean is_flat(float[] arr){
		for(int i=1;i<arr.length;i++){
			if(arr[i]!=arr[0]){
				return false;
			}
		}
		return true;
	}

}
