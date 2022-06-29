/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jalgs.*;
import jguis.*;
import jalgs.jseg.*;

public class find_circular_background_jru_v1 implements PlugIn {
	//here we find a probable circular background region and place it on the image
	//works on the current slice and frame, need to eliminate transmitted light before using
	//start by gaussian smoothing each channel
	//then subtract the min
	//then sum the channels
	//finally find the circle with the minimum average value

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Smooth stdev",1.0,5,15,null);
		gd.addNumericField("Circle Diameter",30,0);
		gd.addNumericField("Border",50,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float smoothstdev=(float)gd.getNextNumber();
		int dia=(int)gd.getNextNumber();
		int border=(int)gd.getNextNumber();
		int currslice=imp.getZ();
		int currframe=imp.getT();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		int nchans=imp.getNChannels();
		Object[] cstack=jutils.get3DCSeries(stack,currslice-1,currframe-1,frames,slices,nchans);
		float[] summed=new float[width*height];
		float[] zeromask=new float[width*height];
		for(int i=0;i<cstack.length;i++){
			float[] smoothed=algutils.convert_arr_float(cstack[i]);
			//start by filling image zeros in with the not0min value from elsewhere (eliminate dark regions)
			float minn0=jstatistics.getstatistic("Not0Min",smoothed,null);
			for(int j=0;j<width*height;j++){
				if(smoothed[j]==0.0f) zeromask[j]=1.0f;
				smoothed[j]=Math.max(minn0,smoothed[j]);
			}
			jsmooth.blur2D(smoothed,smoothstdev,width,height);
			float min=jstatistics.getstatistic("Min",smoothed,null);
			//float max=jstatistics.getstatistic("Max",smoothed,null);
			for(int j=0;j<width*height;j++){
				summed[j]+=(smoothed[j]-min);
			}
		}
		float[] smoothedmask=jsmooth.smooth2DCircle(zeromask,width,height,dia,"Max");
		for(int i=0;i<width*height;i++){
			if(smoothedmask[i]>0.0f) zeromask[i]=1.0f;
		}
		float[] smoothed=jsmooth.smooth2DCircle(summed,width,height,dia,"Avg");
		//find the minimum coordinates
		int minx=border; int miny=border; float minval=jstatistics.getstatistic("Max",smoothed,null);
		for(int i=border;i<(height-border);i++){
			for(int j=border;j<(width-border);j++){
				if(zeromask[j+i*width]==0.0f && smoothed[j+i*width]<minval){
					minval=smoothed[j+i*width];
					minx=j; miny=i;
				}
			}
		}
		//now put the circle on the image
		imp.setRoi(new OvalRoi(minx-dia/2,miny-dia/2,dia,dia));
	}

}
