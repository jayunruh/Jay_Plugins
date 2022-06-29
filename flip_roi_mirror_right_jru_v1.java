/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class flip_roi_mirror_right_jru_v1 implements PlugIn {
	//this plugin takes an roi in the current channel as well as a horizontally mirrored roi
	//if the rightmost roi is higher, the image is left alone
	//otherwise, the image is horizontally flipped
	//assume the original roi is on the right side of the image

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		Roi roi=imp.getRoi();
		Rectangle r=roi.getBounds();
		float center=0.5f*(float)width;
		float ur=(float)(r.x+r.width);
		float neworig=center-(ur-center);
		Rectangle r2=new Rectangle((int)neworig,r.y,r.width,r.height);
		int nch=imp.getNChannels();
		int currchan=imp.getChannel();
		ImageStack stack=imp.getStack();
		int totslices=stack.getSize();
		for(int i=0;i<totslices/nch;i++){
			float right=jstatistics.getstatistic("Avg",stack.getPixels(i*nch+currchan),width,height,r,null);
			float left=jstatistics.getstatistic("Avg",stack.getPixels(i*nch+currchan),width,height,r2,null);
			if(left>right){
				for(int j=0;j<nch;j++){
					flip_hor(stack.getPixels(i*nch+j+1),width,height);
				}
				IJ.log("flipping slice"+i);
			}
			IJ.showProgress(i,totslices/nch);
		}
		imp.updateAndDraw();
	}

	public void flip_hor(Object image,int width,int height){
		if(image instanceof float[]) flip_hor((float[])image,width,height);
		if(image instanceof short[]) flip_hor((short[])image,width,height);
		if(image instanceof byte[]) flip_hor((byte[])image,width,height);
	}

	public void flip_hor(float[] image,int width,int height){
		float[] temp=image.clone();
		for(int i=0;i<height;i++){
			int cnt=0;
			for(int j=(width-1);j>=0;j--){
				image[cnt+i*width]=temp[j+i*width]; cnt++;
			}
		}
	}

	public void flip_hor(short[] image,int width,int height){
		short[] temp=image.clone();
		for(int i=0;i<height;i++){
			int cnt=0;
			for(int j=(width-1);j>=0;j--){
				image[cnt+i*width]=temp[j+i*width]; cnt++;
			}
		}
	}

	public void flip_hor(byte[] image,int width,int height){
		byte[] temp=image.clone();
		for(int i=0;i<height;i++){
			int cnt=0;
			for(int j=(width-1);j>=0;j--){
				image[cnt+i*width]=temp[j+i*width]; cnt++;
			}
		}
	}

	public void flip_hor(int[] image,int width,int height){
		int[] temp=image.clone();
		for(int i=0;i<height;i++){
			int cnt=0;
			for(int j=(width-1);j>=0;j--){
				image[cnt+i*width]=temp[j+i*width]; cnt++;
			}
		}
	}

}
