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
import jguis.*;

public class dilate_color_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<stack.getSize();i++){
			int[] pix1=(int[])stack.getPixels(i+1);
			int[] dil=dilateColorSlice(pix1,width,height);
			retstack.addSlice("",dil);
		}
		new ImagePlus("Color Dilated",retstack).show();
	}

	public int[] dilateColorSlice(int[] pix1,int width,int height){
		int[] pix=pix1.clone();
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				if(pix[i]==0xff000000){
					int[] neighbors=algutils.getNeighbors(pix1,j,i,width,height);
					//count the colors (if any) surrounding our pixel
					//dilate the most common
					int[][] hist=new int[3][256];
					boolean dil=false;
					for(int k=0;k<neighbors.length;k++){
						if(neighbors[k]!=0xff000000){
							int[] rgb=jutils.intval2rgb(neighbors[k]);
							hist[0][rgb[0]]++;
							hist[1][rgb[1]]++;
							hist[2][rgb[2]]++;
							dil=true;
						}
					}
					if(dil){
						int maxr=maxpos(hist[0]);
						int maxg=maxpos(hist[1]);
						int maxb=maxpos(hist[2]);
						pix[j+i*width]=jutils.rgb2intval(maxr,maxg,maxb);
					}
				}
			}
		}
		return pix;
	}

	public int maxpos(int[] arr){
		int max=arr[arr.length-1];
		int maxpos=arr.length-1;
		for(int i=(arr.length-1);i>=0;i--){
			if(arr[i]>max){
				maxpos=i; max=arr[i];
			}
		}
		return maxpos;
	}
}
