/*******************************************************************************
 * Copyright (c) 2021 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.filter.*;
import jalgs.jseg.*;

public class watershed_subset_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin applies watershed segmentation to a subset of objects in an image
		//the subset is based on area and circularity
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_Watershed_Area (pix)",1000,0);
		gd.addNumericField("Max_Watershed_Area (pix)",1000000,0);
		gd.addNumericField("Min_Circ (0-1)",0.0,5,15,null);
		gd.addNumericField("Max_Circ (0-1)",0.9,5,15,null);
		gd.addNumericField("Min_Resulting_Area (pix)",500,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float minarea=(float)gd.getNextNumber();
		float maxarea=(float)gd.getNextNumber();
		float mincirc=(float)gd.getNextNumber();
		float maxcirc=(float)gd.getNextNumber();
		float minresarea=(float)gd.getNextNumber();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack ostack=new ImageStack(width,height);
		for(int i=0;i<size;i++){
			byte[] mask=(byte[])stack.getProcessor(i+1).getPixels();
			byte[] outmask=watershedSubset(mask,width,height,minarea,maxarea,mincirc,maxcirc,minresarea);
			ostack.addSlice("",outmask);
		}
		new ImagePlus("Watershed_Subset",ostack).show();
	}
	
	public byte[] watershedSubset(byte[] mask,int width,int height,float minarea,float maxarea,float mincirc,float maxcirc,float minresarea){
		findblobs3 fb=new findblobs3(width,height);
		float[] objects=fb.dofindblobs(mask);
		float[] filtobjects=objects.clone();
		//int[][] filllims=fb.getallfilllimits(objects);
		fb.filter_area_circ(filtobjects,new float[]{minarea,maxarea,mincirc,maxcirc});
		//new ImagePlus("filtered",new ByteProcessor(width,height,fb.tobinary(filtobjects,false))).show();
		//get the non-filtered objects as well for later
		float[] nonfiltobjects=new float[width*height];
		for(int i=0;i<width*height;i++){
			if(filtobjects[i]==0.0f){
				nonfiltobjects[i]=objects[i];
			}
		}
		//new ImagePlus("non-filtered",new ByteProcessor(width,height,fb.tobinary(nonfiltobjects,false))).show();
		//now run watershed
		ImageProcessor filtip=new ByteProcessor(width,height,fb.tobinary(filtobjects,false));
		(new EDM()).toWatershed(filtip);
		//new ImagePlus("watershed",filtip).show();
		//now copy these changes back to the filtered objects
		byte[] watershed=(byte[])filtip.getPixels();
		float[] filtobjects2=filtobjects.clone();
		for(int i=0;i<watershed.length;i++){
			if((watershed[i]&0xffff)<255){
				filtobjects2[i]=0.0f;
			}
		}
		fb.renumber_objects2(filtobjects2);
		//new ImagePlus("watershed objects",new FloatProcessor(width,height,filtobjects2.clone(),null)).show();
		//need to find the objects that are below our result area threshold
		int[] areas=fb.get_areas(filtobjects2);
		int ntoosmall=0;
		int[] toosmallids=new int[areas.length];
		for(int i=0;i<areas.length;i++){
			if((float)areas[i]<minresarea){
				toosmallids[ntoosmall]=i+1;
				ntoosmall++;
			}
			//IJ.log(""+areas[i]);
		}
		IJ.log(""+ntoosmall+" objects too small after watershed");
		//now go back through and find which old unsplit objects correspond to the new too small split ones
		int[] toosmalloldids=new int[ntoosmall];
		for(int i=0;i<width*height;i++){
			for(int j=0;j<ntoosmall;j++){
				if(toosmalloldids[j]==0 && filtobjects2[i]==toosmallids[j]){
					toosmalloldids[j]=(int)filtobjects[i];
					//IJ.log("unsplit id "+toosmalloldids[j]+","+toosmallids[j]);
					break;
				}
			}
		}
		//finally unsplit the too small ones
		for(int i=0;i<width*height;i++){
			if(filtobjects[i]>0.0f){
				for(int j=0;j<ntoosmall;j++){
					if((int)filtobjects[i]==toosmalloldids[j]){
						filtobjects2[i]=filtobjects[i];
						break;
					}
				}
			}
		}
		//new ImagePlus("corrected",new FloatProcessor(width,height,filtobjects2,null)).show();
		//now copy back into the non-filtered objects image
		for(int i=0;i<width*height;i++){
			if(filtobjects2[i]>0.0f){
				nonfiltobjects[i]+=filtobjects2[i];
			}
		}
		return fb.tobinary(nonfiltobjects,false);
	}

}
