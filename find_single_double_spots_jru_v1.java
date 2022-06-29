/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;
import jguis.*;

public class find_single_double_spots_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Edge_Border",5,0);
		gd.addNumericField("Guess_Stdev",1.5f,5,15,null);
		gd.addNumericField("Fraction_Thresh",0.1f,5,15,null);
		gd.addNumericField("Max_Dist (pix)",15.0f,5,15,null);
		gd.addCheckbox("Find_Best_Slice",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int border=(int)gd.getNextNumber();
		float guessstdev=(float)gd.getNextNumber();
		float thresh=(float)gd.getNextNumber();
		float maxdist=(float)gd.getNextNumber();
		boolean findslice=gd.getNextBoolean();

		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		int currchan=imp.getC()-1;
		int currframe=imp.getT()-1;
		ImageProcessor ip=imp.getProcessor();
		float[] pix=(float[])ip.convertToFloat().getPixels();
		Rectangle r=ip.getRoi();
		if(imp.getRoi()!=null){
			pix=algutils.get_region2(pix,r.x,r.y,r.width,r.height,width,height);
			if(findslice){
				Object[] zstack=jutils.get3DZSeries(stack,currchan,currframe,nframes,nslices,nchans);
				int maxslice=0;
				float maxprofile=0.0f;
				Object maxpix=null;
				for(int i=0;i<nslices;i++){
					Object roipix=algutils.get_region2(zstack[i],r.x,r.y,r.width,r.height,width,height);
					float avg=jstatistics.getstatistic("Avg",roipix,null);
					if(avg>maxprofile){
						maxprofile=avg; maxslice=i; maxpix=roipix;
					}
				}
				pix=algutils.convert_arr_float2(maxpix);
			}
			width=r.width; height=r.height;
		} else {
			if(findslice){
				Object[] zstack=jutils.get3DZSeries(stack,currchan,currframe,nframes,nslices,nchans);
				int maxslice=0;
				float maxprofile=0.0f;
				Object maxpix=null;
				for(int i=0;i<nslices;i++){
					float[] region=algutils.get_region2(zstack[i],border,border,width-2*border,height-2*border,width,height);
					float avg=jstatistics.getstatistic("Avg",region,null);
					if(avg>maxprofile){
						maxprofile=avg; maxslice=i; maxpix=zstack[i];
					}
				}
				pix=algutils.convert_arr_float2(maxpix);
				IJ.log("max slice = "+maxslice);
			}
		}
		float[][] peaks=findPeaks(pix,width,height,guessstdev,border);
		if(r!=null){
			for(int i=0;i<3;i++){
				peaks[i][0]+=r.x; peaks[i][1]+=r.y;
			}
		}
		IJ.log("Peak 1 = "+peaks[0][0]+" , "+peaks[0][1]+" , "+peaks[0][2]);
		IJ.log("Peak 2 = "+peaks[1][0]+" , "+peaks[1][1]+" , "+peaks[1][2]);
		IJ.log("Peak 3 = "+peaks[2][0]+" , "+peaks[2][1]+" , "+peaks[2][2]);
		float fracval=(peaks[1][2]-peaks[2][2])/(peaks[0][2]-peaks[2][2]);
		float dist=(float)Math.sqrt((peaks[0][0]-peaks[1][0])*(peaks[0][0]-peaks[1][0])+(peaks[0][1]-peaks[1][1])*(peaks[0][1]-peaks[1][1]));
		IJ.log("Relative Peak2 Intensity = "+fracval+" , dist = "+dist);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		rman.addRoi(new PointRoi((int)peaks[0][0],(int)peaks[0][1]));
		if(fracval>=thresh && dist<=maxdist){
			rman.addRoi(new PointRoi((int)peaks[1][0],(int)peaks[1][1]));
		}
	}

	public float[][] findPeaks(float[] image,int width,int height,float guessstdev,int border){
		//here we try to figure out if we have one or two peaks
		//use the max not mask approach
		//the third peak should define the noise level
		//3*guessstdev is the occlusion diameter for each peak
		//thresh is the fraction of the difference between the third and first peaks that constitutes a true peak
		//start by finding the first maximum position
		float[] image1=image.clone();
		boolean[] masked=new boolean[width*height];
		float minint=jstatistics.getstatistic("Min",image1,null);
		int[] firstmax=findMaxPos(image1,width,height,border);
		float maskr=1.5f*guessstdev;
		//occlude the first maximum
		float firstint=maskSpot(image1,width,height,firstmax,maskr,minint,masked);
		//new ImagePlus("firstmask",new FloatProcessor(width,height,image1.clone(),null)).show();
		//now find the second max
		int[] secondmax=findMaxPos(image1,width,height,border);
		float secondint=maskSpot(image1,width,height,secondmax,maskr,minint,masked);
		//new ImagePlus("secondmask",new FloatProcessor(width,height,image1.clone(),null)).show();
		//and the third
		int[] thirdmax=findMaxPos(image1,width,height,border);
		float thirdint=maskSpot(image1,width,height,thirdmax,maskr,minint,masked);
		float[][] output=new float[3][];
		output[0]=new float[]{firstmax[0],firstmax[1],firstint};
		output[1]=new float[]{secondmax[0],secondmax[1],secondint};
		output[2]=new float[]{thirdmax[0],thirdmax[1],thirdint};
		return output;
	}

	public int[] findMaxPos(float[] image,int width,int height,int border){
		float max=image[0];
		int[] maxpos=new int[2];
		for(int i=border;i<(height-border);i++){
			for(int j=border;j<(width-border);j++){
				if(image[i*width+j]>max){
					max=image[i*width+j];
					maxpos=new int[]{j,i};
				}
			}
		}
		return maxpos;
	}

	public float maskSpot(float[] image,int width,int height,int[] center,float maskr,float maskval,boolean[] masked){
		int lowerx=(int)((float)center[0]-maskr);
		int upperx=1+lowerx+(int)(2.0f*maskr);
		int lowery=(int)((float)center[1]-maskr);
		int uppery=1+lowery+(int)(2.0f*maskr);
		float maskr2=maskr*maskr;
		float masksum=0.0f;
		int maskcount=0;
		if(lowery<0) lowery=0;
		if(lowerx<0) lowerx=0;
		if(upperx>=width) upperx=width-1;
		if(uppery>=height) uppery=height-1;
		for(int i=lowery;i<=uppery;i++){
			for(int j=lowerx;j<=upperx;j++){
				if(((j-center[0])*(j-center[0])+(i-center[1])*(i-center[1]))<=maskr2){
					if(!masked[j+i*width]){
						masksum+=image[j+i*width];
						image[j+i*width]=maskval;
						masked[j+i*width]=true;
						maskcount++;
					}
				}
			}
		}
		return masksum/(float)maskcount;
	}

}
