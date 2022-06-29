/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.*;
import ij.plugin.filter.GaussianBlur;
import ij.text.*;
import jalgs.*;
import jalgs.jseg.*;
import java.util.Hashtable;
import jguis.*;

public class spb_partition_analysis_jru_v2 implements PlugIn {
	boolean debug;
	ImageStack objstack;

	public void run(String arg) {
		debug=false;
		float threshmult=15.0f;
		int maxnucsize=40;
		int bordersize=30;
		int spbsize=2;
		float nucthresh=0.75f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Debug?",debug);
		gd.addNumericField("SPB_threshold",threshmult,5,15,null);
		gd.addNumericField("Max_Nuc_Size",maxnucsize,0);
		gd.addNumericField("Border_Size",bordersize,0);
		gd.addNumericField("SPB_half_size",spbsize,0);
		gd.addNumericField("Nuclear_threshold",nucthresh,5,15,null);
		gd.addNumericField("SPB_channel",2,0);
		gd.addNumericField("NE_channel",1,0);
		gd.addCheckbox("other_channel",false);
		gd.addNumericField("Other_channel",3,0);
		gd.addCheckbox("rolling_ball_sub?",false);
		gd.addNumericField("rolling_ball_radius",50,0);
		gd.addCheckbox("border_back_sub?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		debug=gd.getNextBoolean();
		threshmult=(float)gd.getNextNumber();
		maxnucsize=(int)gd.getNextNumber();
		bordersize=(int)gd.getNextNumber();
		spbsize=(int)gd.getNextNumber();
		nucthresh=(float)gd.getNextNumber();
		int spbch=(int)gd.getNextNumber()-1;
		int nech=(int)gd.getNextNumber()-1;
		boolean third=gd.getNextBoolean();
		int othch=(int)gd.getNextNumber()-1;
		boolean rollsub=gd.getNextBoolean();
		int rollballrad=(int)gd.getNextNumber();
		boolean bordsub=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int channels=imp.getNChannels();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		float[] maxspbproj=jutils.get3DProjZStat(stack,0,spbch,1,slices,channels,"Max");
		float[] sumspbproj=jutils.get3DProjZStat(stack,0,spbch,1,slices,channels,"Sum");
		float[] sumnpcproj=jutils.get3DProjZStat(stack,0,nech,1,slices,channels,"Sum");
		float[] sumothproj=null;
		if(third) sumothproj=jutils.get3DProjZStat(stack,0,othch,1,slices,channels,"Sum");
		Object[] npcstack=jutils.get3DZSeries(stack,nech,0,1,slices,channels);
		Object[] othstack=null;
		if(third) othstack=jutils.get3DZSeries(stack,othch,0,1,slices,channels);
		Roi roi=imp.getRoi();
		if(roi!=null){
			Polygon polyroi=roi.getPolygon();
			if(polyroi!=null){
				float sumback=jstatistics.getstatistic("Avg",sumnpcproj,width,height,polyroi,null);
				for(int i=0;i<sumnpcproj.length;i++) sumnpcproj[i]-=sumback;
				float sumback2=jstatistics.getstatistic("Avg",sumspbproj,width,height,polyroi,null);
				for(int i=0;i<sumspbproj.length;i++) sumspbproj[i]-=sumback2;
				if(third){
					float sumback3=jstatistics.getstatistic("Avg",sumothproj,width,height,polyroi,null);
					for(int i=0;i<sumothproj.length;i++) sumothproj[i]-=sumback3;
				}
			}
		}
		if(rollsub){
			sumnpcproj=jutils.sub_roll_ball_back(sumnpcproj,rollballrad,width,height);
			sumspbproj=jutils.sub_roll_ball_back(sumspbproj,rollballrad,width,height);
			if(third){
				sumothproj=jutils.sub_roll_ball_back(sumothproj,rollballrad,width,height);
			}
		}
		//now process the spbproj to get the spb positions
		maxspbproj=jsmooth.smooth2D(maxspbproj,width,height);
		//subtract a background
		float[] background=blurImage(maxspbproj,width,height,6.0f);
		for(int i=0;i<width*height;i++) maxspbproj[i]-=background[i];
		float[] tempmsp=maxspbproj.clone();
		//find edges (sobel)
		maxspbproj=(new jsobel(width,height)).do_sobel(maxspbproj)[0];
		//now threshhold at threshmult times the avg intensity
		float avg=jstatistics.getstatistic("Avg",maxspbproj,null);
		byte[] spbs=new byte[width*height];
		float thresh=avg*threshmult;
		for(int i=0;i<width*height;i++) if(maxspbproj[i]>thresh){spbs[i]=(byte)255;}
		//dilate
		(new binary_processing(width,height)).dilate(spbs);
		//find objects
		findblobs3 fb=new findblobs3(width,height);
		float[] objects=fb.dofindblobs(spbs);
		fb.clear_borders(objects,bordersize);
		fb.renumber_objects(objects);
		if(debug){
			new ImagePlus("Spbs",new FloatProcessor(width,height,objects,null)).show();
			objstack=new ImageStack(maxnucsize,maxnucsize);
		}
		//and their coordinates
		float[][] centroids=measure_object.centroids(objects,width,height);
		//now get the areas of the spb's by thresholding at 50% of max-min inside the centroids
		int[] markareas=new int[fb.nobjects];
		int[][] lims=fb.getallfilllimits(objects);
		float[] markmax=fb.get_all_object_stats(objects,tempmsp,lims,"Max");
		float[] markmin=fb.get_all_object_stats(objects,tempmsp,lims,"Min");
		float[][] marksums=new float[fb.nobjects][5];
		float[] markwidths=new float[fb.nobjects];
		for(int i=0;i<objects.length;i++){
			int id=(int)objects[i]-1;
			int y=(int)((float)i/(float)width);
			int x=i-y*width;
			if(id>=0){
				float val=tempmsp[i]-markmin[id];
				if(val>=0.25f*(markmax[id]-markmin[id])){markareas[id]++;}
				marksums[id][0]+=val;
				marksums[id][1]+=val*(float)x;
				marksums[id][2]+=val*(float)y;
				marksums[id][3]+=val*(float)x*(float)x;
				marksums[id][4]+=val*(float)y*(float)y;
			}
		}
		for(int i=0;i<fb.nobjects;i++){
			marksums[i][4]/=marksums[i][0];
			marksums[i][3]/=marksums[i][0];
			marksums[i][2]/=marksums[i][0];
			marksums[i][1]/=marksums[i][0];
			float varx=marksums[i][3]-marksums[i][1]*marksums[i][1];
			float vary=marksums[i][4]-marksums[i][2]*marksums[i][2];
			float stdx=(float)Math.sqrt(varx);
			float stdy=(float)Math.sqrt(vary);
			markwidths[i]=0.5f*2.35f*(stdx+stdy);
			//markwidths[i]=varx;
		}
		//the spb intensities are the sum of the region surrounding the spb centroid in the sum proj image
		float[] spbintensities=new float[centroids.length];
		float[] spbintensitiesspb=new float[centroids.length];
		float[] spbintensitiesoth=new float[centroids.length];
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			rman=new RoiManager();
		}
		rman.runCommand("show all");
		for(int i=0;i<centroids.length;i++){
			centroids[i][0]=(float)Math.round(centroids[i][0]);
			centroids[i][1]=(float)Math.round(centroids[i][1]);
			Rectangle rect=new Rectangle((int)(centroids[i][0]-spbsize),(int)(centroids[i][1]-spbsize),2*spbsize+1,2*spbsize+1);
			spbintensities[i]=jstatistics.getstatistic("Sum",sumnpcproj,width,height,rect,null);
			//might want to subtract the average border intensity from the sbp
			if(bordsub){
				int spbsize2=(int)jstatistics.getstatistic("Count",sumnpcproj,width,height,rect,null);
				Rectangle plusborder=new Rectangle((int)(centroids[i][0]-spbsize-1),(int)(centroids[i][1]-spbsize-1),2*spbsize+3,2*spbsize+3);
				float border=jstatistics.getstatistic("Sum",sumnpcproj,width,height,plusborder,null);
				border-=spbintensities[i];
				int bordercount=(int)jstatistics.getstatistic("Count",sumnpcproj,width,height,plusborder,null);
				bordercount-=spbsize2;
				float borderavg=border/(float)bordercount;
				spbintensities[i]-=borderavg*(float)spbsize2;
			}
			if(third) spbintensitiesoth[i]=jstatistics.getstatistic("Sum",sumothproj,width,height,rect,null);
			spbintensitiesspb[i]=jstatistics.getstatistic("Sum",sumspbproj,width,height,rect,null);
			rman.addRoi(new PointRoi((int)centroids[i][0],(int)centroids[i][1]));
		}
		//now clear the spb to get the nuclear mask
		float[] nucmask=sumnpcproj.clone();
		elim_points(nucmask,width,height,centroids);
		nucmask=blurImage(nucmask,width,height,3.0f);
		//now get the sum of each nucleus
		TextWindow tw=null;
		if(!third) tw=new TextWindow("Intensities","id\tspbint\tnucint\tratio\tspbmarkint\tnucarea\tnucz\tspbarea\tsbpwidth","",400,400);
		else tw=new TextWindow("Intensities","id\tspbint\tnucint\tratio\totherint\totherratio\tspbmarkint\tnucarea\tnucz\tspbarea\tspbwidth","",400,400);
		float[] nucint=new float[centroids.length];
		float[] nucarea=new float[centroids.length];
		float[] nucz=new float[centroids.length];
		float[] othint=new float[centroids.length];
		for(int i=0;i<centroids.length;i++){
			Polygon outline=get_local_object(nucmask,(int)centroids[i][0],(int)centroids[i][1],width,height,maxnucsize,nucthresh);
			nucint[i]=jstatistics.getstatistic("Sum",sumnpcproj,width,height,outline,null);
			if(third) othint[i]=jstatistics.getstatistic("Sum",sumothproj,width,height,outline,null);
			//find the maximum z plane for the nucleus
			nucz[i]=get_max_slice(npcstack,outline,width,height);
			nucarea[i]=measure_object.area(outline);
			if(!third) tw.append(""+(i+1)+"\t"+spbintensities[i]+"\t"+nucint[i]+"\t"+(spbintensities[i]/nucint[i])+"\t"+spbintensitiesspb[i]+"\t"+nucarea[i]+"\t"+nucz[i]+"\t"+markareas[i]+"\t"+markwidths[i]+"\n");
			else tw.append(""+(i+1)+"\t"+spbintensities[i]+"\t"+nucint[i]+"\t"+(spbintensities[i]/nucint[i])+"\t"+othint[i]+"\t"+(spbintensitiesoth[i]/othint[i])+"\t"+spbintensitiesspb[i]+"\t"+nucarea[i]+"\t"+nucz[i]+"\t"+markareas[i]+"\t"+markwidths[i]+"\n");
		}
		if(debug) new ImagePlus("Masks",objstack).show();
	}

	public int get_max_slice(Object[] stack,Polygon outline,int width,int height){
		int maxslice=0;
		float maxval=jstatistics.getstatistic("Sum",stack[0],width,height,outline,null);
		for(int i=1;i<stack.length;i++){
			float val=jstatistics.getstatistic("Sum",stack[i],width,height,outline,null);
			if(val>maxval){
				maxval=val;
				maxslice=i;
			}
		}
		return maxslice;
	}		

	public void elim_points(float[] image,int width,int height,float[][] coords){
		for(int i=0;i<coords.length;i++){
			elim_point(image,width,height,(int)coords[i][0],(int)coords[i][1]);
		}
	}

	public void elim_point(float[] image,int width,int height,int x,int y){
		float perimavg=0.0f;
		int perimsize=16;
		if(x<2) return;
		if(x>(width-3)) return;
		if(y<2) return;
		if(y>(height-3)) return;
		//get the 3x3 perimeter average
		for(int i=(x-2);i<=(x+2);i++){
			perimavg+=image[i+(y-2)*width];
			perimavg+=image[i+(y+2)*width];
		}
		for(int i=(y-1);i<=(y+1);i++){
			perimavg+=image[x-2+i*width];
			perimavg+=image[x+2+i*width];
		}
		perimavg/=(float)perimsize;
		//fill the 3x3 area with the perimeter avg
		for(int i=(x-1);i<=(x+1);i++){
			for(int j=(y-1);j<=(y+1);j++){
				image[i+j*width]=perimavg;
			}
		}
	}

	public Polygon get_local_object(float[] image,int x,int y,int width,int height,int maxsize,float fraction){
		int startx=x-maxsize/2; if(startx<0) startx=0;
		int endx=startx+maxsize; if(endx>=width) endx=width-1;
		int starty=y-maxsize/2; if(starty<0) starty=0;
		int endy=starty+maxsize; if(endy>=height) endy=height-1;
		int newwidth=endx-startx+1;
		int newheight=endy-starty+1;
		byte[] subimage=new byte[newwidth*newheight];
		float thresh=fraction*image[x+y*width];
		int counter=0;
		for(int i=starty;i<=endy;i++){
			for(int j=startx;j<=endx;j++){
				if(image[j+i*width]>thresh) subimage[counter]=(byte)255;
				counter++;
			}
		}
		findblobs3 fb=new findblobs3(newwidth,newheight);
		float[] blobs=fb.dofindblobs(subimage);
		if(debug){
			if(newwidth==maxsize && newheight==maxsize){
				objstack.addSlice("",blobs);
			} else {
				float[] padded=algutils.pad_2D(blobs,newwidth,newheight,maxsize,maxsize,0);
				objstack.addSlice("",padded);
			}
		}
		int newx=x-startx;
		int newy=y-starty;
		//int id=(int)blobs[newx+newy*newwidth];
		//IJ.log("preoutline "+id);
		//new ImagePlus("blob",new FloatProcessor(newwidth,newheight,blobs,null)).show();
		Polygon outline=fb.get_object_outline(blobs,(int)blobs[newx+newy*newwidth]);
		//IJ.log("postoutline");
		//now we have to shift the polygon by startx and starty
		outline.translate(startx,starty);
		return outline;
	}

	public float[] blurImage(float[] image,int width,int height,float sigma){
		GaussianBlur gb=new GaussianBlur();
		FloatProcessor fp=new FloatProcessor(width,height,image.clone(),null);
		gb.blurFloat(fp,(double)sigma,(double)sigma,0.0002);
		gb=null;
		return (float[])fp.getPixels();
	}
		
}
