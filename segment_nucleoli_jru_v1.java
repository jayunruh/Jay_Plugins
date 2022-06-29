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
import ij.plugin.frame.RoiManager;
import jalgs.*;
import jalgs.jseg.*;
import jguis.*;

public class segment_nucleoli_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Nuclear_Channel",1,0);
		gd.addNumericField("Nucleolar_Channel",3,0);
		gd.addNumericField("Nuclear_Thresh_Fraction",0.1,5,15,null);
		gd.addNumericField("Min_Nuclear_Size (pixels)",1000,0);
		gd.addNumericField("Nucleolar_Thresh_Fraction",0.17,5,15,null);
		gd.addNumericField("Min_Nucleolar_Size (pixels)",5,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int nucchan=(int)gd.getNextNumber();
		int nucleolarchan=(int)gd.getNextNumber();
		float nucthresh=(float)gd.getNextNumber();
		int minnucsize=(int)gd.getNextNumber();
		float nucleolarthresh=(float)gd.getNextNumber();
		int minnucleolarsize=(int)gd.getNextNumber();
		float[] nucpix=algutils.convert_arr_float(stack.getPixels(nucchan));
		float[] nucleolarpix=algutils.convert_arr_float(stack.getPixels(nucleolarchan));
		//start by thresholding the nuclei
		//first subtract a 50 pixel rolling ball background
		nucpix=jutils.sub_roll_ball_back(nucpix,50.0f,width,height);
		//blur 1 pixel
		jsmooth.blur2D(nucpix,1.0f,width,height);
		float nucmax=jstatistics.getstatistic("Max",nucpix,null);
		//now find the nuclear objects
		findblobs3 fb=new findblobs3(width,height);
		float[] nucobj=fb.dofindblobs(nucpix,nucthresh*nucmax);
		//eliminate the edges
		fb.clear_edges(nucobj,false);
		//fill holes
		fb.fill_holes(nucobj);
		//filter on size
		fb.filter_area(nucobj,new int[]{minnucsize,10000000},true);
		new ImagePlus("Nuclear Mask",new ByteProcessor(width,height,fb.tobinary(nucobj,false))).show();
		//now do the nucleoli
		nucleolarpix=jutils.sub_roll_ball_back(nucleolarpix,15.0f,width,height);
		//find the maxima in each nucleus
		int[][] lims=fb.getallfilllimits(nucobj);
		float[] nucmaxs=fb.get_all_object_stats(nucobj,nucleolarpix,lims,"Max");
		float[] nucavgs=fb.get_all_object_stats(nucobj,nucleolarpix,lims,"Avg");
		float[] nucstdev=fb.get_all_object_stats(nucobj,nucleolarpix,lims,"StDev");
		int[] nucareas=fb.get_areas(nucobj);
		byte[] nucleolarmask=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(nucobj[i]>0.0f){
				float temp=nucmaxs[(int)(nucobj[i])-1];
				float thresh=nucleolarthresh*temp;
				if(nucleolarpix[i]>=thresh) nucleolarmask[i]=(byte)255;
			}
		}
		//new ImagePlus("Nucleolar Mask",new ByteProcessor(width,height,nucleolarmask)).show();
		int nnuclei=fb.nobjects;
		int[] nuccount=new int[nnuclei];
		float[] nucleolarobj=fb.dofindblobs(nucleolarmask);
		fb.filter_area(nucleolarobj,new int[]{minnucleolarsize,10000000},true);
		int[] clusterids=new int[fb.nobjects];
		for(int i=0;i<nnuclei;i++){
			int[] ids=fb.get_cluster_ids(nucobj,nucleolarobj,i+1);
			nuccount[i]=ids.length;
			for(int j=0;j<ids.length;j++){
				clusterids[ids[j]-1]=i+1;
			}
		}
		//nuclear measurements: 0nucleolar_id, 1id, 2area, 3avg, 4stdev, 5nucleolar count; nucleolar measurements: 6area, 7avg, 8stdev, 9circularity
		int[][] nlims=fb.getallfilllimits(nucleolarobj);
		float[] nucleolaravgs=fb.get_all_object_stats(nucleolarobj,nucleolarpix,nlims,"Avg");
		float[] nucleolarstdev=fb.get_all_object_stats(nucleolarobj,nucleolarpix,nlims,"StDev");
		Polygon[] nucleolaroutlines=fb.get_object_outlines(nucleolarobj);
		float[][] nucleolarmeas=fb.get_area_perim_circ(nucleolarobj,nucleolaroutlines);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		float[][] measurements=new float[fb.nobjects][10];
		for(int i=0;i<fb.nobjects;i++){
			rman.addRoi(new PolygonRoi(nucleolaroutlines[i],Roi.FREEROI));
			int nucid=clusterids[i];
			measurements[i][0]=(float)(i+1);
			measurements[i][1]=(float)nucid;
			measurements[i][2]=(float)nucareas[nucid-1];
			measurements[i][3]=nucavgs[nucid-1];
			measurements[i][4]=nucstdev[nucid-1];
			measurements[i][5]=nuccount[nucid-1];
			measurements[i][6]=nucleolarmeas[0][i];
			measurements[i][7]=nucleolaravgs[i];
			measurements[i][8]=nucleolarstdev[i];
			measurements[i][9]=nucleolarmeas[2][i];
		}
		String[] collabels={"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity"};
		table_tools.create_table("Nucleolar Measurements",measurements,collabels);
	}

}
