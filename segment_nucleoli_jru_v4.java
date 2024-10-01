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
import ij.io.*;
import java.io.*;
import java.util.concurrent.*;

public class segment_nucleoli_jru_v4 implements PlugIn {
	//this version utilizes nuclear minimum and measures in a third channel

	public void run(String arg) {
		/*ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		//assume that this is a two channel image with dapi first
		ImageStack stack=imp.getStack();
		float[] nucpix=algutils.convert_arr_float(stack.getPixels(1));
		float[] nucleolipix=algutils.convert_arr_float(stack.getPixels(2));
		float[] nucobj=segmentNuclei(nucpix,width,height,0.1f,1000,4000);
		Object[] measurements=segmentNucleoli(nucleolipix,nucobj,width,height,0.1f,4);
		String[] collabels={"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity"};
		table_tools.create_table("Nucleolar Measurements",(float[][])measurements[0],collabels);
		byte[] nucmask=findblobs3.threshimage(nucobj,0.5f);
		new ImagePlus("Nuclei",new ByteProcessor(width,height,nucmask)).show();
		Polygon[] polys=(Polygon[])measurements[2];
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		for(int i=0;i<polys.length;i++){
			rman.addRoi(new PolygonRoi(polys[i],Roi.FREEROI));
		}*/
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Nuclear threshold",0.1,5,15,null);
		//gd.addNumericField("Nucleolar threshold",0.15,5,15,null);
		gd.addNumericField("Nucleolar threshold",0.4,5,15,null);
		gd.addNumericField("Scale factor",1,0);
		gd.addNumericField("N_Threads",5,0);
		gd.addNumericField("Min_Nuc_Size",1000,0);
		gd.addNumericField("Max_Nuc_Size",4000,0);
		gd.addNumericField("Min_Nucleolar_Size",4,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float nucthresh=(float)gd.getNextNumber();
		float nucleolithresh=(float)gd.getNextNumber();
		int scale=(int)gd.getNextNumber();
		int nthreads=(int)gd.getNextNumber();
		int nucmin=(int)gd.getNextNumber();
		int nucmax=(int)gd.getNextNumber();
		int nuclmin=(int)gd.getNextNumber();
		DirectoryChooser dc=new DirectoryChooser("Choose Directory");
		String dir=dc.getDirectory();
		if(dir==null) return;
		String[] list=new File(dir).list();
		jsort.javasort_order(list);
		ExecutorService executor=Executors.newFixedThreadPool(nthreads);
		for(int i=0;i<list.length;i++){
			if(list[i].endsWith("1.tif")){//the dapi image
				String nucleoliname=list[i].substring(0,list[i].length()-5)+"3.tif"; //the nucleolar image
				String thirdname=list[i].substring(0,list[i].length()-5)+"2.tif"; //the third image
				Runnable worker=new measure_nucleoli2(dir+list[i],dir+nucleoliname,dir+thirdname,nucthresh,nucleolithresh,scale,nucmin,nucmax,nuclmin);
				executor.execute(worker);
				//IJ.log(nucleoliname);
			}
		}
	}

}

class measure_nucleoli2 implements Runnable{
	String dapiname,nucleoliname,thirdname;
	float nucthresh,nucleolithresh;
	int scale,nucmin,nucmax,nuclmin;

	public measure_nucleoli2(String dapiname,String nucleoliname,String thirdname,float nucthresh,float nucleolithresh,int scale,int nucmin,int nucmax,int nuclmin){
		this.dapiname=dapiname;
		this.nucleoliname=nucleoliname;
		this.thirdname=thirdname;
		this.nucthresh=nucthresh;
		this.nucleolithresh=nucleolithresh;
		this.scale=scale;
		this.nucmin=nucmin;
		this.nucmax=nucmax;
		this.nuclmin=nuclmin;
	}

	//public static void measure_images(String dapiname,String nucleoliname,float nucthresh,float nucleolithresh){
	public void run(){
		//first open the images
		ImagePlus dapiimp=IJ.openImage(dapiname);
		int width=dapiimp.getWidth(); int height=dapiimp.getHeight();
		ImagePlus nucleoliimp=IJ.openImage(nucleoliname);
		ImagePlus thirdimp=IJ.openImage(thirdname);
		float[] nucpix=(float[])dapiimp.getProcessor().convertToFloat().getPixels();
		float[] nucleolipix=(float[])nucleoliimp.getProcessor().convertToFloat().getPixels();
		float[] thirdpix=(float[])thirdimp.getProcessor().convertToFloat().getPixels();
		float[] nucobj=segmentNuclei(nucpix,width,height,nucthresh,nucmin*scale,nucmax*scale,scale);
		Object[] measurements=segmentNucleoli(nucleolipix,nucobj,width,height,nucleolithresh,nuclmin*scale,scale,thirdpix);
		String[] collabels={"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity","third_nucavg","third_nucstdev","third_nuclavg","third_nuclstdev"};
		String maskname=dapiname.substring(0,dapiname.length()-4)+"_mask.tif";
		String roiname=dapiname.substring(0,dapiname.length()-4)+".zip";
		String tabname=dapiname.substring(0,dapiname.length()-4)+".xls";
		byte[] nucmask=findblobs3.threshimage(nucobj,0.5f);
		ImagePlus maskimp=new ImagePlus("Nuclei",new ByteProcessor(width,height,nucmask));
		//now save the mask, rois, and table
		IJ.save(maskimp,maskname);
		Polygon[] polys=(Polygon[])measurements[2];
		if(polys.length>0) multi_roi_writer.writeRois(polys,roiname);
		float[][] measvals=(float[][])measurements[0];
		if(measvals.length>0) table_tools.writeTableToFile(tabname,table_tools.print_string_array(collabels),table_tools.print_float_array(measvals));
		System.out.println(nucleoliname);
	}

	public static Object[] segmentNucleoli(float[] nucleolipix,float[] nucobj,int width,int height,float thresh,int minarea,int scale,float[] thirdpix){
		//start with a rolling ball subtraction
		float[] backsub=jutils.sub_roll_ball_back(nucleolipix,15*scale,width,height);
		float[] thirdbacksub=jutils.sub_roll_ball_back(thirdpix,15*scale,width,height);
		//a minimal blur to decrease noise
		jsmooth.blur2D(backsub,0.7f*(float)scale,width,height);
		jsmooth.blur2D(thirdbacksub,0.7f*(float)scale,width,height);
		findblobs3 fb=new findblobs3(width,height);
		fb.set_objects(nucobj);
		//now get the maxima within each nucleus
		int[][] filllims=fb.getallfilllimits(nucobj);
		float[] nucmaxs=fb.get_all_object_stats(nucobj,backsub,filllims,"Max");
		float[] nucmins=fb.get_all_object_stats(nucobj,backsub,filllims,"Min");
		float[] nucavgs=fb.get_all_object_stats(nucobj,backsub,filllims,"Avg");
		float[] nucstdevs=fb.get_all_object_stats(nucobj,backsub,filllims,"StDev");
		float[] thirdnucavgs=fb.get_all_object_stats(nucobj,thirdbacksub,filllims,"Avg");
		float[] thirdnucstdevs=fb.get_all_object_stats(nucobj,thirdbacksub,filllims,"StDev");
		int[] nucareas=fb.get_areas(nucobj);
		//now segment within each nucleus based on the nucleolar maximum for that nucleus
		byte[] nucleolarmask=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(nucobj[i]>0.0f){
				//float temp=nucmaxs[(int)(nucobj[i])-1];
				//float thresh2=thresh*temp;
				float temp=nucmaxs[(int)(nucobj[i])-1]-nucmins[(int)(nucobj[i])-1];
				float thresh2=thresh*temp+nucmins[(int)(nucobj[i])-1];
				if(backsub[i]>=thresh2) nucleolarmask[i]=(byte)255;
			}
		}
		int nnuclei=fb.nobjects;
		int[] nuccount=new int[nnuclei];
		float[] nucleolarobj=fb.dofindblobs(nucleolarmask);
		fb.filter_area(nucleolarobj,new int[]{minarea,10000000},true);
		int[] clusterids=new int[fb.nobjects];
		for(int i=0;i<nnuclei;i++){
			int[] ids=fb.get_cluster_ids(nucobj,nucleolarobj,i+1);
			nuccount[i]=ids.length;
			for(int j=0;j<ids.length;j++){
				clusterids[ids[j]-1]=i+1;
			}
		}
		//nuclear measurements: 0nucleolar_id, 1nucid, 2area, 3avg, 4stdev, 5nucleolar count; nucleolar measurements: 6area, 7avg, 8stdev, 9circularity, 10thirdnucavg, 11thirdnucstdev, 12thirdavg, 13thirdstdev
		int[][] nlims=fb.getallfilllimits(nucleolarobj);
		float[] nucleolaravgs=fb.get_all_object_stats(nucleolarobj,backsub,nlims,"Avg");
		float[] nucleolarstdev=fb.get_all_object_stats(nucleolarobj,backsub,nlims,"StDev");
		float[] thirdnucleolaravgs=fb.get_all_object_stats(nucleolarobj,thirdbacksub,nlims,"Avg");
		float[] thirdnucleolarstdev=fb.get_all_object_stats(nucleolarobj,thirdbacksub,nlims,"StDev");
		Polygon[] nucleolaroutlines=fb.get_object_outlines(nucleolarobj);
		float[][] nucleolarmeas=fb.get_area_perim_circ(nucleolarobj,nucleolaroutlines);
		float[][] measurements=new float[fb.nobjects][14];
		for(int i=0;i<fb.nobjects;i++){
			int nucid=clusterids[i];
			measurements[i][0]=(float)(i+1);
			measurements[i][1]=(float)nucid;
			measurements[i][2]=(float)nucareas[nucid-1];
			measurements[i][3]=nucavgs[nucid-1];
			measurements[i][4]=nucstdevs[nucid-1];
			measurements[i][5]=nuccount[nucid-1];
			measurements[i][6]=nucleolarmeas[0][i];
			measurements[i][7]=nucleolaravgs[i];
			measurements[i][8]=nucleolarstdev[i];
			measurements[i][9]=nucleolarmeas[2][i];
			measurements[i][10]=thirdnucavgs[nucid-1];
			measurements[i][11]=thirdnucstdevs[nucid-1];
			measurements[i][12]=thirdnucleolaravgs[i];
			measurements[i][13]=thirdnucleolarstdev[i];
		}
		return new Object[]{measurements,nucleolarobj,nucleolaroutlines};
	}

	public static float[] segmentNuclei(float[] nucpix,int width,int height,float thresh,int minarea,int maxarea,int scale){
		//do a rolling ball subtraction
		float[] backsub=jutils.sub_roll_ball_back(nucpix,50*scale,width,height);
		//a blur to even out the nuclei
		jsmooth.blur2D(backsub,2.0f*(float)scale,width,height);
		//divide by blurred to make the nuclei all the same intensity
		div_by_smooth(backsub,width,height,50.0f*(float)scale,1.0f,1.0f);
		findblobs3 fb=new findblobs3(width,height);
		float max=jstatistics.getstatistic("Max",backsub,null);
		//threshold at thresh*max
		float[] objects=fb.dofindblobs(backsub,thresh*max);
		//finally filter the objects and clear the edges
		fb.clear_edges(objects,true);
		fb.filter_area(objects,new int[]{minarea,maxarea},true);
		fb.fill_holes(objects);
		return objects;
	}

	public static void div_by_smooth(float[] image,int width,int height,float blurstdev,float mindiv,float multiplier){
		float[] smoothed=image.clone();
		jsmooth.blur2D(smoothed,blurstdev,width,height);
		for(int i=0;i<width*height;i++){
			float temp=smoothed[i];
			if(temp<mindiv) temp=mindiv;
			image[i]/=(temp/multiplier);
		}
		smoothed=null;
	}

}
