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
import java.awt.Frame;
import java.awt.Polygon;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.RoiManager;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class segment_nucleoli_jru_v3 implements PlugIn {
	//this version is back to lower throughput for easier troubleshooting
	//also tried to streamline some of the background subtraction for larger scale images

	public void run(String arg) {
		float nucblurstdev=12.0f; float nucrbr=200.0f; float nucdbs=200.0f; float nucthresh=0.15f; int nucmina=5000; int nucmaxa=18000; float nucmincirc=0.8f;
		int nuclmina=4; float nuclblurstdev=2.0f; float nuclrbr=25.0f; float nuclmeasrbr=200.0f; float nucldbs=30.0f; float nuclthresh=0.2f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Nuclear_Channel",3,0);
		gd.addNumericField("Nucleolar_Channel",2,0);
		gd.addNumericField("Nucleolar_Meas_Channel",1,0);
		gd.addNumericField("Third_Meas_Channel",-1,0);
		gd.addNumericField("Nuclear_Blur_Stdev",nucblurstdev,5,15,null);
		gd.addNumericField("Nuclear_Roll_Ball_Rad",nucrbr,5,15,null);
		gd.addNumericField("Nuclear_Div_By_Stdev",nucdbs,5,15,null);
		gd.addNumericField("Nuclear_Thresh",nucthresh,5,15,null);
		gd.addNumericField("Nuclear_Min_Area",nucmina,0);
		gd.addNumericField("Nuclear_Max_Area",nucmaxa,0);
		gd.addNumericField("Nuclear_Min_Circularity",nucmincirc,5,15,null);
		gd.addNumericField("Nucleolar_Min_Area",nuclmina,0);
		gd.addNumericField("Nucleolar_Blur_Stdev",nuclblurstdev,5,15,null);
		gd.addNumericField("Nucleolar_Roll_Ball_Rad",nuclrbr,5,15,null);
		gd.addNumericField("Nucleolar_Meas_Roll_Ball_Rad",nuclmeasrbr,5,15,null);
		gd.addNumericField("Nucleolar_Div_By_Stdev",nucldbs,5,15,null);
		gd.addNumericField("Nucleolar_Thresh",nuclthresh,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int nucchan=(int)gd.getNextNumber();
		int nuclchan=(int)gd.getNextNumber();
		int measchan=(int)gd.getNextNumber();
		int measchan3=(int)gd.getNextNumber();
		nucblurstdev=(float)gd.getNextNumber();
		nucrbr=(float)gd.getNextNumber();
		nucdbs=(float)gd.getNextNumber();
		nucthresh=(float)gd.getNextNumber();
		nucmina=(int)gd.getNextNumber();
		nucmaxa=(int)gd.getNextNumber();
		nucmincirc=(float)gd.getNextNumber();
		nuclmina=(int)gd.getNextNumber();
		nuclblurstdev=(float)gd.getNextNumber();
		nuclrbr=(float)gd.getNextNumber();
		nuclmeasrbr=(float)gd.getNextNumber();
		nucldbs=(float)gd.getNextNumber();
		nuclthresh=(float)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();

		//start by getting the channel images
		float[] nuc=algutils.convert_arr_float2(stack.getPixels(nucchan));
		float[] nucl=algutils.convert_arr_float2(stack.getPixels(nuclchan));
		float[][] meas=new float[][]{nucl};
		if(measchan!=nuclchan) meas=new float[][]{algutils.convert_arr_float2(stack.getPixels(measchan))};
		float[] third=null;
		if(measchan3>0){
			third=algutils.convert_arr_float2(stack.getPixels(measchan3));
			float[] tmeas=meas[0];
			meas=new float[][]{tmeas,third};
		}

		//segment the nuclei: gaus blur (try 8), roll ball sub (try 200), div by blurred (try 200), threshold (try 0.1), filter size (try 5000 to 18000)

		float[] nucobj=segmentNuclei(nuc,width,height,nucblurstdev,nucrbr,nucdbs,nucthresh,nucmina,nucmaxa,nucmincirc);
		new ImagePlus("nuc objects",new FloatProcessor(width,height,nucobj,null)).show();

		//segment the nucleoli

		Object[] temp=segmentNucleoli(nucl,meas,nucobj,width,height,nuclblurstdev,nuclrbr,nuclmeasrbr,nucldbs,nuclthresh,nuclmina);
		//this has measurements, objects, and outlines--output them
		float[] nuclobj=(float[])temp[1];
		//new ImagePlus("nucleolar objects",new FloatProcessor(width,height,nuclobj,null)).show();
		Polygon[] nucloutlines=(Polygon[])temp[2];
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		for(int i=0;i<nucloutlines.length;i++){
			rman.addRoi(new PolygonRoi(nucloutlines[i],Roi.POLYGON));
		}
		String[] collabels={"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","nuclear_circularity","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity"};
		if(measchan3>0){
			collabels=new String[]{"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","nuclear_circularity","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity","thirdnucavg","thirdnucstdev","thirdnuclavg","thirdnuclstdev"};
		}
		float[][] measvals=(float[][])temp[0];
		table_tools.create_table("Nucleolar_Measurements",measvals,collabels);
		List<List<String>> meastable=table_tools.table2listtable(measvals);
		//now summarize the nuclear parameters
		//need to generate nuclear integral, nucleolar integral, nucleoplasmic integral, f nucleolar, tot nucleolar area
		//firstly sort the table by nuclear_id
		table_tools.sort_listtable(meastable,1);
		//now add a nucleolar sum column
		for(int i=0;i<meastable.size();i++){
			float nuclarea=table_tools.get_number(meastable,i,7);
			float nuclavg=table_tools.get_number(meastable,i,8);
			meastable.get(i).add(""+nuclarea*nuclavg);
			if(measchan3>0){
				float thirdavg=table_tools.get_number(meastable,i,13);
				meastable.get(i).add(""+nuclarea*thirdavg);
			}
		}
		//now get the sums
		List<List<String>> sumtable=table_tools.get_cell_stat_list(meastable,1,"Sum",true);
		//and the avgs
		List<List<String>> avgtable=table_tools.get_cell_stat_list(meastable,1,"Avg",true);
		//for the summary table, use the average nuclear parameters (first 7) and the sum nucleolar parameters (except the stdev and circularity)
		String[] sumcollabels={"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","nuclear_circularity","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity","nucleolar_sum","nuclear_sum"};
		if(measchan3>0){
			sumcollabels=new String[]{"id","nuclear_id","nuclear_area","nuclear_avg","nuclear_stdev","nuclear_circularity","number_nucleoli","nucleolar_area","nucleolar_avg","nucleolar_stdev","nucleolar_circularity","thirdnucavg","thirdnucstdev","thirdnuclavg","thirdnuclstdev","nucleolar_sum","thirdnuclsum","nuclear_sum","thirdnucsum"};
		}
		for(int i=0;i<avgtable.size();i++){
			if(measchan3>0){
				List<String> avgrow=avgtable.get(i);
				List<String> sumrow=sumtable.get(i);
				avgrow.set(7,sumrow.get(7)); //the total nucleolar area
				avgrow.set(15,sumrow.get(15)); //the total nucleolar sum
				float nuclsum=Float.parseFloat(sumrow.get(15)); //the total nucleolar sum
				float nuclarea=Float.parseFloat(sumrow.get(7)); //the total nucleolar area
				avgrow.set(8,""+nuclsum/nuclarea); //the revised nucleolar avg
				avgrow.set(16,sumrow.get(16));
				float thirdnuclsum=Float.parseFloat(sumrow.get(16));
				avgrow.set(13,sumrow.get(13));
				float nucavg=Float.parseFloat(avgrow.get(3)); //the nuclear avg
				float nucarea=Float.parseFloat(avgrow.get(2)); //the nuclear area
				avgrow.add(""+nucarea*nucavg);
				float thirdnucavg=Float.parseFloat(avgrow.get(11));
				avgrow.add(""+nucarea*thirdnucavg);
			} else {
				List<String> avgrow=avgtable.get(i);
				List<String> sumrow=sumtable.get(i);
				avgrow.set(7,sumrow.get(7)); //the total nucleolar area
				avgrow.set(11,sumrow.get(11)); //the total nucleolar sum
				float nuclsum=Float.parseFloat(sumrow.get(11)); //the total nucleolar sum
				float nuclarea=Float.parseFloat(sumrow.get(7)); //the total nucleolar area
				avgrow.set(8,""+nuclsum/nuclarea); //the revised nucleolar avg
				float nucavg=Float.parseFloat(avgrow.get(3)); //the nuclear avg
				float nucarea=Float.parseFloat(avgrow.get(2)); //the nuclear area
				avgrow.add(""+nucarea*nucavg);
			}
		}
		table_tools.create_table("Nuclear_Summaries",avgtable,sumcollabels);
	}

	public static Object[] segmentNucleoli(float[] nucleolipix,float[][] nucmeaspix,float[] nucobj,int width,int height,float blurstdev,float rollballrad,float measrollballrad,float divblurstdev,float threshfrac,int nuclmina){
		//the workflow here is blur (try 2), roll ball sub (try 25), div by blurred (try 30), thresh within nuclei (try 0.15)
		float[] backsub=nucleolipix.clone();
		//a minimal blur to decrease noise
		jsmooth.blur2D(backsub,blurstdev,width,height);
		//now a rolling ball subtraction
		backsub=jutils.sub_roll_ball_back(backsub,rollballrad,width,height);
		float[][] meas=new float[][]{backsub};
		if(nucmeaspix!=null){
			meas=new float[nucmeaspix.length][];
			for(int i=0;i<nucmeaspix.length;i++){
				meas[i]=nucmeaspix[i].clone();
				jsmooth.blur2D(meas[i],blurstdev,width,height);
				meas[i]=jutils.sub_roll_ball_back(meas[i],measrollballrad,width,height);
			}
		}
		findblobs3 fb=new findblobs3(width,height);
		fb.set_objects(nucobj);
		//now get the maxima within each nucleus
		int[][] filllims=fb.getallfilllimits(nucobj);
		float[] nucmaxs=fb.get_all_object_stats(nucobj,backsub,filllims,"Max");
		float[] nucavgs=fb.get_all_object_stats(nucobj,backsub,filllims,"Avg");
		float[] nucstdevs=fb.get_all_object_stats(nucobj,backsub,filllims,"StDev");
		float[][] nucmeasavgs=new float[][]{nucavgs};
		float[][] nucmeasstdevs=new float[][]{nucstdevs};
		if(nucmeaspix!=null){
			nucmeasavgs=new float[nucmeaspix.length][];
			nucmeasstdevs=new float[nucmeaspix.length][];
			for(int i=0;i<nucmeaspix.length;i++){
				nucmeasavgs[i]=fb.get_all_object_stats(nucobj,meas[i],filllims,"Avg");
				nucmeasstdevs[i]=fb.get_all_object_stats(nucobj,meas[i],filllims,"StDev");
			}
		}
		Polygon[] nuclearoutlines=fb.get_object_outlines(nucobj);
		float[][] nuclearmeas=fb.get_area_perim_circ(nucobj,nuclearoutlines);
		int[] nucareas=fb.get_areas(nucobj);
		//now segment within each nucleus based on the nucleolar maximum for that nucleus
		byte[] nucleolarmask=new byte[width*height];
		for(int i=0;i<width*height;i++){
			if(nucobj[i]>0.0f){
				float temp=nucmaxs[(int)(nucobj[i])-1];
				float thresh2=threshfrac*temp;
				if(backsub[i]>=thresh2) nucleolarmask[i]=(byte)255;
			} else {
				backsub[i]=0.0f;
			}
		}
		//new ImagePlus("blurred nucleoli",new FloatProcessor(width,height,backsub,null)).show();
		int nnuclei=fb.nobjects;
		int[] nuccount=new int[nnuclei];
		float[] nucleolarobj=fb.dofindblobs(nucleolarmask);
		fb.filter_area(nucleolarobj,new int[]{nuclmina,10000000},true);
		int[] clusterids=new int[fb.nobjects];
		for(int i=0;i<nnuclei;i++){
			int[] ids=fb.get_cluster_ids(nucobj,nucleolarobj,i+1);
			nuccount[i]=ids.length;
			for(int j=0;j<ids.length;j++){
				clusterids[ids[j]-1]=i+1;
			}
		}
		//nuclear measurements: 0nucleolar_id, 1nucid, 2area, 3avg, 4stdev, 5nuccirc, 6nucleolar count; nucleolar measurements: 7area, 8avg, 9stdev, 10circularity
		int[][] nlims=fb.getallfilllimits(nucleolarobj);
		float[][] nucleolaravgs=new float[meas.length][];
		float[][] nucleolarstdev=new float[meas.length][];
		for(int i=0;i<meas.length;i++){
			nucleolaravgs[i]=fb.get_all_object_stats(nucleolarobj,meas[i],nlims,"Avg");
			nucleolarstdev[i]=fb.get_all_object_stats(nucleolarobj,meas[i],nlims,"StDev");
		}
		Polygon[] nucleolaroutlines=fb.get_object_outlines(nucleolarobj);
		float[][] nucleolarmeas=fb.get_area_perim_circ(nucleolarobj,nucleolaroutlines);
		float[][] measurements=new float[fb.nobjects][11+(meas.length-1)*4];
		for(int i=0;i<fb.nobjects;i++){
			int nucid=clusterids[i];
			measurements[i][0]=(float)(i+1);
			measurements[i][1]=(float)nucid;
			measurements[i][2]=(float)nucareas[nucid-1];
			measurements[i][3]=nucmeasavgs[0][nucid-1];
			measurements[i][4]=nucmeasstdevs[0][nucid-1];
			measurements[i][5]=nuclearmeas[2][nucid-1];
			measurements[i][6]=nuccount[nucid-1];
			measurements[i][7]=nucleolarmeas[0][i];
			measurements[i][8]=nucleolaravgs[0][i];
			measurements[i][9]=nucleolarstdev[0][i];
			measurements[i][10]=nucleolarmeas[2][i];
			for(int j=0;j<(meas.length-1);j++){
				measurements[i][j*4+11]=nucmeasavgs[j][nucid-1];
				measurements[i][j*4+12]=nucmeasstdevs[j][nucid-1];
				measurements[i][j*4+13]=nucleolaravgs[j][nucid-1];
				measurements[i][j*4+14]=nucleolarstdev[j][nucid-1];
			}
		}
		return new Object[]{measurements,nucleolarobj,nucleolaroutlines};
	}

	public static float[] segmentNuclei(float[] nucpix,int width,int height,float blurstdev,float rollballrad,float divblurstdev, float threshfrac,int minarea,int maxarea,float mincirc){
		//here is the workflow: gaus blur (try 8), roll ball sub (try 200), div by blurred (try 200), threshold (try 0.15), filter size (try 10000 to 30000), filter circularity (try 0.8)
		float[] backsub=nucpix.clone();
		//a blur to even out the nuclei
		//jsmooth.blur2D(backsub,blurstdev,width,height);
		if(blurstdev>0.0f) backsub=jutils.gaussian_blur(backsub,blurstdev,width,height);
		//do a rolling ball subtraction
		if(rollballrad>0.0f) backsub=jutils.sub_roll_ball_back(backsub,rollballrad,width,height);
		//divide by blurred to make the nuclei all the same intensity
		if(divblurstdev>0.0f) div_by_smooth(backsub,width,height,divblurstdev,1.0f,1.0f);
		findblobs3 fb=new findblobs3(width,height);
		float max=jstatistics.getstatistic("Max",backsub,null);
		//float[] percentile={99.9f};
		//float max=jstatistics.getstatistic("percentile",backsub,percentile);
		//max=percentile[0];
		//IJ.log(""+max);
		//threshold at thresh*max
		float[] objects=fb.dofindblobs(backsub,threshfrac*max);
		//finally filter the objects and clear the edges
		fb.clear_edges(objects,true);
		fb.filter_area_circ(objects,new float[]{(float)minarea,(float)maxarea,mincirc,1.0f});
		//fb.filter_area(objects,new int[]{minarea,maxarea},true);
		fb.renumber_objects(objects);
		fb.fill_holes(objects);
		return objects;
	}

	public static void div_by_smooth(float[] image,int width,int height,float blurstdev,float mindiv,float multiplier){
		//float[] smoothed=image.clone();
		//jsmooth.blur2D(smoothed,blurstdev,width,height);
		float[] smoothed=jutils.gaussian_blur(image,blurstdev,width,height);
		for(int i=0;i<width*height;i++){
			float temp=smoothed[i];
			if(temp<mindiv) temp=mindiv;
			image[i]/=(temp/multiplier);
		}
		smoothed=null;
	}

}
