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
import java.awt.Frame;
import java.awt.Polygon;
import ij.plugin.*;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;
import ij.text.*;
import ij.io.*;
import java.util.*;

public class flowsight_fret_map_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open CIF Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		DirectoryChooser dc=new DirectoryChooser("Save Directory");
		String outdir=dc.getDirectory();
		if(outdir==null) return;
		LOCI_series_reader lsr=new LOCI_series_reader(directory,fname,false);
		ImagePlus imp1=(ImagePlus)lsr.getNextFrame();
		int nch=imp1.getStack().getSize();
		ImagePlus imp2=(ImagePlus)lsr.getNextFrame();
		//imp1.show();
		//imp2.show();
		IJ.log(imp1.getTitle());
		IJ.log(imp2.getTitle());
		String labels="object\tarea\ttransavg\tgradcv\tgradavg\tgradmax\tx\ty\tangle\tmajor\taspectratio\tcirc";
		for(int i=0;i<nch;i++){
			labels+="\tch"+(i+1);
		}
		labels+="\tfretratio";
		float[][] limits=null;
		StringBuffer sb=new StringBuffer();

		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Trans_channel (1-6)",3,0);
		gd.addNumericField("Mask_channel (1-6)",3,0);
		gd.addNumericField("numerator_channel (1-6)",2,0);
		gd.addNumericField("denominator_channel (1-6)",5,0);
		gd.addNumericField("smoothing stdev",2.0,5,15,null);
		gd.addNumericField("threshold (denominator, times average)",0.75f,5,15,null);
		gd.addNumericField("number of images to show",10,0);
		gd.addNumericField("min_ratio_to_show",0.0,5,15,null);
		gd.addNumericField("max_ratio_to_show",0.35,5,15,null);
		gd.addCheckbox("Output_Raw?",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int transch=(int)gd.getNextNumber()-1;
		int maskch=(int)gd.getNextNumber()-1;
		int numch=(int)gd.getNextNumber()-1;
		int dench=(int)gd.getNextNumber()-1;
		float smoothstdev=(float)gd.getNextNumber();
		float threshmult=(float)gd.getNextNumber();
		int numshow=(int)gd.getNextNumber();
		float minratio=(float)gd.getNextNumber();
		float maxratio=(float)gd.getNextNumber();
		boolean raw=gd.getNextBoolean();
		Object[] temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,threshmult,smoothstdev);
		int counter=0;
		sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
		List<ImagePlus> selected=new ArrayList<ImagePlus>();
		//if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
		ImagePlus tempimp=new ImagePlus("HeatMap1.tif",processIP((FloatProcessor)temp[1],minratio,maxratio,raw));
		if(numshow>0) tempimp.show();
		saveImp(tempimp,outdir,numshow>0);
		for(int i=1;i<lsr.nseries/2;i++){
			imp1=(ImagePlus)lsr.getNextFrame();
			imp2=(ImagePlus)lsr.getNextFrame();
			IJ.log(imp1.getTitle());
			IJ.log(imp2.getTitle());
			temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,threshmult,smoothstdev);
			counter++;
			sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
			//if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
			ImagePlus tempimp2=new ImagePlus("HeatMap"+(i+1)+".tif",processIP((FloatProcessor)temp[1],minratio,maxratio,raw));
			if(i<=numshow) tempimp2.show();
			saveImp(tempimp2,outdir,i<=numshow);
			IJ.showProgress(i,lsr.nseries/2);
			if(IJ.escapePressed()) break;
		}
		lsr.dispose();
		new TextWindow("FlowSight Measurements",labels,sb.toString(),400,200);
	}

	public ImageProcessor processIP(FloatProcessor fp,float minval,float maxval,boolean raw){
		if(raw){
			fp.setMinAndMax(0.0,(double)maxval);
			return fp;
		} else {
			//need to convert to the "nice" lookup table
			byte[][] lut=lututils.nice_lut(false);
			float[] pix=(float[])fp.getPixels();
			int[] newpix=lututils.applyLUT(pix,minval,maxval,lut);
			return new ColorProcessor(fp.getWidth(),fp.getHeight(),newpix);
		}
	}

	public void saveImp(ImagePlus imp,String outdir,boolean keep){
		FileSaver fs=new FileSaver(imp);
		fs.saveAsTiff(outdir+imp.getTitle());
		if(!keep) imp.close();
	}

	public boolean checkImage(float[] params,float[][] limits){
		if(limits==null) return false;
		boolean selected=true;
		for(int i=0;i<params.length;i++){
			if(params[i]<limits[i][0] || params[i]>limits[i][1]){selected=false; break;}
		}
		return selected;
	}

	public boolean hasBud(float[] params,float[] aspectlims,float[] arealims){
		float area=params[0];
		float aspect=params[9];
		if(aspectlims==null) return false;
		if(arealims==null) return false;
		if(area>arealims[0] && area<arealims[1] && aspect>aspectlims[0] && aspect<aspectlims[1]) return true;
		return false;
	}

	public Object[] analyzeImage(ImagePlus imp1,ImagePlus imp2,int maskch,int transch,int numch,int dench,float threshmult,float smoothstdev){
		//always use the last channel for morphometrics
		Object[] stack=jutils.stack2array(imp1.getStack());
		int width=imp1.getWidth(); int height=imp1.getHeight();
		Object[] masks=jutils.stack2array(imp2.getStack());
		findblobs3 fb=new findblobs3(width,height);
		float[] object=fb.dofindblobs((byte[])masks[maskch]);
		float[] back=new float[object.length];
		for(int i=0;i<back.length;i++) if(object[i]==0.0f) back[i]=1.0f;
		Polygon outline=fb.get_object_outline(object,1);
		float[] grmsd=gradRMSD(stack[transch],outline,width,height); //5 params with area first
		float[] ellipse=ellipseParams(outline); //6 parameters with circ last
		float[] params=(float[])algutils.combine_arrays(grmsd,ellipse);
		float[] backs=new float[stack.length];
		for(int i=0;i<stack.length;i++) backs[i]=fb.get_object_stats(back,1,stack[i],"Avg");
		float maskarea=fb.get_object_stats(object,1,stack[0],"Count");
		int offset=params.length;
		float[] temp=(float[])algutils.expand_array(params,offset+stack.length+1);
		for(int i=0;i<stack.length;i++){
			float sum=fb.get_object_stats(object,1,stack[i],"Avg");
			temp[i+offset]=sum-backs[i];
		}
		temp[offset+stack.length]=temp[offset+numch]/temp[offset+dench];
		float denthresh=threshmult*temp[offset+dench];
		float[] heatmap=getHeatMap(stack[numch],stack[dench],object,width,height,backs[numch],backs[dench],denthresh,smoothstdev);
		return new Object[]{temp,new FloatProcessor(width,height,heatmap,null)};
	}

	public float[] getHeatMap(Object numimg1,Object denimg1,float[] object,int width,int height,float numback,float denback,float denthresh,float smoothstdev){
		float[] numimg=algutils.convert_arr_float(numimg1);
		float[] denimg=algutils.convert_arr_float(denimg1);
		if(denthresh<=0.0f) return new float[width*height];
		//first subtrack the backgrounds
		for(int i=0;i<numimg.length;i++) numimg[i]-=numback;
		for(int j=0;j<denimg.length;j++) denimg[j]-=denback;
		//now smooth
		jsmooth.blur2D(numimg,smoothstdev,width,height);
		jsmooth.blur2D(denimg,smoothstdev,width,height);
		float[] heatmap=new float[width*height];
		for(int i=0;i<numimg.length;i++){
			if(object[i]>0.0f && denimg[i]>=denthresh) heatmap[i]=numimg[i]/denimg[i];
		}
		return heatmap;
	}

	public float[] ellipseParams(Polygon outline){
		//find the centroid
		//float[] centroid=measure_object.centroid(outline);
		if(outline==null) return new float[6];
		float[] ellipseparams=measure_object.get_ellipse_parameters(outline);
		//these are x,y,angle,major,minor
		//change minor to aspect ratio
		ellipseparams[4]/=ellipseparams[3];
		float circularity=measure_object.circularity(outline);
		float[] params=(float[])algutils.expand_array(ellipseparams,ellipseparams.length+1);
		params[params.length-1]=circularity;
		return params;
	}

	public float[] gradRMSD(Object pix,Polygon outline1,int width,int height){
		float[] pix2=algutils.convert_arr_float2(pix);
		float[] sobel=(new jsobel(width,height)).do_sobel(pix2)[0];
		Polygon outline=new Polygon(new int[]{0,width-1,width-1,0},new int[]{0,0,height-1,height-1},4);
		if(outline1!=null) outline=outline1;
		float cnt=jstatistics.getstatistic("Count",pix,width,height,outline,null);
		float avg=jstatistics.getstatistic("Avg",pix,width,height,outline,null);
		float sobelcv=jstatistics.getstatistic("stdev",sobel,width,height,outline,null);
		sobelcv/=avg;
		float sobelmax=jstatistics.getstatistic("Max",pix,width,height,outline,null);
		float sobelavg=jstatistics.getstatistic("Avg",pix,width,height,outline,null);
		//note pixel size is around 330 nm per pixel (area is 0.11 um^2 per pixel)
		return new float[]{cnt*0.11f,avg,sobelcv,sobelavg,sobelmax};
	}

}
