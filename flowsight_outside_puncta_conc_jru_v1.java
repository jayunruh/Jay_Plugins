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

public class flowsight_outside_puncta_conc_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open CIF Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		//DirectoryChooser dc=new DirectoryChooser("Save Directory");
		//String outdir=dc.getDirectory();
		//if(outdir==null) return;
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
		labels+="\tpunctaarea";
		labels+="\tpunctaavg";
		labels+="\tcytoavg";	
		labels+="\tfretratio";
		labels+="\tvolume";
		labels+="\ttotconc";
		labels+="\tcytoconc";
		labels+="\tpunctafretratio";
		labels+="\tcytofretratio";
		float[][] limits=null;
		StringBuffer sb=new StringBuffer();

		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Trans_channel (1-6)",3,0);
		gd.addNumericField("Mask_channel (1-6)",3,0);
		gd.addNumericField("numerator_channel (1-6)",2,0);
		gd.addNumericField("denominator_channel (1-6)",6,0);
		gd.addNumericField("Div_by_stdev",10.0f,5,15,null);
		gd.addNumericField("Roll_ball_rad",3.0f,5,15,null);
		gd.addNumericField("Processed_threshold",3.0f,5,15,null);
		gd.addNumericField("Cyto_circ_radius (away from puncta)",2,0);
		gd.addNumericField("Cyto_circ_gap (away from puncta)",1,0);
		gd.addNumericField("number of images to show",10,0);
		gd.addCheckbox("Show log",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int transch=(int)gd.getNextNumber()-1;
		int maskch=(int)gd.getNextNumber()-1;
		int numch=(int)gd.getNextNumber()-1;
		int dench=(int)gd.getNextNumber()-1;
		float divbystdev=(float)gd.getNextNumber();
		float rollballrad=(float)gd.getNextNumber();
		float procthresh=(float)gd.getNextNumber();
		int circrad=(int)gd.getNextNumber();
		int circgap=(int)gd.getNextNumber();
		int numshow=(int)gd.getNextNumber();
		boolean showlog=gd.getNextBoolean();
		Object[] temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,divbystdev,rollballrad,procthresh,circrad,circgap);
		int counter=0;
		sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
		List<ImagePlus> selected=new ArrayList<ImagePlus>();
		//if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
		//ImagePlus tempimp=new ImagePlus("Mask1.tif",new FloatProcessor(imp1.getWidth(),imp1.getHeight(),(float[])temp[1],null));
		ImagePlus tempimp=makeImage("Mask1.tif",(float[])temp[1],(float[])temp[2],(float[])temp[3],imp1.getWidth(),imp1.getHeight());
		if(numshow>0) tempimp.show();
		//saveImp(tempimp,outdir,numshow>0);
		for(int i=1;i<lsr.nseries/2;i++){
			imp1=(ImagePlus)lsr.getNextFrame();
			imp2=(ImagePlus)lsr.getNextFrame();
			if(showlog){	
				IJ.log(imp1.getTitle());
				IJ.log(imp2.getTitle());
			} else {
				System.out.println(imp1.getTitle());
				System.out.println(imp2.getTitle());
			}
			temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,divbystdev,rollballrad,procthresh,circrad,circgap);
			counter++;
			sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
			//if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
			//ImagePlus tempimp2=new ImagePlus("Mask"+(i+1)+".tif",new FloatProcessor(imp1.getWidth(),imp1.getHeight(),(float[])temp[1],null));
			ImagePlus tempimp2=makeImage("Mask"+(i+1)+".tif",(float[])temp[1],(float[])temp[2],(float[])temp[3],imp1.getWidth(),imp1.getHeight());
			if(i<numshow) tempimp2.show();
			//saveImp(tempimp2,outdir,i<=numshow);
			IJ.showProgress(i,lsr.nseries/2);
			if(IJ.escapePressed()) break;
		}
		lsr.dispose();
		new TextWindow("FlowSight Measurements",labels,sb.toString(),400,200);
	}

	public ImagePlus makeImage(String name,float[] objects,float[] image,float[] objects2,int width,int height){
		ImageStack stack=new ImageStack(width,height);
		stack.addSlice("",image);
		stack.addSlice("",objects);
		stack.addSlice("",objects2);
		ImagePlus imp=new ImagePlus(name,stack);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(3,1,1);
		return new CompositeImage(imp,CompositeImage.COLOR);
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

	public Object[] analyzeImage(ImagePlus imp1,ImagePlus imp2,int maskch,int transch,int numch,int dench,float divbystdev,float rollballrad,float procthresh,int circrad,int circgap){
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
		float[] temp=(float[])algutils.expand_array(params,offset+stack.length+9);
		for(int i=0;i<stack.length;i++){
			float sum=fb.get_object_stats(object,1,stack[i],"Avg");
			temp[i+offset]=sum-backs[i];
		}
		Object[] punctastats=getPunctaStats(stack[dench],stack[numch],object,width,height,backs[dench],backs[numch],divbystdev,rollballrad,procthresh,circrad,circgap);
		float[] pstats=(float[])punctastats[0];
		temp[offset+stack.length]=pstats[0]; //puncta area
		temp[offset+stack.length+1]=pstats[1]; //puncta avg
		temp[offset+stack.length+2]=pstats[2]; //cyto avg
		temp[offset+stack.length+3]=temp[offset+numch]/temp[offset+dench]; //fretratio
		temp[offset+stack.length+4]=(float)((4.0/3.0)*Math.PI*Math.pow(temp[0]/Math.PI,1.5)); //volume in microns cubed
		temp[offset+stack.length+5]=temp[offset+dench]*temp[0]/temp[offset+stack.length+4]; //totconc
		temp[offset+stack.length+6]=temp[offset+stack.length+2]*temp[0]/temp[offset+stack.length+4]; //cytoconc
		temp[offset+stack.length+7]=pstats[3]/pstats[1]; //punctafretratio
		temp[offset+stack.length+8]=pstats[4]/pstats[2]; //punctafretratio
		return new Object[]{temp,punctastats[1],punctastats[2],object};
	}

	public Object[] getPunctaStats(Object pix,Object numpix,float[] object,int width,int height,float back,float numback,float divbystdev,float rollballrad,float procthresh,int circrad,int circgap){
		float[] backsub=algutils.convert_arr_float(pix);
		for(int i=0;i<backsub.length;i++) backsub[i]-=back;
		float[] numbacksub=algutils.convert_arr_float(numpix);
		for(int i=0;i<numbacksub.length;i++) numbacksub[i]-=numback;
		float[] blurred=backsub.clone();
		jsmooth.blur2D(blurred,divbystdev,width,height);
		for(int i=0;i<width*height;i++) blurred[i]=backsub[i]/blurred[i];
		blurred=jutils.sub_roll_ball_back(blurred,rollballrad,width,height);
		for(int i=0;i<width*height;i++){
			if(object[i]<=0.0f) blurred[i]=0.0f;
		}
		findblobs3 fb=new findblobs3(width,height);
		float[] objects1=fb.dofindblobs(blurred,procthresh);
		//dilate by the gap radius
		float[] tempobj=objects1.clone();
		for(int i=0;i<circgap;i++) fb.dilateobjects(tempobj,false);
		float[] objects2=fb.get_circ(tempobj,circrad);
		float punctarea=0.0f;
		float circarea=0.0f;
		float punctavg=0.0f;
		float circavg=0.0f;
		float numpunctavg=0.0f;
		float numcircavg=0.0f;
		for(int i=0;i<backsub.length;i++){
			if(object[i]>0.0f){
				if(objects1[i]>0.0f){punctavg+=backsub[i]; numpunctavg+=numbacksub[i]; punctarea+=1.0f;}
				if(objects2[i]>0.0f){circavg+=backsub[i]; numcircavg+=numbacksub[i]; circarea+=1.0f;}
			} else {
				objects1[i]=0.0f;
				objects2[i]=0.0f;
			}
		}
		punctavg/=punctarea;
		circavg/=circarea;
		numpunctavg/=punctarea;
		numcircavg/=punctarea;
		float[] temp={punctarea,punctavg,circavg,numpunctavg,numcircavg};
		return new Object[]{temp,objects2,backsub,numbacksub};
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
		//return new float[]{cnt*0.11f,avg,sobelcv,sobelavg,sobelmax};
		return new float[]{cnt,avg,sobelcv,sobelavg,sobelmax};
	}

}
