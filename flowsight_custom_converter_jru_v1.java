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

public class flowsight_custom_converter_jru_v1 implements PlugIn {
	public boolean showbud;

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
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
		labels+="\tbudarea";
		for(int i=0;i<nch;i++){
			labels+="\tbudch"+(i+1);
		}
		labels+="\tbudfretratio";
		labels+="\tmotherarea";
		for(int i=0;i<nch;i++){
			labels+="\tmotherch"+(i+1);
		}
		labels+="\tmotherfretratio";
		float[][] limits=null;
		StringBuffer sb=new StringBuffer();
		float[] budaspectlims={0.58f,0.8f};
		float[] budarealims={60.0f,105.0f};
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_Diploid_Area",budarealims[0],5,15,null);
		gd.addNumericField("Max_Diploid_Area",budarealims[1],5,15,null);
		gd.addNumericField("Min_Diploid_Aspect",budaspectlims[0],5,15,null);
		gd.addNumericField("Max_Diploid_Aspect",budaspectlims[1],5,15,null);
		gd.addNumericField("Trans_channel (1-6)",6,0);
		gd.addNumericField("Mask_channel (1-6)",6,0);
		gd.addNumericField("numerator_channel (1-6)",2,0);
		gd.addNumericField("denominator_channel (1-6)",5,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		budarealims[0]=(float)gd.getNextNumber();
		budarealims[1]=(float)gd.getNextNumber();
		budaspectlims[0]=(float)gd.getNextNumber();
		budaspectlims[1]=(float)gd.getNextNumber();
		int transch=(int)gd.getNextNumber()-1;
		int maskch=(int)gd.getNextNumber()-1;
		int numch=(int)gd.getNextNumber()-1;
		int dench=(int)gd.getNextNumber()-1;
		showbud=true;
		float[] temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,budaspectlims,budarealims);
		int counter=0;
		sb.append(""+(counter+1)+"\t"+table_tools.print_float_array(temp)+"\n");
		List<ImagePlus> selected=new ArrayList<ImagePlus>();
		if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
		for(int i=1;i<lsr.nseries/2;i++){
			imp1=(ImagePlus)lsr.getNextFrame();
			imp2=(ImagePlus)lsr.getNextFrame();
			IJ.log(imp1.getTitle());
			IJ.log(imp2.getTitle());
			temp=analyzeImage(imp1,imp2,maskch,transch,numch,dench,budaspectlims,budarealims);
			counter++;
			sb.append(""+(counter+1)+"\t"+table_tools.print_float_array(temp)+"\n");
			if(checkImage(temp,limits)){selected.add(imp1); selected.add(imp2);}
			IJ.showProgress(i,lsr.nseries/2);
			if(IJ.escapePressed()) break;
		}
		lsr.dispose();
		new TextWindow("FlowSight Measurements",labels,sb.toString(),400,200);
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

	public float[] analyzeImage(ImagePlus imp1,ImagePlus imp2,int maskch,int transch,int numch,int dench,float[] aspectlims,float[] arealims){
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
		float[] params=(float[])algutils.combine_arrays(grmsd,ellipse); //this now has 11 parameters
		float[] backs=new float[stack.length];
		for(int i=0;i<stack.length;i++) backs[i]=fb.get_object_stats(back,1,stack[i],"Avg");
		float maskarea=fb.get_object_stats(object,1,stack[0],"Count");
		int offset=params.length;
		float[] temp=(float[])algutils.expand_array(params,offset+stack.length+1); //now we add stack.length+1 (sums and ratio) parameters
		for(int i=0;i<stack.length;i++){
			float sum=fb.get_object_stats(object,1,stack[i],"Sum");
			temp[i+offset]=sum-backs[i]*maskarea;
		}
		temp[offset+stack.length]=temp[offset+numch]/temp[offset+dench];
		offset=temp.length;
		temp=(float[])algutils.expand_array(temp,offset+stack.length+2); //add another stack.length+2 (budarea+sums+ratio) parameters
		if(hasBud(temp,aspectlims,arealims)){
			//calculate the bud parameters
			float[] temp2=measureBud(stack,ellipse,object,algutils.convert_arr_float(stack[transch]),width,height,backs);
			System.arraycopy(temp2,0,temp,offset,temp2.length);
			temp[offset+1+stack.length]=temp2[1+numch]/temp2[1+dench];
			//calculate the mother parameters
			offset=temp.length;
			temp=(float[])algutils.expand_array(temp,offset+stack.length+2); //now add another stack.length+2 (motherarea+sums+ratio) parameters
			temp[offset]=temp[0]-temp[11+stack.length+1]; //area-budarea
			for(int i=0;i<stack.length;i++){
				temp[offset+1+i]=temp[11+i]-temp[11+stack.length+2+i]; //int-budint
			}
			temp[offset+1+stack.length]=temp[offset+1+numch]/temp[offset+1+dench];
		} else {
			offset=temp.length;
			temp=(float[])algutils.expand_array(temp,offset+stack.length+2); //now add another stack.length+2 (motherarea+sums+ratio) parameters (zeros in this case)
		}
		return temp;
	}

	public float[] measureBud(Object[] stack,float[] ellipseparams,float[] object,float[] transimg,int width,int height,float[] backavg){
		Object[] budmeas=measure_object.getBudObject(ellipseparams,object,transimg,width,height);
		if(budmeas==null) return new float[stack.length+1];
		float[] budobj=(float[])budmeas[0];
		float[] budmeas2=(float[])budmeas[1];
		IJ.log(table_tools.print_float_array(budmeas2,1));
		if(showbud){
			Object[] temp={object.clone(),transimg.clone(),budobj.clone()};
			new ImagePlus("bud object",jutils.array2stack(temp,width,height)).show();
			showbud=false;
		}
		float[] measurements=new float[stack.length+1];
		findblobs3 fb=new findblobs3(width,height);
		measurements[0]=fb.get_object_stats(budobj,1,stack[0],"count");
		float budarea=measurements[0];
		measurements[0]*=0.11f;
		for(int i=0;i<stack.length;i++) measurements[i+1]=fb.get_object_stats(budobj,1,stack[i],"Sum")-budarea*backavg[i];
		return measurements;
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
