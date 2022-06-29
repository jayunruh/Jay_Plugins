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

public class flowsight_stack_importer_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open CIF Image...", arg);
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
		float[][] limits=null;
		StringBuffer sb=new StringBuffer();

		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Mask_channel (consecutive)",3,0);
		gd.addNumericField("Trans_channel (consecutive)",3,0);
		gd.addStringField("Channels_to_import","1,2");
		gd.addCheckbox("Subtract_Background (outside mask)",true);
		gd.addNumericField("Pad_Image_Values",0,5,15,null);
		int maxdim=imp1.getHeight();
		if(imp1.getWidth()>maxdim) maxdim=imp1.getWidth();
		gd.addNumericField("Max_Image_Dimension",maxdim,0);
		gd.addCheckbox("Dont_Import_Images",false);
		gd.addNumericField("Max Images (0 for all)",0,0);
		gd.addCheckbox("Show_Mask (last channel)",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int maskch=(int)gd.getNextNumber()-1;
		int transch=(int)gd.getNextNumber()-1;
		String channels=gd.getNextString();
		String[] channels2=table_tools.split(channels,",");
		int[] selchan=new int[channels2.length];
		for(int i=0;i<selchan.length;i++) selchan[i]=(int)Float.parseFloat(channels2[i])-1;
		boolean backsub=gd.getNextBoolean();
		float padval=(float)gd.getNextNumber();
		maxdim=(int)gd.getNextNumber();
		boolean noimages=gd.getNextBoolean();
		int maximgs=(int)gd.getNextNumber();
		if(maximgs==0) maximgs=lsr.nseries/2;
		boolean showmask=gd.getNextBoolean();
		int showchan=selchan.length;
		if(showmask) showchan+=1;

		ImageStack outstack=new ImageStack(maxdim,maxdim);

		Object[] temp=analyzeImage(imp1,imp2,maskch,transch,selchan,backsub,padval,maxdim,noimages);
		int counter=0;
		sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
		if(!noimages) for(int j=0;j<showchan;j++) outstack.addSlice("",((float[][])temp[1])[j]);
		for(int i=1;i<maximgs;i++){
			imp1=(ImagePlus)lsr.getNextFrame();
			imp2=(ImagePlus)lsr.getNextFrame();
			if(i%500==0){
				IJ.log(imp1.getTitle());
				IJ.log(imp2.getTitle());
			}
			temp=analyzeImage(imp1,imp2,maskch,transch,selchan,backsub,padval,maxdim,noimages);
			counter++;
			sb.append(""+(counter+1)+"\t"+table_tools.print_float_array((float[])temp[0])+"\n");
			if(!noimages) for(int j=0;j<showchan;j++) outstack.addSlice("",((float[][])temp[1])[j]);
			IJ.showProgress(i,lsr.nseries/2);
			if(IJ.escapePressed()) break;
		}
		lsr.dispose();
		new TextWindow("FlowSight Measurements",labels,sb.toString(),400,200);
		if(!noimages){
			ImagePlus outimp=new ImagePlus(fname,outstack);
			outimp.setOpenAsHyperStack(true);
			outimp.setDimensions(showchan,1,counter+1);
			new CompositeImage(outimp,CompositeImage.COLOR).show();
		}
	}

	public Object[] analyzeImage(ImagePlus imp1,ImagePlus imp2,int maskch,int transch,int[] outchans,boolean backsub,float fillval,int maxdim,boolean noimages){
		//produces a background subtracted (optional) multi-channel image of maxdim x maxdim dimension with the last channel being the mask
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
		float[] temp=(float[])algutils.expand_array(params,offset+stack.length);
		for(int i=0;i<stack.length;i++){
			float avg=fb.get_object_stats(object,1,stack[i],"Avg");
			temp[i+offset]=avg-backs[i];
		}
		if(!noimages){
			float[][] images=new float[outchans.length+1][maxdim*maxdim];
			float xshift=0.5f*(float)(maxdim-width);
			float yshift=0.5f*(float)(maxdim-height);
			for(int i=0;i<outchans.length;i++){
				for(int j=0;j<maxdim*maxdim;j++) images[i][j]=fillval;
				int selch=outchans[i];
				float[] tempch=algutils.convert_arr_float(stack[selch]);
				if(backsub) for(int j=0;j<tempch.length;j++) tempch[j]-=backs[selch];
				interpolation.shift_copy_image(tempch,width,height,images[i],maxdim,maxdim,xshift,yshift);
			}
			for(int j=0;j<maxdim*maxdim;j++) images[outchans.length][j]=fillval;
			interpolation.shift_copy_image(object,width,height,images[outchans.length],maxdim,maxdim,xshift,yshift);
			return new Object[]{temp,images};
		} else {
			return new Object[]{temp};
		}
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
