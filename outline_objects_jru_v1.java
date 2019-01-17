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
import java.awt.event.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import jalgs.*;
import jalgs.jseg.*;
import jguis.*;
import java.util.*;

public class outline_objects_jru_v1 implements PlugIn, DialogListener{
	public ImagePlus display;
	public float[] pixels;
	public findblobs3 fb;
	public int width,height,minarea,maxarea;
	public float dispmin,dispmax,thresh,minint,maxint;
	public boolean showmask,editobj,segstack,measure_immediate;
	public float[][] offmult;
	public float[][] chstats;
	public boolean combine;
	//private TextField nobjects;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		pixels=(float[])imp.getProcessor().convertToFloat().getPixels();
		dispmin=(float)imp.getProcessor().getMin();
		dispmax=(float)imp.getProcessor().getMax();
		width=imp.getWidth();
		height=imp.getHeight();
		int nchannels=imp.getNChannels();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int currchan=imp.getChannel()-1;
		if(nchannels>1){
			int slice=imp.getSlice()-1;
			int frame=imp.getFrame()-1;
			Object[] cstack=jutils.get3DCSeries(imp.getStack(),slice,frame,frames,slices,nchannels);
			pixels=collapse_multichannel(cstack,currchan,width*height);
			if(pixels==null) return;
		}
		GenericDialog gd1=new GenericDialog("Options");
		gd1.addCheckbox("Smooth Image",false);
		gd1.addNumericField("StDev",1.0,5,15,null);
		gd1.showDialog(); if(gd1.wasCanceled()){return;}
		boolean smooth=gd1.getNextBoolean();
		float stdev=(float)gd1.getNextNumber();
		if(smooth) jsmooth.blur2D(pixels,stdev,width,height);
		thresh=0.5f;
		minint=jstatistics.getstatistic("Min",pixels,null);
		maxint=jstatistics.getstatistic("Max",pixels,null);
		fb=new findblobs3(width,height);
		ImageStack dispstack=new ImageStack(width,height);
		dispstack.addSlice("",pixels);
		dispstack.addSlice("",new float[width*height]);
		display=new ImagePlus("Outlined Objects",dispstack);
		display.setOpenAsHyperStack(true);
		display.setDimensions(2,1,1);
		display.copyScale(imp);
		display=new CompositeImage(display,CompositeImage.COMPOSITE);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		graylut.min=dispmin;
		graylut.max=dispmax;
		((CompositeImage)display).setChannelLut(graylut,1);
		((CompositeImage)display).setDisplayRange(dispmin,dispmax);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		redlut.min=0.0;
		redlut.max=255.0;
		((CompositeImage)display).setChannelLut(redlut,2);
		display.show();
		showmask=false;
		update_image(thresh,10,10000);
		//GenericDialog gd=new NonBlockingGenericDialog("Threshholding Options");
		GenericDialog gd=new GenericDialog("Threshholding Options");
		gd.addNumericField("Threshold",thresh,5,15,null);
		gd.addNumericField("Min_Area",10,0);
		gd.addNumericField("Max_Area",10000,0);
		gd.addCheckbox("Show_Mask",showmask);
		gd.addCheckbox("Edit_Objects",true);
		gd.addCheckbox("Segment_Stack",false);
		gd.addCheckbox("Measure_Objects_Immediately",false);
		//gd.addNumericField("#_of_objects",fb.nobjects,0);
		//Vector numberfields=gd.getNumericFields();
		//nobjects=(TextField)numberfields.elementAt(3);
		gd.addDialogListener(this);
		gd.showDialog();
		if(gd.wasCanceled()){
			display.close();
			return;
		}
		thresh=(float)gd.getNextNumber();
		minarea=(int)gd.getNextNumber();
		maxarea=(int)gd.getNextNumber();
		showmask=gd.getNextBoolean();
		editobj=gd.getNextBoolean();
		segstack=gd.getNextBoolean();
		measure_immediate=gd.getNextBoolean();
		update_image(thresh,minarea,maxarea);
		if(segstack){
			ImageStack threshstack=new ImageStack(width,height);
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					Object[] cstack=jutils.get3DCSeries(imp.getStack(),j,i,frames,slices,nchannels);
					float[] pix=combine_channels(cstack,currchan,width*height);
					if(smooth) jsmooth.blur2D(pix,stdev,width,height);
					minint=jstatistics.getstatistic("Min",pix,null);
					maxint=jstatistics.getstatistic("Max",pix,null);
					float[] temp=get_mask(pix,thresh,minarea,maxarea);
					threshstack.addSlice("",temp);
				}
			}
			jutils.create_hyperstack("Stack Threshold",threshstack,frames,slices,1,false,null).show();
		}
		if(measure_immediate){
			float thresh2=thresh*(maxint-minint)+minint;
			float[] objects=fb.dofindblobs(pixels,thresh2);
			int[] arealims={minarea,maxarea};
			fb.clear_edges(objects);
			fb.filter_area(objects,arealims,true);
			ImagePlus[] measimp=jutils.selectImages(false,1,new String[]{"Measurement_Image"});
			if(measimp!=null){
				ImageStack stack=measimp[0].getStack();
				Object[] stack2=jutils.stack2array(stack);
				float[][] stats=fb.get_all_object_stats(objects,stack2,"Avg");
				String newlabels=table_tools.createcollabels(stack2.length,"Avg");
				table_tools.create_table("Object_Measurements",stats,table_tools.split_string_tab(newlabels));
			}
		}
		if(editobj){
			threshold_panel tp=new threshold_panel();
			float thresh2=thresh*(maxint-minint)+minint;
			float[] objects=fb.dofindblobs(pixels,thresh2);
			int[] arealims={minarea,maxarea};
			fb.clear_edges(objects);
			fb.filter_area(objects,arealims);
			tp.init(display,objects,showmask);
			threshold_panel.launch_frame(tp);
		}
	}

	public float[] collapse_multichannel(Object[] stack,int currchan,int length){
		//here we combine channel images by matching 25th and 75th percentiles
		get_multichan_offsets(stack,currchan);
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Combine Channels?",true);
		for(int i=0;i<stack.length;i++){
			gd.addNumericField("Ch"+(i+1)+"_offset",offmult[i][0],5,15,null);
			gd.addNumericField("Ch"+(i+1)+"_multiplier",offmult[i][1],5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		combine=gd.getNextBoolean();
		for(int i=0;i<stack.length;i++){
			offmult[i][0]=(float)gd.getNextNumber();
			offmult[i][1]=(float)gd.getNextNumber();
		}
		return combine_channels(stack,currchan,length);
	}

	public float[] combine_channels(Object[] stack,int currchan,int length){
		float[] original=algutils.convert_arr_float2(stack[currchan]);
		if(!combine){return original;}
		float[] temp=new float[length];
		for(int i=0;i<stack.length;i++){
			float[] temp2=null;
			if(i==currchan) temp2=original;
			temp2=algutils.convert_arr_float2(stack[i]);
			for(int j=0;j<length;j++){
				float sub=temp2[j]-chstats[i][0];
				temp[j]+=(offmult[i][1]*sub+offmult[i][0]+chstats[currchan][0])/(float)stack.length;
			}
		}
		return temp;
	}

	public void get_multichan_offsets(Object[] stack,int currchan){
		//here we combine channel images by matching 25th and 75th percentiles
		/*chstats=new float[stack.length][];
		for(int i=0;i<stack.length;i++){
			chstats[i]=new float[]{25.0f,75.0f};
			jstatistics.getstatistic("Percentile",stack[i],chstats[i]);
		}
		offmult=new float[stack.length][2];
		for(int i=0;i<stack.length;i++){
			if(i==currchan){
				offmult[i][1]=1.0f;
			} else {
				offmult[i][0]=chstats[currchan][0]-chstats[i][0];
				offmult[i][1]=(chstats[i][1]-chstats[i][0])/(chstats[currchan][1]-chstats[currchan][0]);
			}
		}*/
		//lets try matching the minimum value and the mean intensity
		chstats=new float[stack.length][2];
		for(int i=0;i<stack.length;i++){
			chstats[i][0]=jstatistics.getstatistic("Min",stack[i],null);
			chstats[i][1]=jstatistics.getstatistic("Max",stack[i],null);
		}
		offmult=new float[stack.length][2];
		for(int i=0;i<stack.length;i++){
			if(i==currchan){
				offmult[i][1]=1.0f;
			} else {
				offmult[i][0]=chstats[i][0];
				offmult[i][1]=(chstats[currchan][1]-chstats[currchan][0])/(chstats[i][1]-chstats[i][0]);
			}
		}
	}

	public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		thresh=(float)gd.getNextNumber();
		minarea=(int)gd.getNextNumber();
		maxarea=(int)gd.getNextNumber();
		showmask=gd.getNextBoolean();
		editobj=gd.getNextBoolean();
		segstack=gd.getNextBoolean();
		measure_immediate=gd.getNextBoolean();
		update_image(thresh,minarea,maxarea);
		//nobjects.setText(""+fb.nobjects);
		return true;
	}

	public void update_image(float thresh,int minarea,int maxarea){
		float thresh2=thresh*(maxint-minint)+minint;
		float[] objects=fb.dofindblobs(pixels,thresh2);
		int[] arealims={minarea,maxarea};
		fb.clear_edges(objects);
		fb.filter_area(objects,arealims);
		ImageProcessor ip=display.getStack().getProcessor(2);
		float[] disppix=(float[])ip.getPixels();
		if(showmask){
			float[] temp=(float[])jutils.convert_array(fb.tobinary(objects,true),2);
			System.arraycopy(temp,0,disppix,0,disppix.length);
		} else {
			for(int i=0;i<disppix.length;i++){
				disppix[i]=0.0f;
			}
			Polygon[] objects2=fb.get_object_outlines(objects);
			ip.setColor(Color.white);
			for(int i=0;i<objects2.length;i++){
				jutils.draw_polygon(ip,objects2[i],true);
			}
		}
		display.updateAndDraw();
	}

	public float[] get_mask(float[] pix,float thresh,int minarea,int maxarea){
		float thresh2=thresh*(maxint-minint)+minint;
		float[] objects=fb.dofindblobs(pix,thresh2);
		int[] arealims={minarea,maxarea};
		fb.clear_edges(objects);
		fb.filter_area(objects,arealims);
		return jutils.convert_arr_float(fb.tobinary(objects,true));
	}

}
