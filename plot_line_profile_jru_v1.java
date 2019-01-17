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
import ij.text.*;
import ij.util.*;
import jguis.*;
import jalgs.*;
import ij.measure.*;

public class plot_line_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		int width = imp.getWidth();
		int height=imp.getHeight();
		int channels=imp.getNChannels();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int currz=imp.getZ()-1;
		int currt=imp.getT()-1;
		Object[] cstack=jutils.get3DCSeries(imp.getStack(),currz,currt,frames,slices,channels);
		//FloatProcessor fp = (FloatProcessor)imp.getProcessor();
		//float[] pixels = (float[])fp.getPixels();
		GenericDialog gd = new GenericDialog("Options");
		boolean vertline=false;
		gd.addCheckbox("Vertical?",vertline);
		int linenum=(int)((float)height/2.0f);
		gd.addNumericField("Line number",linenum,0);
		gd.addNumericField("Thickness (odd)",1,0);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		vertline=gd.getNextBoolean();
		linenum=(int)gd.getNextNumber();
		int thickness=(int)gd.getNextNumber();
		Calibration cal=imp.getCalibration();
		float psize=(float)cal.pixelWidth;
		if(vertline){
			psize=(float)cal.pixelHeight;
		}
		if(!vertline){
			float[][] xvals=new float[channels][width];
			float[][] linevals=new float[channels][];
			for(int i=0;i<width;i++){
				xvals[0][i]=psize*(float)i;
			}
			linevals[0]=getHorLine(cstack[0],linenum,thickness,width,height);
			for(int i=1;i<channels;i++){
				xvals[i]=xvals[0];
				linevals[i]=getHorLine(cstack[i],linenum,thickness,width,height);
			}
			PlotWindow4 plot = new PlotWindow4("Profile "+linenum,"x","Intensity",xvals,linevals,null);
			plot.draw();
		}
		else{
			float[][] xvals=new float[channels][height];
			float[][] linevals=new float[channels][];
			for(int i=0;i<height;i++){
				xvals[0][i]=psize*(float)i;
			}
			linevals[0]=getVertLine(cstack[0],linenum,thickness,width,height);
			for(int i=1;i<channels;i++){
				xvals[i]=xvals[0];
				linevals[i]=getVertLine(cstack[i],linenum,thickness,width,height);
			}
			PlotWindow4 plot = new PlotWindow4("Profile "+linenum,"x","Intensity",xvals,linevals,null);
			plot.draw();
		}
	}

	public float[] getHorLine(Object pix,int linenum,int thickness,int width,int height){
		if(thickness<=1) return algutils.convert_arr_float2(algutils.get_image_row(pix,width,height,linenum));
		int start=linenum-thickness/2;
		int end=start+thickness-1;
		if(start<0) start=0;
		if(end>(height-1)) end=height-1;
		float newthickness=(float)(end-start+1);
		float[][] rows=new float[(int)newthickness][];
		for(int i=start;i<=end;i++){
			rows[i-start]=algutils.convert_arr_float2(algutils.get_image_row(pix,width,height,i));
		}
		float[] avg=new float[width];
		for(int i=0;i<width;i++){
			for(int j=0;j<(int)newthickness;j++){
				avg[i]+=rows[j][i];
			}
		}
		for(int i=0;i<width;i++) avg[i]/=newthickness;
		return avg;
	}

	public float[] getVertLine(Object pix,int linenum,int thickness,int width,int height){
		if(thickness<=1) return algutils.convert_arr_float2(algutils.get_image_col(pix,width,height,linenum));
		int start=linenum-thickness/2;
		int end=start+thickness-1;
		if(start<0) start=0;
		if(end>(width-1)) end=height-1;
		float newthickness=(float)(end-start+1);
		float[][] cols=new float[(int)newthickness][];
		for(int i=start;i<=end;i++){
			cols[i-start]=algutils.convert_arr_float2(algutils.get_image_col(pix,width,height,i));
		}
		float[] avg=new float[width];
		for(int i=0;i<height;i++){
			for(int j=0;j<(int)newthickness;j++){
				avg[i]+=cols[j][i];
			}
		}
		for(int i=0;i<height;i++) avg[i]/=newthickness;
		return avg;
	}

}
