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
import jalgs.*;
import jalgs.jfft.*;
import jalgs.jseg.*;
import jguis.*;
import ij.measure.*;
import ij.text.*;

public class STICS_map_jru_v2 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		Calibration cal=imp.getCalibration();
		GenericDialog gd=new GenericDialog("Options");
		int subsize=32;
		gd.addNumericField("Subregion Size (pixels)?",subsize,0);
		int stepsize=16;
		gd.addNumericField("Step Size?",stepsize,0);
		int shift=3;
		gd.addNumericField("STICS temporal Shift?",shift,0);
		float xoffset=0.0f;
		gd.addNumericField("X_Offset",xoffset,5,15,null);
		float yoffset=0.0f;
		gd.addNumericField("Y_Offset",yoffset,5,15,null);
		float multiplier=8.0f;
		gd.addNumericField("Velocity Multiplier",multiplier,5,15,null);
		float ftime=1.0f;
		gd.addNumericField("Frame_Time(min)",ftime,5,15,null);
		float scaling=(float)cal.pixelWidth;
		gd.addNumericField("Pixel_Size(um)",scaling,5,15,null);
		boolean norm=true;
		gd.addCheckbox("Normalize_Vector_lengths?",norm);
		boolean centered=true;
		gd.addCheckbox("Center_Vectors?",centered);
		float magthresh=0.0f;
		gd.addNumericField("Magnitude_Threshhold?",magthresh,5,15,null);
		boolean outvel=false;
		gd.addCheckbox("Output_Velocities?",outvel);
		boolean maskroi=false;
		gd.addCheckbox("Use_Movie_Mask",maskroi);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		subsize=(int)gd.getNextNumber();
		stepsize=(int)gd.getNextNumber();
		shift=(int)gd.getNextNumber();
		xoffset=(float)gd.getNextNumber();
		yoffset=(float)gd.getNextNumber();
		multiplier=(float)gd.getNextNumber();
		ftime=(float)gd.getNextNumber();
		scaling=(float)gd.getNextNumber();
		norm=gd.getNextBoolean();
		centered=gd.getNextBoolean();
		magthresh=(float)gd.getNextNumber();
		outvel=gd.getNextBoolean();
		maskroi=gd.getNextBoolean();

		int width=imp.getWidth();
		int xregions=1+(int)(((float)width-(float)subsize)/(float)stepsize);
		int newwidth=xregions*subsize;
		int height=imp.getHeight();
		int yregions=1+(int)(((float)height-(float)subsize)/(float)stepsize);
		int newheight=yregions*subsize;
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		int frames=imp.getNFrames();
		if(frames==1){frames=slices; slices=1;}

		Roi roi=imp.getRoi();
		if(roi==null) roi=new Roi(0,0,width,height);
		STICS_map map=new STICS_map(subsize,stepsize);
		map.linewidth=2;
		if(maskroi){
			ImagePlus maskimp=jutils.selectImages(false,1,new String[]{"Mask"})[0];
			ImageStack maskstack=maskimp.getStack();
			boolean[][] mask=new boolean[maskstack.getSize()][width*height];
			for(int i=0;i<mask.length;i++){
				float[] pix=algutils.convert_arr_float2(maskstack.getPixels(i+1));
				for(int j=0;j<pix.length;j++) mask[i][j]=(pix[j]>0.0f);
			}
			map.update_STICS_map(jutils.get3DTSeries(stack,0,0,frames,slices,channels),width,height,0,frames,mask,shift);
		} else {
			map.update_STICS_map(jutils.get3DTSeries(stack,0,0,frames,slices,channels),width,height,0,frames,roi.getPolygon(),shift);
		}
		FloatProcessor fp=map.get_map(scaling,ftime,stepsize,centered,norm,multiplier,stepsize,magthresh);
		ImageStack vector_stack=new ImageStack(fp.getWidth(),fp.getHeight());
		ImageStack velocity_stack=null;
		vector_stack.addSlice("",fp);
		float[][] scaled_velocities=null;
		TextWindow tw=null;
		if(outvel){
			velocity_stack=new ImageStack(map.xregions,map.yregions);
			scaled_velocities=map.get_scaled_velocities(scaling,ftime,stepsize);
			for(int i=0;i<2;i++) velocity_stack.addSlice("",scaled_velocities[i]);
			tw=new TextWindow("STICS Values","xpos\typos\txvel\tyvel","",400,200);
			for(int i=0;i<map.xregions*map.yregions;i++)
				tw.append(scaled_velocities[3][i]+"\t"+scaled_velocities[4][i]+"\t"+scaled_velocities[0][i]+"\t"+scaled_velocities[1][i]+"\n");
		}
		for(int i=1;i<channels*slices;i++){
			map.update_STICS_map(jutils.get3DTSeries(stack,0,i,frames,1,channels*slices),width,height,frames,roi.getPolygon(),shift);
			FloatProcessor fp2=map.get_map(scaling,ftime,stepsize,centered,norm,multiplier,stepsize,magthresh);
			vector_stack.addSlice("",fp2);
			if(outvel){
				scaled_velocities=map.get_scaled_velocities(scaling,ftime,stepsize);
				for(int j=0;j<2;j++) velocity_stack.addSlice("",map.velocities[j]);
				for(int j=0;j<map.xregions*map.yregions;j++)
					tw.append(scaled_velocities[3][j]+"\t"+scaled_velocities[4][j]+"\t"+scaled_velocities[0][j]+"\t"+scaled_velocities[1][j]+"\n");
			}
		}
		(new ImagePlus("STICS Vectors",vector_stack)).show();
		if(outvel){
			ImagePlus imp3=new ImagePlus("Velocities",velocity_stack);
			imp3.setOpenAsHyperStack(true);
			imp3.setDimensions(2,slices*channels,1);
			new CompositeImage(imp3,CompositeImage.COLOR).show();
		}
	}

}
