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

public class STICS_map_jru_v3 implements PlugIn {

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
		int rlength=10;
		gd.addNumericField("Running_avg_length",rlength,0);
		int inc=5;
		gd.addNumericField("Start_frame_increment",inc,0);
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
		rlength=(int)gd.getNextNumber();
		inc=(int)gd.getNextNumber();

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
		if(roi==null){
			roi=new Roi(0,0,width,height);
		}

		STICS_map map=new STICS_map(subsize,stepsize);
		Object[] tseries=jutils.get3DTSeries(stack,0,0,frames,slices,channels);
		map.update_STICS_map(tseries,width,height,0,rlength,roi.getPolygon(),shift);
		FloatProcessor fp=map.get_map(scaling,ftime,stepsize,centered,norm,multiplier,stepsize,magthresh);
		ImageStack vector_stack=new ImageStack(fp.getWidth(),fp.getHeight());
		vector_stack.addSlice("",fp);
		float[][] vel=map.get_scaled_velocities(scaling,ftime,stepsize);
		ImageStack velstack=new ImageStack(map.xregions,map.yregions);
		velstack.addSlice("",vel[0]);
		velstack.addSlice("",vel[1]);
		int velframes=2;
		IJ.showStatus("frame "+0+" calculated");
		for(int i=inc;i<(frames-rlength);i+=inc){
			map.update_STICS_map(tseries,width,height,i,rlength,roi.getPolygon(),shift);
			FloatProcessor fp2=map.get_map(scaling,ftime,stepsize,centered,norm,multiplier,stepsize,magthresh);
			vector_stack.addSlice("",fp2);
			vel=map.get_scaled_velocities(scaling,ftime,stepsize);
			velstack.addSlice("",vel[0]);
			velstack.addSlice("",vel[1]);
			velframes+=2;
			IJ.showStatus("frame "+i+" calculated");
		}
		(new ImagePlus("STICS Vectors",vector_stack)).show();
		ImagePlus imp3=new ImagePlus("Velocities",velstack);
		imp3.setOpenAsHyperStack(true);
		imp3.setDimensions(2,1,velframes/2);
		new CompositeImage(imp3,CompositeImage.COLOR).show();
	}

}
