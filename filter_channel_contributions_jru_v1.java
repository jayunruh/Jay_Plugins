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
import java.awt.*;
import ij.plugin.*;
import jalgs.*;
import jguis.*;

public class filter_channel_contributions_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Channel to Filter",1,0);
		gd.addNumericField("Fraction Contribution Threshold",0.05,5,15,null);
		gd.addCheckbox("Output Fraction Contribution",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int filtchan=(int)gd.getNextNumber()-1;
		float fracthresh=(float)gd.getNextNumber();
		boolean outcontr=gd.getNextBoolean();
		int dtype=algutils.get_array_type(stack.getPixels(1));
		float[][] contr=new float[nslices*nframes][];
		for(int i=0;i<nslices*nframes;i++){
			Object[] chanstack=jutils.get3DCSeries(stack,i,0,1,nslices*nframes,nchans);
			if(outcontr) contr[i]=new float[width*height];
			for(int j=0;j<width*height;j++){
				Object spectrum=algutils.get_stack_col(chanstack,width,height,j,nchans);
				float[] spec2=algutils.convert_arr_float(spectrum);
				float sum=jstatistics.getstatistic("Sum",spec2,null);
				if(sum>0.0f){
					float frac=spec2[filtchan]/sum;
					if(outcontr) contr[i][j]=frac;
					if(frac<fracthresh) algutils.setArrVal(chanstack[filtchan],0.0f,j,dtype);
				}
			}
			IJ.showProgress(i,nslices*nframes);
		}
		imp.updateAndDraw();
		if(outcontr){
			new ImagePlus("Fractional Contribution",jutils.array2stack(contr,width,height)).show();
		}
	}

}
