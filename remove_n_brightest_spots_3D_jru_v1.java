/*******************************************************************************
 * Copyright (c) 2021 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jseg.*;
import jguis.*;

public class remove_n_brightest_spots_3D_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin delete's the N brightest spots in a 3D image
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Number of spots",10,0);
		gd.addNumericField("Threshold (fraction of max)",0.05,5,15,null);
		gd.addNumericField("Spheroid_xy_diameter (pixels)",15,0);
		gd.addNumericField("Spheroid_z_diameter (slices)",5,0);
		gd.addNumericField("Replacement value",0.0,5,15,null);
		gd.addCheckbox("Output spots",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int nspots=(int)gd.getNextNumber();
		float threshfrac=(float)gd.getNextNumber();
		int minsep=(int)gd.getNextNumber();
		int zsep=(int)gd.getNextNumber();
		float fillval=(float)gd.getNextNumber();
		boolean output=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int currchan=imp.getC()-1;
		int currframe=imp.getT()-1;
		ImageStack stack=imp.getStack();
		Object[] stack2=jutils.get3DZSeries(stack,currchan,currframe,frames,slices,channels);
		//criteria are 0minarea, 1maxarea, 2searchd(nu), 3maxblobs, 4thresh, 5minsep, 6edgebuf
		float[] criteria={0.0f,10000.0f,0.0f,(float)nspots,0.05f,(float)minsep,0.0f};
		//FindBlobsFast fb=new FindBlobsFast(width,height,256,criteria);
		findblobs fb=new findblobs(width,height,criteria);
		fb.usemaxpt=true;
		float[][] pixels=algutils.convert_arr_float2(stack2);
		String fracstat="max";
		float[] spectrum=jstatistics.getspectrum(fracstat,pixels,null);
		//assume that the statistic on the spectrum is the same as on the whole stack
		float max=jstatistics.getstatistic(fracstat,spectrum,null);
		fb.thresh=threshfrac*max;
		float zedgebuf=0.0f;
		float searchrz=0.5f*(float)zsep;
		float[][] outobjects=new float[slices][width*height];
		float[][] blobstats=fb.dofindblobs3D(pixels,outobjects,searchrz,zedgebuf);
		//now mask out the objects based on the blob stats
		int dtype=algutils.get_array_type(stack2[0]);
		for(int i=0;i<outobjects.length;i++){
			for(int j=0;j<outobjects[i].length;j++){
				if(outobjects[i][j]>0.0f){
					algutils.setArrVal(stack2[i],fillval,j,dtype);
				}
			}
		}
		imp.updateAndDraw();
		if(output){
			new ImagePlus("Objects",jutils.array2stack(outobjects,width,height)).show();
		}
	}

}
