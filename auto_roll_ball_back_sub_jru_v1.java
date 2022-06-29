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

public class auto_roll_ball_back_sub_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Roll_Ball_Rad",10.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float rollballrad=(float)gd.getNextNumber();
		float[] rollballmin=new float[slices];
		float[][] subtracted=new float[slices][];
		for(int i=0;i<slices;i++){
			subtracted[i]=algutils.convert_arr_float(stack.getPixels(i+1));
			float[] backsub=subtracted[i].clone();
			backsub=jutils.sub_roll_ball_back(backsub,rollballrad,width,height);
			for(int j=0;j<width*height;j++) backsub[j]=subtracted[i][j]-backsub[j];
			rollballmin[i]=jstatistics.getstatistic("Min",backsub,null);
			for(int j=0;j<backsub.length;j++) subtracted[i][j]-=rollballmin[i];
			IJ.showProgress(i,slices);
		}
		new PlotWindow4("Roll Ball Min Profile","frame","intensity",rollballmin).draw();
		jutils.create_hyperstack("Subtracted",jutils.array2stack(subtracted,width,height),imp).show();
	}

}
