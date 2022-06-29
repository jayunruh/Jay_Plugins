/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jseg.*;
import ij.plugin.frame.RoiManager;

public class contour_selection_jru_v1 implements PlugIn {

	//this plugin creates a contour for the current channel at the appropriate fractional level and adds it to the roi manager
	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		Roi roi=imp.getRoi();
		if(roi==null) roi=new Roi(0,0,width,height);
		float[] pix=(float[])imp.getProcessor().convertToFloat().getPixels();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Fractional Contour Level (0-1)",0.75f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float fcont=(float)gd.getNextNumber();
		//clear all of the regions outside the roi
		float[] pix2=pix.clone();
		float max=pix2[0];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(!roi.contains(j,i)) pix2[j+i*width]=0.0f;
				else{
					if(pix2[j+i*width]>max) max=pix2[j+i*width];
				}
			}
		}
		//now segment the selected region at the specified level relative to max
		byte[] seg=new byte[width*height];
		for(int i=0;i<width*height;i++) if(pix2[i]>=fcont*max) seg[i]=(byte)255;
		//now outline that segment and add to roi manager
		//if there are multiple objects, add them all to the manager
		findblobs3 fb=new findblobs3(width,height);
		float[] objects=fb.dofindblobs(seg);
		Polygon[] outlines=fb.get_object_outlines(objects);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		for(int i=0;i<outlines.length;i++) rman.addRoi(new PolygonRoi(outlines[i],Roi.FREEROI));
	}

}
