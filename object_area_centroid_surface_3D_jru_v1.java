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
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;
import ij.text.*;

public class object_area_centroid_surface_3D_jru_v1 implements PlugIn, gui_interface {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Output_Indexed",false);
		gd.addNumericField("Number of surface edge voxels",3,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean outindexed=gd.getNextBoolean();
		int surfacethresh=(int)gd.getNextNumber();
		ImageStack objstack=imp.getStack();
		Object[] objpix=jutils.stack2array(objstack);
		boolean indexed=(objpix[0] instanceof float[]);
		findblobs3D fb=new findblobs3D(imp.getWidth(),imp.getHeight(),objstack.getSize(),this);
		float[][] objects=null;
		if(!indexed){
			IJ.showStatus("indexing objects");
			objects=fb.dofindblobs(objpix);
			if(outindexed){
				new ImagePlus("Indexed Objects 3D",jutils.array2stack(objects,imp.getWidth(),imp.getHeight())).show();
			}
		} else {
			objects=new float[objpix.length][];
			for(int i=0;i<objpix.length;i++) objects[i]=(float[])objpix[i];
			fb.set_objects(objects);
		}
		IJ.showStatus("measuring centroid, area, surface");
		float[][] centareas=fb.getCentroidsAreasSurface(objects,surfacethresh);
		//IJ.showStatus("measuring longest dimension");
		//float[] longestdims=fb.getLongestDimensions(objects);
		IJ.showStatus("outputting measurements");
		String headings="object\tx\ty\tz\tvolume\tsurface";
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<fb.nobjects;i++){
			sb.append(""+i);
			sb.append("\t"+centareas[i][0]+"\t"+centareas[i][1]+"\t"+centareas[i][2]+"\t"+centareas[i][3]+"\t"+centareas[i][4]+"\n");
		}
		new TextWindow("Object Measurements",headings,sb.toString(),400,200);
		IJ.showStatus("finished");
	}

	public void showMessage(String message){}

	public void showProgress(int currpos,int finalpos){IJ.showProgress(currpos,finalpos);}

}
