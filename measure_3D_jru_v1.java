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
import ij.plugin.frame.*;
import ij.measure.*;
import jguis.*;

public class measure_3D_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		float rsize=(float)jutils.get_psize(imp);
		float zsize=(float)jutils.get_pdepth(imp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Slice_Spacing",zsize,5,15,null);
		gd.addNumericField("Pixel_Spacing",rsize,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		zsize=(float)gd.getNextNumber();
		rsize=(float)gd.getNextNumber();
		new WaitForUserDialog("Select First Point").show();
		Roi roi1=imp.getRoi();
		int[] pt1=null;
		if(roi1 instanceof PointRoi){
			Rectangle r=roi1.getBounds();
			int[] temp={r.x,r.y,imp.getZ()};
			pt1=temp;
		} else {
			IJ.showMessage("Roi must be a point");
			return;
		}
		new WaitForUserDialog("Select Second Point").show();
		Roi roi2=imp.getRoi();
		int[] pt2=null;
		if(roi2 instanceof PointRoi){
			Rectangle r=roi2.getBounds();
			int[] temp={r.x,r.y,imp.getZ()};
			pt2=temp;
		} else {
			IJ.showMessage("Roi must be a point");
			return;
		}
		//IJ.log("pt1 = \t"+table_tools.print_int_array(pt1));
		//IJ.log("pt2 = \t"+table_tools.print_int_array(pt2));
		float dist=(float)Math.sqrt(((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+(pt2[1]-pt1[1])*(pt2[1]-pt1[1]))*rsize*rsize+(pt2[2]-pt1[2])*(pt2[2]-pt1[2])*zsize*zsize);
		IJ.log("Distance = \t"+dist);
	}

}
