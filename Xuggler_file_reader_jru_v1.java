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
import ij.io.*;
import jguis.*;

public class Xuggler_file_reader_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Movie...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		Xuggler_file_reader reader=new Xuggler_file_reader(directory,fname);
		double size=((double)reader.width*(double)reader.height*4.0*(double)reader.nframes);
		IJ.log("width = "+reader.width);
		IJ.log("height = "+reader.height);
		IJ.log("packets = "+reader.npackets);
		IJ.log("duration = "+reader.duration);
		IJ.log("est frame interval = "+reader.frameInterval);
		IJ.log("est # frames = "+reader.nframes);
		IJ.log("Est Size (MB) = "+(float)(size/1000000.0));
		GenericDialog gd=new GenericDialog("Options");
		gd.addMessage("Estimated # Frames = "+reader.nframes);
		gd.addMessage("Estimated Size (MB) = "+(float)(size/1000000.0));
		gd.addCheckbox("Read_Start_To_Finish",true);
		gd.addNumericField("Start_Frame",1,0);
		gd.addNumericField("End_Frame",reader.nframes,0);
		gd.addNumericField("Frame_Interval",10,0);
		gd.addCheckbox("Convert_to_greyscale?",true);
		gd.addCheckbox("Crop?",false);
		gd.addNumericField("Crop_X_Start",0,0);
		gd.addNumericField("Crop_Y_Start",0,0);
		gd.addNumericField("Crop_Width",reader.width,0);
		gd.addNumericField("Crop_Height",reader.height,0);
		gd.addCheckbox("Show_Crop_GUI",false);
		gd.addNumericField("Spatially_Bin_By?",1,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean readall=gd.getNextBoolean();
		int start=(int)gd.getNextNumber()-1;
		int end=(int)gd.getNextNumber();
		int interval=(int)gd.getNextNumber();
		boolean grey=gd.getNextBoolean();
		boolean crop=gd.getNextBoolean();
		int xstart=(int)gd.getNextNumber();
		int ystart=(int)gd.getNextNumber();
		int width=(int)gd.getNextNumber();
		int height=(int)gd.getNextNumber();
		boolean cropgui=gd.getNextBoolean();
		int bin=(int)gd.getNextNumber();
		if(readall) end=0;
		Rectangle r=new Rectangle(xstart,ystart,width,height);
		if(crop && cropgui){
			reader.firstframe.show();
			(new WaitForUserDialog("Select Crop Region with Rectangle and Click OK")).show();
			Roi roi=reader.firstframe.getRoi();
			if(roi!=null){
				r=roi.getBounds();
			}
			reader.firstframe.hide();
		}
		ImagePlus imp=reader.subopen(start,end,interval,grey,r,bin);
		imp.show();
		reader.close_xuggler();
	}

}
