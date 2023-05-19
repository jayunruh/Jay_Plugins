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
import jalgs.*;

public class LOCI_file_reader_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		int nseries=(new LOCI_file_reader()).getNSeries(directory,fname);
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Avoid_Metadata",false);
		gd2.addCheckbox("Simple_Import",false);
		gd2.addCheckbox("Large_Image_Canvas (lower resolution)",false);
		gd2.addNumericField("Max_Display_Pixels (for large canvas)",1024,0);
		gd2.addCheckbox("Dynamic_Resolution (for large canvas)",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		boolean nometa=gd2.getNextBoolean();
		boolean simple=gd2.getNextBoolean();
		boolean large=gd2.getNextBoolean();
		int maxpix=(int)gd2.getNextNumber();
		boolean dynamicscale=gd2.getNextBoolean();
		LOCI_file_reader lfr=new LOCI_file_reader();
		lfr.nometa=nometa;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Output Metadata?",false);
		if(nseries>1){
			String[] seriesnames=lfr.getSeriesNames(directory,fname);
			if(seriesnames==null){
				seriesnames=new String[nseries];
				for(int i=0;i<nseries;i++) seriesnames[i]=""+i;
			}
			gd.addChoice("Choose_Series",seriesnames,seriesnames[0]);
		}
		gd.addCheckbox("Z_Project?",false);
		gd.addChoice("Proj_Stat",jstatistics.stats,jstatistics.stats[0]);
		gd.addNumericField("Ref_Chan (0 if none)",0,0);
		gd.addCheckbox("Open_all_frames",true);
		gd.addNumericField("Frames_to_open",1,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean outmeta=gd.getNextBoolean();
		int series=0;
		if(nseries>1) series=gd.getNextChoiceIndex();
		boolean zproj=gd.getNextBoolean();
		int refchan=(int)gd.getNextNumber()-1;
		boolean allframes=gd.getNextBoolean();
		int nframes=(int)gd.getNextNumber();
		int[] lims={0,-1,0,-1,0,nframes-1};
		if(allframes) lims[5]=-1;
		String projstat=jstatistics.stats[gd.getNextChoiceIndex()];
		ImagePlus imp=null;
		if(!simple) imp=lfr.get_loci_subimp(directory,fname,outmeta,series,zproj,projstat,refchan,lims);
		else imp=lfr.get_loci_imp_simple(directory,fname,series);
		//ImagePlus imp=(new LOCI_file_reader()).get_loci_imp(directory,fname,outmeta,series);
		if(!large){
			imp.show();
		} else {
			largeImageCanvas canvas=new largeImageCanvas(imp);
			canvas.maxpixels=maxpix;
			canvas.dynamicscale=dynamicscale;
			ImageWindow win=new StackWindow(imp,canvas);
		}
		IJ.log(""+lims[0]+" , "+lims[1]+" , "+lims[2]+" , "+lims[3]+" , "+lims[4]+" , "+lims[5]);
	}

}
