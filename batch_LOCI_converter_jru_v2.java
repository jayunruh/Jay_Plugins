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

public class batch_LOCI_converter_jru_v2 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		DirectoryChooser dc=new DirectoryChooser("Save Directory");
		String outdir=dc.getDirectory();
		if(outdir==null) return;
		LOCI_file_reader lfr=new LOCI_file_reader();
		boolean success=lfr.batchExportSeries(directory,fname,outdir,false);
		if(!success) IJ.log("something went terribly wrong");
	}

}
