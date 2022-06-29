/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import java.awt.Polygon;
import ij.plugin.*;
import ij.io.*;
import jguis.*;
import java.util.*;
import jalgs.*;
import ij.text.*;

public class get_loci_objective_info_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String dir = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		LOCI_file_reader r=new LOCI_file_reader();
		//String[] names=r.getSeriesNames(dir,fname);
		String[][] vals=r.getOMEXMLObjectiveInfo(dir,fname,0);
		String[] col_labels={"Description","Value"};
		table_tools.create_table("Objective Info",vals,col_labels);
	}

}
