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
import java.io.*;
import jalgs.*;
import jguis.*;
import ij.text.*;

public class import_flowcyte_jru_v1 implements PlugIn {

	public void run(String arg) {
		//start by reading the header
		//jdataio jdio=new jdataio();
		OpenDialog od = new OpenDialog("Open File",arg);
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Show Metadata?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean showmetadata=gd.getNextBoolean();
		//Object[] data=getFCSFile(directory+name);
		Object[] data=(new import_flowcyte()).getFCSFile(directory+name);
		if(data!=null){
			if(showmetadata){
				String[][] meta=(String[][])data[2];
				for(int i=0;i<meta[0].length;i++){
					IJ.log(meta[0][i]+" , "+meta[1][i]);
				}
			}
			String[] chnames=(String[])data[0];
			String headings=table_tools.print_string_array(chnames);
			TextWindow tw=new TextWindow(name,headings,"",400,200);
			if(data[1]!=null) tw.append(table_tools.print_float_array((float[][])data[1]));
		}
	}

}
