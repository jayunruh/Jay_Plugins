/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import java.util.*;
import ij.plugin.*;
import ij.io.*;
import java.io.*;
import jalgs.*;
import ij.text.*;
import jguis.*;

public class copy_table_files_jru_v1 implements PlugIn {
	int updatecount=0;

	public void run(String arg) {
		//this plugin copies files from a table from one director into another
		DirectoryChooser dc=new DirectoryChooser("Select_Source Folder");
		String src=dc.getDirectory();
		if(src==null) return;
		DirectoryChooser dc2=new DirectoryChooser("Select_Destination Folder");
		String dest=dc2.getDirectory();
		if(dest==null) return;
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"File_Table"});
		if(tw==null) return;
		TextPanel tp=tw[0].getTextPanel();
		List<List<String>> listtable=table_tools.table2listtable(tp);
		String[] col_labels=table_tools.getcollabels(tp);
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Name Column",col_labels,col_labels[0]);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		int col=gd2.getNextChoiceIndex();
		for(int i=0;i<listtable.size();i++){
			String fname=listtable.get(i).get(col);
			byte[] contents=getFile(src+fname);
			if(contents==null){IJ.log("error reading "+fname); continue;}
			if(!writeFile(contents,dest+fname)){IJ.log("error writing "+fname); continue;}
			IJ.log(fname);
		}
	}

	public boolean writeFile(byte[] contents,String fname){
		return (new jdataio()).writeFile(contents,fname);
	}

	public byte[] getFile(String fname){
		return (new jdataio()).readentirebytefile(new File(fname));
	}

}
