/******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
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
import java.io.*;
import ij.plugin.*;
import ij.text.*;
import jguis.*;
import jalgs.*;
import ij.io.*;

public class export_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		if(tw==null) return;
		GenericDialog gd=new GenericDialog("Options");
		String[] formats={"xls(tab)","txt(tab)","txt(space)","txt(comma)","csv","other(tab)","other(space)","other(comma)"};
		String[] exts={".xls",".txt",".txt",".txt",".csv","","",""};
		gd.addChoice("Format",formats,formats[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int format=gd.getNextChoiceIndex();
		SaveDialog sd = new SaveDialog("Save Text...", tw[0].getTitle(), exts[format]);
       		String fileName = sd.getFileName();
		if (fileName == null) return;
		String fileDir = sd.getDirectory();
		int delim=0;
		if(format==2 || format==6) delim=2;
		if(format==3 || format==4 || format==7) delim=1;
		try{
			TextPanel tp=tw[0].getTextPanel();
			List<List<String>> listtable=table_tools.table2listtable(tp);
			String[] collabels=table_tools.getcollabels(tp);
			BufferedWriter b=new BufferedWriter(new FileWriter(new File(fileDir+fileName)));
			//IJ.log(table_tools.print_string_array(collabels,delim));
			b.write(table_tools.print_string_array(collabels,delim)+"\n");
			String s=table_tools.print_listtable(listtable,delim);
			b.write(s);
			b.close();
			//(new jdataio()).writestringfile(b,s);
		} catch(IOException e){
			IJ.error("error writing file");
		}
	}

}
