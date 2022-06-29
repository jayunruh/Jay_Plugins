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
import java.awt.Frame;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.io.*;
import jguis.*;
import jalgs.*;
import java.io.*;
import java.util.*;

public class batch_file_open_jru_v1 implements PlugIn {
	String directory,extension;
	String[] flist;
	int nfiles;

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image Sequence...", arg);
        		directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Sort_Names_Numerically",false);
		gd.addCheckbox("Sort_Date_Modified",false);
		gd.addStringField("File_Name_Contains:","",10);
		gd.addCheckbox("Use_ImageJ_Opener",false);
		//gd.addNumericField("Scale Factor",1.0,5,15,null);
		//gd.addCheckbox("Force Z Stack",true);
		//gd.addCheckbox("Max_Project",false);
		//gd.addCheckbox("Open_Virtual_After",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean sortnum=gd.getNextBoolean();
		boolean sortdatemod=gd.getNextBoolean();
		String filter=gd.getNextString();
		boolean ijopen=gd.getNextBoolean();
		//scale=gd.getNextNumber();
		//boolean forcez=gd.getNextBoolean();
		//maxproj=gd.getNextBoolean();
		//boolean openafter=gd.getNextBoolean();
		int period=fname.lastIndexOf('.');
		extension=fname.substring(period+1);
		String[] mask=null;
		if(filter!="" && filter!=null && filter.length()>0){
			mask=new String[2]; mask[1]=filter; mask[0]=extension;
			IJ.log("masks: "+mask[0]+" , "+mask[1]);
		} else {
			mask=new String[1]; mask[0]=extension;
		}
		if(sortnum){
			flist=(new jdataio()).get_numeric_sorted_string_list(directory,mask,null);
		} else {
			if(sortdatemod){
				flist=(new jdataio()).get_datemod_sorted_string_list(directory,mask,null);
			} else {
				flist=(new jdataio()).get_sorted_string_list(directory,mask,null);
			}
		}
		for(int i=0;i<flist.length;i++){
			IJ.log(""+flist[i]);
		}
		if(flist==null){return;}
		nfiles=flist.length;
		//IJ.showMessage("test");
		Opener opener=new Opener(); //opener.setSilentMode(true);
		for(int i=0;i<flist.length;i++){
			if(extension.equals("xls")){
				List<List<String>> table=table_tools.getTableFromFile(new File(directory+flist[i]),"\t",false);
				String[] headings=table_tools.list2stringarray(table.get(0));
				table.remove(0);
				table_tools.create_table(flist[i],table,headings);
			} else {
				ImagePlus currimp=null;
				if(!extension.equals("lsm")){
					if(ijopen) currimp=opener.openImage(directory+flist[i]);
					else currimp=(new LOCI_file_reader()).get_loci_imp(directory,flist[i]);
				} else {
					currimp=(new LSM_file_reader()).open(directory,flist[i]);
				}
				if(currimp!=null) currimp.show();
			}
			if(IJ.escapePressed()) break;
		}
	}

}
