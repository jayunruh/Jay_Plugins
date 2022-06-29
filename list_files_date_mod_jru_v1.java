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
import ij.text.*;

public class list_files_date_mod_jru_v1 implements PlugIn {
	TextWindow tw;

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image Sequence...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Sort_Names_Numerically",false);
		gd.addCheckbox("Sort_Date_Modified",false);
		gd.addStringField("File_Name_Contains:","",10);
		gd.addCheckbox("Full_Path?",true);
		gd.addCheckbox("Include_Subfolders",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean sortnum=gd.getNextBoolean();
		boolean sortdatemod=gd.getNextBoolean();
		String filter=gd.getNextString();
		boolean path=gd.getNextBoolean();
		boolean subfolders=gd.getNextBoolean();
		int period=fname.lastIndexOf('.');
		String extension=fname.substring(period+1);
		String[] mask=null;
		tw=new TextWindow("File_Names","FileName","",200,400);
		if(filter!="" && filter!=null && filter.length()>0){
			mask=new String[2]; mask[1]=filter; mask[0]=extension;
			//IJ.log("masks: "+mask[0]+" , "+mask[1]);
		} else {
			mask=new String[1]; mask[0]=extension;
		}
		if(!subfolders){
			String[] flist=null;
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
				if(!path) tw.append(""+flist[i]);
				else tw.append(directory+flist[i]);
			}
		} else {
			listDir(directory,mask,sortnum,sortdatemod,path);
		}
	}

	public void listDir(String dirpath,String[] mask,boolean sortnum,boolean sortdatemod,boolean path){
		String[] flist=null;
		if(sortnum){
			flist=(new jdataio()).get_numeric_sorted_string_list(dirpath,mask,null);
		} else {
			if(sortdatemod){
				flist=(new jdataio()).get_datemod_sorted_string_list(dirpath,mask,null);
			} else {
				flist=(new jdataio()).get_sorted_string_list(dirpath,mask,null);
			}
		}
		for(int i=0;i<flist.length;i++){
			if(!path) tw.append(""+flist[i]);
			else tw.append(dirpath+flist[i]);
			if((new File(dirpath+flist[i])).isDirectory()){
				listDir(dirpath+flist[i]+File.separator,mask,sortnum,sortdatemod,path);
			}
		}
		//now list the directories
		String[] flist2=(new File(dirpath)).list();
		for(int i=0;i<flist2.length;i++){
			if((new File(dirpath+flist2[i])).isDirectory()){
				listDir(dirpath+flist2[i]+File.separator,mask,sortnum,sortdatemod,path);
			}
		}
		return;
	}

}
