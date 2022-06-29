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

public class get_metadata_tag_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String dir = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		LOCI_file_reader r=new LOCI_file_reader();
		String[] names=r.getSeriesNames(dir,fname);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("No_of_Keys",1,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int nkeys=(int)gd.getNextNumber();
		String[] keys=new String[nkeys];
		GenericDialog gd2=new GenericDialog("Options");
		for(int i=0;i<nkeys;i++) gd2.addStringField("Key_"+(i+1),"");
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		for(int i=0;i<nkeys;i++) keys[i]=gd2.getNextString();
		String[][] vals=r.batch_get_series_metadata_value(dir,fname,keys);
		String[][] newvals=new String[names.length][nkeys+2];
		for(int i=0;i<names.length;i++){
			newvals[i][0]=names[i]; 
			newvals[i][1]=""+(i+1);
			for(int j=0;j<nkeys;j++) newvals[i][j+2]=vals[j][i];
		}
		String[] col_labels=new String[nkeys+2];
		col_labels[0]="name"; col_labels[1]="series"; for(int j=0;j<nkeys;j++) col_labels[j+2]=keys[j];
		TextWindow tw=jutils.selectTable("Metadata tags");
		if(tw==null) table_tools.create_table("Metadata tags",newvals,col_labels);
		else tw.append(table_tools.print_string_array(newvals));
		/*for(int i=0;i<names.length;i++){
			IJ.log(names[i]+" , "+vals[0][i]+" , "+vals[1][i]);
		}*/
		//Hashtable<String,Object> meta=r.getSeriesMetaData(dir,fname,3);
		//r.dumpMetaData(meta);
		/*String[] names=r.getSeriesNames(dir,fname);
		for(int i=0;i<names.length;i++){
			Hashtable<String,Object> meta=r.getSeriesMetaData(dir,fname,i);
			String xval=meta.get("X Location").toString();
			String yval=meta.get("Y Location").toString();
			IJ.log(""+names[i]+","+xval+","+yval);
			IJ.showProgress(i,names.length);
		}*/
		/*String[] names=r.getSeriesNames(dir,fname);
		for(int i=0;i<names.length;i++){
			String xval=r.get_metadata_value(dir, fname, i, "X Location",false);
			String yval=r.get_metadata_value(dir, fname, i, "Y Location",false);
			IJ.log(""+names[i]+","+xval+","+yval);
			IJ.showProgress(i,names.length);
		}*/
		//String meta=r.get_metadata_value(dir, fname, 3, "X Location",false);
		//Hashtable<String,Object> meta=r.getMetaDataValue(dir,fname);
		//String temp="2 (raw tile 187) X Location";
		//String temp2=meta.get(temp).toString();
		//IJ.log(meta);
		//String temp="[Acquisition Parameters Common] Number Of Point";
		//int points=Integer.parseInt(meta.get(temp).toString());
		//temp="[Acquisition Parameters Common] Number Of Cycle";
	}

}
