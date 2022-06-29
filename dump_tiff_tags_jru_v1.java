/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.io.*;
import ij.plugin.*;
import jguis.*;
import jalgs.*;
import loci.formats.gui.XMLWindow;
import ij.text.*;

public class dump_tiff_tags_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open tiff File...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		IJ.log(fname);
		TiffReader tr=new TiffReader(directory+fname);
		List<List<String>> tags=tr.readTags();
		tr.closeStream();
		String[] collabels={"tag","data_type","count","value","tag_name"};
		for(int i=0;i<tags.size();i++){
			String val=tags.get(i).get(3);
			String tag=tags.get(i).get(0);
			if(val.length()>50){
				if(val.indexOf("xml")>=0){
					try{
						XMLWindow xmlWindow=new XMLWindow("tag = "+tag);
						xmlWindow.setXML(val.substring(0,val.length()));
						xmlWindow.setVisible(true);
					} catch(Exception e){
						IJ.log((new jdataio()).getExceptionTrace(e));
						new TextWindow("tag = "+tag,val.substring(0,val.length()),400,200);
					}
				}
				val=val.substring(0,50);
				val=val.replace('\n',' ');
				val=val.replace('\t',' ');
				tags.get(i).set(3,val);
			}
		}
		table_tools.create_table(fname+"_tags",tags,collabels);
	}

}
