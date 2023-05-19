/*******************************************************************************
 * Copyright (c) 2021 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import ij.text.*;

public class view_screen_folder_jru_v1 implements PlugIn {
	//this plugin starts an interactive viewer for plate screening data
	//data is initialized from a folder and optionally a table with layout information
	//the table should have columns: name,row,col,field,frame,slice,channel
	//images should be single planes

	public void run(String arg) {
		//start by getting the plate directory
		DirectoryChooser dc=new DirectoryChooser("Select Plate Folder");
		String dir=dc.getDirectory();
		if(dir==null) return;
		plate_viewer_panel pvp=new plate_viewer_panel();
		//now optionally get the layout table (otherwise initialize from folder)
		List<List<String>> flist=null;
		TextWindow[] tw=jutils.selectTables(true,1);
		if(tw[0]==null){
			flist=plate_viewer_panel.makePEFlist(dir);
			String[] col_labels=null;
			table_tools.create_table("file table",flist,col_labels);
		} else {
			flist=table_tools.table2listtable(tw[0].getTextPanel());
		}
		pvp.init(dir,flist);
		plate_viewer_panel.launch_frame(pvp);
	}

}
