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
import java.util.*;
import java.awt.Frame;
import ij.plugin.*;
import ij.text.*;
import java.io.*;
import ij.io.*;
import jguis.*;

public class batch_import_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		DirectoryChooser od=new DirectoryChooser("Choose Directory");
		String dir=od.getDirectory();
		if(dir==null) return;
		measure_directory(dir);
	}

	public void measure_directory(String dir){
		String[] list=new File(dir).list();
		int nfiles=list.length;
		for(int i=0;i<nfiles;i++){
			if(list[i].endsWith(".xls")){
				List<List<String>> listtable=(new import_table_jru_v1()).getTableFromFile(dir,list[i],"\t",false);
				if(listtable!=null){
					table_tools.create_table(list[i],listtable);
				}
			} else {
				File temp=new File(dir+list[i]);
				if(temp.isDirectory()){
					measure_directory(temp.getAbsolutePath()+File.separator);
				}
			}
		}
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
