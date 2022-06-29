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
import jalgs.*;

public class import_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		/*OpenDialog od = new OpenDialog("Open File","",".txt");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null || name.equals("")) return;*/
		String args2=Macro.getOptions();
		File[] filearray=null;
		if(args2!=null && args2.length()>0){
			if(args2.indexOf("open=")>=0){
				filearray=new File[1];
				filearray[0]=new File(Macro.getValue(args2,"open",""));
			}
		}
		if(filearray==null) filearray=(new jdataio()).openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
		if(filearray.length==0) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Has_Column_Labels",true);
		String[] delims={"Tab","Comma","Space"};
		String[] delims2={"\t",","," "};
		gd.addChoice("Delimiter",delims,delims[0]);
		gd.addNumericField("ASCII number(optional)",0,0);
		gd.addCheckbox("Treat Consecutive Delims as 1",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean haslabels=gd.getNextBoolean();
		int delimindex=gd.getNextChoiceIndex();
		int ascii=(int)gd.getNextNumber();
		boolean noconsec=gd.getNextBoolean();
		String delim;
		if(ascii!=0){delim=""+(char)ascii;}
		else {
			delim=delims2[delimindex];
		}
		for(int i=0;i<filearray.length;i++){
			List<List<String>> listtable=getTableFromFile(filearray[i],delim,noconsec);
			if(listtable!=null){
				String[] headings=null;
				if(haslabels){
					headings=table_tools.list2stringarray(listtable.get(0));
					listtable.remove(0);
				} else {
					String temp=table_tools.createcollabels(listtable.get(0).size());
					headings=table_tools.split_string_tab(temp);
				}
				table_tools.create_table(filearray[i].getName(),listtable,headings);
			}
		}
	}

	public List<List<String>> getTableFromFile(File infile,String delim,boolean noconsec){
		try{
			BufferedReader b=new BufferedReader(new FileReader(infile));
			//String tempheadings=b.readLine();
			//String[] headings=table_tools.split(tempheadings,delim);
			List<String> lines=new ArrayList<String>();
			String temp=b.readLine();
			while(temp!=null && temp.length()>1){
				lines.add(temp);
				temp=b.readLine();
			}
			//(new WaitForUserDialog("Check Memory")).show();
			List<List<String>> listtable=table_tools.table2listtable(lines,delim,noconsec);
			lines=null;
			//(new WaitForUserDialog("Check Memory")).show();
			//table_tools.create_table(name,listtable,headings);
			return listtable;
		}
		catch(IOException e){
			showErrorMessage(e);
			return null;
			//return;
		}
	}

	public List<List<String>> getTableFromFile(String directory,String name,String delim,boolean noconsec){
		return getTableFromFile(new File(directory+name),delim,noconsec);
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
