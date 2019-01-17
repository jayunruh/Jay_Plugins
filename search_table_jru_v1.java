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
import jguis.*;
import jalgs.*;
import ij.text.*;
import java.util.*;

public class search_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		//first get the table window
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length];
		for(int i=0;i<niframes.length;i++){
			titles[i]=niframes[i].getTitle();
		}
		GenericDialog gd=new GenericDialog("Name List");
		gd.addTextAreas("",null,10,20);
		gd.addChoice("Table Window",titles,titles[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String input=gd.getNextText();
		int index=gd.getNextChoiceIndex();
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String[] col_labels=table_tools.getcollabels(tp);
			GenericDialog gd2=new GenericDialog("Choose Column");
			gd2.addChoice("Search_column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int colindex=gd2.getNextChoiceIndex();
			List<List<String>> table=table_tools.table2listtable(tp);
			String[] list=(new delimit_string('\t')).getrows(input);
			for(int i=0;i<list.length;i++){
				list[i]=list[i].trim();
				//IJ.log(""+list[i]);
			}
			table_tools.sort_listtable(table,colindex);
			List<List<String>> results=new ArrayList<List<String>>();
			for(int j=0;j<list.length;j++){
				int foundindex=table_tools.find_sorted_listtable_string(table,colindex,list[j]);
				if(foundindex<0 || foundindex>=table.size()){
					IJ.log(list[j]+" not found");
				} else {
					results.add(table.get(foundindex));
				}
			}
			new TextWindow("Search_Results",tp.getColumnHeadings(),table_tools.print_listtable(results),400,200);
		} else {
			IJ.showMessage("wrong window type");
		}
		
	}

}
