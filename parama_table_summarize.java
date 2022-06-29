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
import jguis.*;
import java.io.*;
import ij.text.*;
import java.util.*;

public class parama_table_summarize implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		if(tw==null) return;
		//start by getting the listtable
		TextPanel tp=tw[0].getTextPanel();
		List<List<String>> listtable=table_tools.table2listtable(tp);
		String[] col_labels=table_tools.getcollabels(tp);
		int fileindex=find_string(col_labels,"file");
		//now sort the listtable
		table_tools.sort_listtable(listtable,fileindex);
		//assume we already have corrected ratio and row, col, and wellid columns
		//now count the cell stats for all nuclei
		List<List<String>> nuclear_counts=table_tools.get_cell_stat_list(listtable,fileindex,"Count");
		List<List<String>> avgs=table_tools.get_cell_stat_list(listtable,fileindex,"Avg");
		//and for the high nuclei
		int ratio2index=find_string(col_labels,"ratio2");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Ratio_Threshold",1.2,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float thresh=(float)gd.getNextNumber();
		float[][] filter={{ratio2index,2,thresh}};
		List<List<String>> filtered=table_tools.get_filtered_listtable(listtable,filter);
		List<List<String>> high_nuc_counts=table_tools.get_cell_stat_list(filtered,fileindex,"Count");
		List<String> celllist=table_tools.get_cell_list(listtable,fileindex);
		//check if any images are unlisted after filtering
		List<List<String>> extralist=table_tools.get_missing_cells(high_nuc_counts,fileindex,celllist,"0");
		for(int i=0;i<extralist.size();i++) high_nuc_counts.add(extralist.get(i));
		table_tools.sort_listtable(high_nuc_counts,fileindex);
		//find the percent high as well as row,col
		int rowindex=find_string(col_labels,"row");
		int colindex=find_string(col_labels,"col");
		int idindex=find_string(col_labels,"id");
		List<List<String>> finallist=new ArrayList<List<String>>();
		for(int i=0;i<nuclear_counts.size();i++){
			int nnucs=(int)table_tools.get_number(nuclear_counts,i,idindex);
			int nhigh=(int)table_tools.get_number(high_nuc_counts,i,idindex);
			float frachigh=(float)nhigh/(float)nnucs;
			int row=(int)table_tools.get_number(avgs,i,rowindex);
			int col=(int)table_tools.get_number(avgs,i,colindex);
			String file=nuclear_counts.get(i).get(fileindex);
			String wellid=padinteger(row)+padinteger(col);
			List<String> temp=new ArrayList<String>();
			temp.add(file);
			temp.add(wellid);
			temp.add(""+row);
			temp.add(""+col);
			temp.add(""+frachigh);
			finallist.add(temp);
		}
		String headings="file\twellid\trow\tcol\tfhigh";
		//output the processed table for each image
		new TextWindow("Image Table",headings,table_tools.print_listtable(finallist),400,200);
		//now average the well replicates
		List<List<String>> wellavg=table_tools.get_cell_stat_list(finallist,1,"Avg");
		List<List<String>> wellerr=table_tools.get_cell_stat_list(finallist,1,"StDev");
		List<List<String>> welltable=new ArrayList<List<String>>();
		String headings2="wellid\trow\tcol\tfhigh\tstdev";
		for(int i=0;i<wellavg.size();i++){
			List<String> temp=new ArrayList<String>();
			temp.add(wellavg.get(i).get(1));
			temp.add(wellavg.get(i).get(2));
			temp.add(wellavg.get(i).get(3));
			temp.add(wellavg.get(i).get(4));
			temp.add(wellerr.get(i).get(4));
			welltable.add(temp);
		}
		new TextWindow("Well Table",headings2,table_tools.print_listtable(welltable),400,200);
	}

	public String padinteger(int val){
		String temp=""+val;
		if(val<100) temp="0"+temp;
		if(val<10) temp="0"+temp;
		return temp;
	}

	public int find_string(String[] arr,String searchstring){
		for(int i=0;i<arr.length;i++){
			if(arr[i].equalsIgnoreCase(searchstring)) return i;
		}
		return -1;
	}

}
