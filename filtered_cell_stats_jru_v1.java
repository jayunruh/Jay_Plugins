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
import ij.util.*;
import ij.text.*;
import jalgs.*;
import jguis.*;

public class filtered_cell_stats_jru_v1 implements PlugIn {

	public void run(String arg) {
		//first get the table window
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length];
		for(int i=0;i<niframes.length;i++){
			titles[i]=niframes[i].getTitle();
		}
		GenericDialog gd=new GenericDialog("Windows");
		gd.addChoice("Windows",titles,titles[0]);
		String[] stats=jstatistics.stats;
		gd.addChoice("Statistic?",stats,stats[0]);
		gd.addStringField("String for Empty Cells","NaN");
		gd.addCheckbox("Add_Cell_Count_Column",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		int statindex=gd.getNextChoiceIndex();
		String emptystring=gd.getNextString();
		boolean cellctcol=gd.getNextBoolean();
		String stat=stats[statindex];
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String[] col_labels=table_tools.getcollabels(tp);
			int ncols=col_labels.length;
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Cell Name Column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int cellnameindex=(int)gd2.getNextChoiceIndex();
			float[][] filters=get_filters(col_labels,cellnameindex);
			List<List<String>> listtable=table_tools.table2listtable(tp);
			table_tools.sort_listtable(listtable,cellnameindex);
			//new TextWindow("Sorted Table",tp.getColumnHeadings(),table_tools.print_listtable(listtable),400,200);
			List<List<String>> filtered=null;
			if(filters!=null){
				filtered=table_tools.get_filtered_listtable(listtable,filters);
				if(filtered==null){
					return;
				}
			} else {
				filtered=listtable;
			}
			//new TextWindow("Filtered Table",tp.getColumnHeadings(),table_tools.print_listtable(filtered),400,200);
			List<List<String>> stattable=table_tools.get_cell_stat_list(filtered,cellnameindex,stat,cellctcol);
			//now we need to find the cells that were completely filtered and return NaN for those
			//both tables are already sorted
			if(filters!=null){
				List<String> celllist=table_tools.get_cell_list(listtable,cellnameindex);
				List<List<String>> extralist=new ArrayList<List<String>>();
				for(int i=0;i<celllist.size();i++){
					String cellname=celllist.get(i);
					int foundindex=table_tools.find_sorted_listtable_string(stattable,cellnameindex,cellname);
					if(foundindex<0){
						List<String> temp=new ArrayList<String>();
						for(int j=0;j<ncols;j++){
							if(j==cellnameindex) temp.add(cellname);
							else temp.add(emptystring);
						}
						if(cellctcol) temp.add("0");
						extralist.add(temp);
						//stattable.add(temp);
					}
				}
				for(int i=0;i<extralist.size();i++) stattable.add(extralist.get(i));
				table_tools.sort_listtable(stattable,cellnameindex);
			}
			String newcollabels=tp.getColumnHeadings();
			if(cellctcol) newcollabels+="\tcellct";
			TextWindow tw2=new TextWindow("Table Data Stats",newcollabels,table_tools.print_listtable(stattable),400,200);	
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public float[][] get_filters(String[] col_labels,int cellnameindex){
		String[] columnlabels={"Use Filter?","Filter Parameter","Filter Type","Filter Value"};
		String[] filterlabels=new String[col_labels.length-1];
		for(int i=0;i<cellnameindex;i++){
			filterlabels[i]=col_labels[i];
		}
		for(int i=(cellnameindex+1);i<col_labels.length;i++){
			filterlabels[i-1]=col_labels[i];
		}
		String[] filters={"<","=",">"};
		Object[][] tabledata={{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"}};
		int[][] options=new int[8][4];
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Filters",columnlabels,tabledata,options);
		if(retvals==null){return null;}
		int nfilters=0;
		boolean[] usefilter=new boolean[8];
		for(int i=0;i<8;i++){
			usefilter[i]=((Boolean)retvals[i][0]).booleanValue();
			if(usefilter[i]){nfilters++;}
		}
		if(nfilters==0){return null;}
		float[][] outvals=new float[nfilters][3];
		int counter=0;
		for(int i=0;i<8;i++){
			if(usefilter[i]){
				outvals[counter][0]=(float)table_tools.find_string_array_index(col_labels,(String)retvals[i][1]);
				if(outvals[counter][0]<0.0f){return null;}
				outvals[counter][1]=(float)table_tools.find_string_array_index(filters,(String)retvals[i][2]);
				if(outvals[counter][1]<0.0f){return null;}
				try{
					outvals[counter][2]=Float.parseFloat((String)retvals[i][3]);
				} catch(NumberFormatException e){
					return null;
				}
				counter++;
			}
		}
		return outvals;
	}

}
