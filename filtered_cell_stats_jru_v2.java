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

public class filtered_cell_stats_jru_v2 implements PlugIn {
	//this version doesn't use the filter selection table (only one filter allowed)

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
		float[] statoptions=jutils.getStatsOptions(stat);
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String[] col_labels=table_tools.getcollabels(tp);
			int ncols=col_labels.length;
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Cell Name Column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int cellnameindex=(int)gd2.getNextChoiceIndex();
			float[][] filters=get_filters(col_labels,cellnameindex,0);
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
			List<List<String>> stattable=table_tools.get_cell_stat_list(filtered,cellnameindex,stat,cellctcol,statoptions);
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

	public float[][] get_filters(String[] col_labels,int cellnameindex,int recursion){
		GenericDialog gd=new GenericDialog("Select Filter");
		if(recursion==0){
			gd.addCheckbox("Use Filter?",false);
			String[] labels2=table_tools.make_labels_unique(col_labels);
			labels2[cellnameindex]="null";
			gd.addChoice("Filter Parameter",labels2,labels2[0]);
			String[] filters={"<","=",">"};
			gd.addChoice("Filter_Type",filters,filters[0]);
			gd.addNumericField("Filter_Value",0.0,5,15,null);
			gd.addCheckbox("Add_Another Filter?",false);
		} else {
			gd.addCheckbox("Use"+recursion+" Filter?",false);
			String[] labels2=table_tools.make_labels_unique(col_labels);
			labels2[cellnameindex]="null";
			gd.addChoice("Filter"+recursion+" Parameter",labels2,labels2[0]);
			String[] filters={"<","=",">"};
			gd.addChoice("Filter_Type"+recursion,filters,filters[0]);
			gd.addNumericField("Filter_Value"+recursion,0.0,5,15,null);
			gd.addCheckbox("Add_Another"+recursion+" Filter?",false);
		}
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		if(!gd.getNextBoolean()) return null;
		int colindex=gd.getNextChoiceIndex();
		if(colindex==cellnameindex) return null;
		float[][] outvals=new float[1][3];
		outvals[0][0]=colindex;
		outvals[0][1]=gd.getNextChoiceIndex();
		outvals[0][2]=(float)gd.getNextNumber();
		boolean add_another=gd.getNextBoolean();
		if(add_another && recursion<6){
			float[][] temp=get_filters(col_labels,cellnameindex,recursion+1);
			outvals=appendFilters(outvals,temp);
		}
		return outvals;
	}

	public float[][] appendFilters(float[][] filters1,float[][] filters2){
		float[][] temp=new float[filters1.length+filters2.length][];
		for(int i=0;i<filters1.length;i++) temp[i]=filters1[i];
		for(int i=0;i<filters2.length;i++) temp[i+filters1.length]=filters2[i];
		return temp;
	}

}
