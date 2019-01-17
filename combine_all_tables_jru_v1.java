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

public class combine_all_tables_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Add_Titles",true);
		gd.addCheckbox("Close Original",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean addtitles=gd.getNextBoolean();
		boolean closeorig=gd.getNextBoolean();
		//first get the table window
		Frame[] niframes=WindowManager.getNonImageWindows();
		boolean first=true;
		TextWindow tw2=null;
		int ncols=0;
		String[] col_labels=null;
		for(int i=0;i<niframes.length;i++){
			if(niframes[i] instanceof TextWindow && !niframes[i].getTitle().equals("Log") && !niframes[i].getTitle().equals("Debug")){
				TextWindow tw=(TextWindow)niframes[i];
				TextPanel tp=tw.getTextPanel();
				List<List<String>> listtable=table_tools.table2listtable(tp);
				if(listtable.size()==0){
					col_labels=table_tools.getcollabels(tp);
					ncols=col_labels.length;
					ArrayList<String> temp=new ArrayList<String>();
					for(int j=0;j<ncols;j++) temp.add("");
					listtable.add(temp);
				}
				if(first){
					col_labels=table_tools.getcollabels(tp);
					ncols=col_labels.length;
					String headings=tp.getColumnHeadings();
					if(addtitles){
						String[] titles=repeated(tw.getTitle(),listtable.size());
						table_tools.add_listtable_column(listtable,titles,0);
						headings="name\t"+headings;
					}
					tw2=new TextWindow("Combined Table",headings,table_tools.print_listtable(listtable),400,200);
					first=false;
				} else {
					if(addtitles){
						String[] titles=repeated(tw.getTitle(),listtable.size());
						table_tools.add_listtable_column(listtable,titles,0);
					}
					tw2.append(table_tools.print_listtable(listtable));
				}
				if(closeorig) tw.close();
			}
		}
	}

	public String[] repeated(String elem,int length){
		String[] temp=new String[length];
		for(int i=0;i<temp.length;i++) temp[i]=elem;
		return temp;
	}

}
