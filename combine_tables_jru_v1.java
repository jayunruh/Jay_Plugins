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

public class combine_tables_jru_v1 implements PlugIn {

	public void run(String arg) {
		Object[] windowlist=jutils.getTableWindowList(false);
		String[] titles=(String[])windowlist[1];
		Frame[] frames=(Frame[])windowlist[0];
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Combine_Vertically",false);
		gd.addChoice("Table1",titles,titles[0]);
		gd.addChoice("Table2",titles,titles[1]);
		gd.addCheckbox("Replace Table1",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean vertical=gd.getNextBoolean();
		TextWindow tw1=(TextWindow)frames[gd.getNextChoiceIndex()];
		TextWindow tw2=(TextWindow)frames[gd.getNextChoiceIndex()];
		boolean replace=gd.getNextBoolean();
		TextPanel tp1=tw1.getTextPanel();
		TextPanel tp2=tw2.getTextPanel();
		List<List<String>> listtable1=table_tools.table2listtable(tp1);
		List<List<String>> listtable2=table_tools.table2listtable(tp2);
		if(vertical){
			String headings=tp1.getColumnHeadings();
			if(!replace){
				TextWindow tw3=new TextWindow("Combined Table",headings,table_tools.print_listtable(listtable1),400,200);
				tw3.append(table_tools.print_listtable(listtable2));
			} else {
				tp1.append(table_tools.print_listtable(listtable2));
			}
		} else {
			String headings=tp1.getColumnHeadings()+"\t"+tp2.getColumnHeadings();
			for(int i=0;i<listtable2.get(0).size();i++){
				String[] col=table_tools.get_listtable_column(listtable2,i);
				table_tools.add_listtable_column(listtable1,col,1000);
			}
			if(replace) table_tools.replace_table(tp1,listtable1,headings);
			else new TextWindow("Combined Table",headings,table_tools.print_listtable(listtable1),400,200);
		}
	}

}
