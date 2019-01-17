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
import ij.text.*;
import java.util.*;
import jguis.*;

public class delete_table_column_jru_v1 implements PlugIn {

	public void run(String arg) {
		//first get the table window
		TextWindow[] tw=jutils.selectTables(false,1);
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Column_to_Delete",col_labels,col_labels[0]);
		gd.addCheckbox("Replace Original Table",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int col=gd.getNextChoiceIndex();
		boolean replace=gd.getNextBoolean();
		String headings=tp.getColumnHeadings();
		String newheadings=table_tools.delete_col_label(headings,col);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		table_tools.delete_listtable_column(listtable,col);
		if(!replace) new TextWindow(tw[0].getTitle()+"_1",newheadings,table_tools.print_listtable(listtable),400,200);
		else table_tools.replace_table(tp,listtable,newheadings);
	}
}
