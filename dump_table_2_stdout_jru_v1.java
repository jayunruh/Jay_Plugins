/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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

public class dump_table_2_stdout_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		GenericDialog gd=new GenericDialog("Options");
		String[] delims={"tab","comma","space","newline"};
		gd.addChoice("Delimiter",delims,delims[0]);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int delim=gd.getNextChoiceIndex();
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		//System.out.println("Printing "+tw[0].getTitle());
		System.out.println(table_tools.print_string_array(col_labels,delim));
		for(int i=0;i<listtable.size();i++){
			System.out.println(table_tools.print_string_array(listtable.get(i),delim));
		}
	}

}
