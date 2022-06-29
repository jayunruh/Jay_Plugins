/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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

public class make_subtable_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("no_columns",1,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int ncols=(int)gd.getNextNumber();
		GenericDialog gd2=new GenericDialog("Select Columns");
		for(int i=0;i<ncols;i++) gd2.addChoice("Col"+(i+1),col_labels,col_labels[0]);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		int[] selcols=new int[ncols];
		for(int i=0;i<ncols;i++) selcols[i]=gd2.getNextChoiceIndex();
		String[] startcol=table_tools.get_listtable_column(listtable,selcols[0]);
		List<List<String>> subtable=table_tools.table2listtable(startcol,"\t");
		String[] sublabels=new String[ncols];
		sublabels[0]=col_labels[selcols[0]];
		for(int i=1;i<ncols;i++){
			String[] tempcol=table_tools.get_listtable_column(listtable,selcols[i]);
			table_tools.add_listtable_column(subtable,tempcol,ncols);
			sublabels[i]=col_labels[selcols[i]];
		}
		table_tools.create_table("Subtable",subtable,sublabels);
	}

}
