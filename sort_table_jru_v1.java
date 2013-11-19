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
import java.util.*;
import ij.plugin.*;
import ij.util.*;
import ij.text.*;
import jguis.*;

public class sort_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we take a table column and convert it to a plot
		//first get the table window
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length];
		for(int i=0;i<niframes.length;i++){
			titles[i]=niframes[i].getTitle();
		}
		GenericDialog gd=new GenericDialog("Windows");
		gd.addChoice("Windows",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String[] col_labels=table_tools.getcollabels(tp);
			GenericDialog gd2=new GenericDialog("Choose Column");
			gd2.addChoice("Sort Column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int colindex=gd2.getNextChoiceIndex();
			List<List<String>> list=table_tools.table2listtable(tp);
			table_tools.sort_listtable(list,colindex);
			TextWindow tw2=new TextWindow("Sorted Table",tp.getColumnHeadings(),table_tools.print_listtable(list),400,200);
		} else {
			IJ.showMessage("wrong window type");
		}
	}
}
