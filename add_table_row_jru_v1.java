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
import java.awt.*;
import ij.plugin.*;
import ij.text.*;
import jguis.*;

public class add_table_row_jru_v1 implements PlugIn {

	public void run(String arg) {
		//first get the table window
		TextWindow[] tw=jutils.selectTables(false,1);
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.make_labels_unique(table_tools.getcollabels(tp));
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<col_labels.length;i++){
			gd.addStringField(col_labels[i],"0");
		}
		gd.showDialog(); if(gd.wasCanceled()) return;
		String[] newcols=new String[col_labels.length];
		for(int i=0;i<col_labels.length;i++){
			newcols[i]=gd.getNextString();
		}
		tp.append(table_tools.print_string_array(newcols));
	}

}
