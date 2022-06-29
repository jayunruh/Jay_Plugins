/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
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

public class edit_column_labels_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		if(tw==null) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] cl=table_tools.getcollabels(tp);
		GenericDialog gd2=new GenericDialog("Edit Labels");
		for(int i=0;i<cl.length;i++){
			gd2.addStringField("Col"+(i+1),cl[i]);
		}
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		String[] cl2=new String[cl.length];
		for(int i=0;i<cl.length;i++) cl2[i]=gd2.getNextString();
		tp.updateColumnHeadings(table_tools.print_string_array(cl2));
		tp.updateDisplay();
	}

}
