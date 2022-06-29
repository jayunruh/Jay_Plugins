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

public class col_to_uppercase_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Select Column",col_labels,col_labels[0]);
		gd.addCheckbox("Upper Case",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int col=gd.getNextChoiceIndex();
		boolean upper=gd.getNextBoolean();
		for(int i=0;i<listtable.size();i++){
			String val=listtable.get(i).get(col);
			if(upper) val=val.toUpperCase();
			else val=val.toLowerCase();
			listtable.get(i).set(col,val);
		}
		table_tools.replace_table(tp,listtable,table_tools.print_string_array(col_labels));
	}

}
