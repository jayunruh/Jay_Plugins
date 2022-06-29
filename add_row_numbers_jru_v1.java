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
import ij.plugin.*;
import ij.util.*;
import ij.text.*;
import java.awt.Frame;
import java.util.*;
import jguis.*;

public class add_row_numbers_jru_v1 implements PlugIn {

	public void run(String arg) {
		//first get the table window
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length];
		for(int i=0;i<niframes.length;i++){
			titles[i]=niframes[i].getTitle();
		}
		GenericDialog gd=new GenericDialog("Windows");
		gd.addChoice("Windows",titles,titles[0]);
		gd.addCheckbox("Replace Original Table",true);
		gd.addCheckbox("Start at 1",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		boolean replace=gd.getNextBoolean();
		boolean start1=gd.getNextBoolean();
		int addval=0;
		if(start1) addval=1;
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String newheadings=tp.getColumnHeadings()+"\trow";
			List<List<String>> table=table_tools.table2listtable(tp);
			for(int i=0;i<table.size();i++){
				table.get(i).add(""+(i+addval));
			}
			if(replace) table_tools.replace_table(tp,table,newheadings);
			else new TextWindow(tw.getTitle()+"_1",newheadings,table_tools.print_listtable(table),200,400);
		} else {
			IJ.showMessage("wrong window type");
		}
	}

}
