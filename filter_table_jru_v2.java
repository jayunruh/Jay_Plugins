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

public class filter_table_jru_v2 implements PlugIn {

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
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		boolean replace=gd.getNextBoolean();
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String[] col_labels=table_tools.getcollabels(tp);
			int ncols=col_labels.length;
			float[][] filters=get_filters(col_labels);
			//IJ.log(""+filters[0][0]+" , "+filters[0][1]+" , "+filters[0][2]);
			List<List<String>> listtable=table_tools.table2listtable(tp);
			//IJ.log(""+listtable.size());
			List<List<String>> filtered=null;
			if(filters!=null){
				filtered=table_tools.get_filtered_listtable(listtable,filters);
				if(filtered==null){
					return;
				}
				//IJ.log(""+filtered.size());
				if(!replace) new TextWindow("filtered",tp.getColumnHeadings(),table_tools.print_listtable(filtered),400,200);
				else table_tools.replace_table(tp,filtered,col_labels);
			}
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public float[][] get_filters(String[] col_labels){
		GenericDialog gd=new GenericDialog("Select Filter");
		gd.addCheckbox("Use Filter?",false);
		String[] labels2=table_tools.make_labels_unique(col_labels);
		gd.addChoice("Filter Parameter",labels2,labels2[0]);
		String[] filters={"<","=",">"};
		gd.addChoice("Filter_Type",filters,filters[0]);
		gd.addNumericField("Filter_Value",0.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		if(!gd.getNextBoolean()) return null;
		int colindex=gd.getNextChoiceIndex();
		float[][] outvals=new float[1][3];
		outvals[0][0]=colindex;
		outvals[0][1]=gd.getNextChoiceIndex();
		outvals[0][2]=(float)gd.getNextNumber();
		return outvals;
	}

}
