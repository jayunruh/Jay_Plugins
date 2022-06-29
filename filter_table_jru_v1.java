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

public class filter_table_jru_v1 implements PlugIn {

	public void run(String arg) {
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
				new TextWindow("filtered",tp.getColumnHeadings(),table_tools.print_listtable(filtered),400,200);
			}
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public float[][] get_filters(String[] col_labels){
		String[] columnlabels={"Use Filter?","Filter Parameter","Filter Type","Filter Value"};
		String[] filterlabels=col_labels;
		String[] filters={"<","=",">"};
		Object[][] tabledata={{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"},
		{new Boolean(false),filterlabels,filters,"0.0"}};
		int[][] options=new int[8][4];
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Filters",columnlabels,tabledata,options);
		if(retvals==null){return null;}
		int nfilters=0;
		boolean[] usefilter=new boolean[8];
		for(int i=0;i<8;i++){
			usefilter[i]=((Boolean)retvals[i][0]).booleanValue();
			if(usefilter[i]){nfilters++;}
		}
		if(nfilters==0){return null;}
		float[][] outvals=new float[nfilters][3];
		int counter=0;
		for(int i=0;i<8;i++){
			if(usefilter[i]){
				outvals[counter][0]=(float)table_tools.find_string_array_index(col_labels,(String)retvals[i][1]);
				if(outvals[counter][0]<0.0f){return null;}
				outvals[counter][1]=(float)table_tools.find_string_array_index(filters,(String)retvals[i][2]);
				if(outvals[counter][1]<0.0f){return null;}
				try{
					outvals[counter][2]=Float.parseFloat((String)retvals[i][3]);
				} catch(NumberFormatException e){
					return null;
				}
				counter++;
			}
		}
		return outvals;
	}

}
