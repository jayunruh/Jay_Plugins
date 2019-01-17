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
import javax.script.*;

public class add_table_column_jru_v1 implements PlugIn {
	ScriptEngineManager manager;
	ScriptEngine engine;

	public void run(String arg) {
		manager=new ScriptEngineManager();
		engine=manager.getEngineByName("js");
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
			String[] new_col_labels=new String[col_labels.length];
			for(int i=0;i<col_labels.length;i++){
				new_col_labels[i]=(col_labels[i].trim()).replace(' ','_');
				if(new_col_labels[i].equals("")) new_col_labels[i]="_";
				new_col_labels[i]=new_col_labels[i].replace('*','_');
				new_col_labels[i]=new_col_labels[i].replace('.','_');
				new_col_labels[i]=new_col_labels[i].replace('-','_');
				new_col_labels[i]=new_col_labels[i].replace('/','_');
				new_col_labels[i]=new_col_labels[i].replace('#','_');
				new_col_labels[i]=new_col_labels[i].replace('^','_');
				//IJ.log(new_col_labels[i]);
			}
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Column_Labels (for reference)",new_col_labels,new_col_labels[0]);
			gd2.addStringField("Column_Name","",20);
			gd2.addStringField("Equation","",50);
			gd2.addCheckbox("Replace Original Table",true);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			gd2.getNextChoiceIndex();
			String colname=gd2.getNextString();
			String formula=gd2.getNextString();
			boolean replace=gd2.getNextBoolean();
			int nlines=tp.getLineCount();
			List<List<String>> listtable=table_tools.table2listtable(tp);
			for(int i=0;i<nlines;i++){
				List<String> arr=listtable.get(i);
				double temp=func(new_col_labels,arr,formula);
				arr.add(""+temp);
			}
			if(!replace) new TextWindow(tw.getTitle()+"_1",table_tools.print_string_array(col_labels)+"\t"+colname,table_tools.print_listtable(listtable),200,400);
			else table_tools.replace_table(tp,listtable,table_tools.print_string_array(col_labels)+"\t"+colname);
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public double func(String[] column_names,List<String> column_values,String function){
		StringBuffer script=new StringBuffer();
		for(int i=0;i<column_names.length;i++){
			try{
				double temp=Double.parseDouble(column_values.get(i));
				script.append(column_names[i]+"="+temp+"; ");
			}catch(NumberFormatException e){
			}
		}
		script.append("retval="+function+";");
		//IJ.log(script.toString());
		Double temp=new Double(0.0);
		try{
			temp=(Double)engine.eval(script.toString());
		}catch(Exception e){
			IJ.log(e.getMessage());
		}
		return temp.doubleValue();
	}
}
