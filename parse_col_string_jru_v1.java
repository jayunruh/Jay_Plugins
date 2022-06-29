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

public class parse_col_string_jru_v1 implements PlugIn {
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
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Parse_Column",col_labels,col_labels[0]);
			gd2.addStringField("New_Column Name","",20);
			String[] meth={"Delimited","Fixed_Width"};
			gd2.addChoice("Parse_Method",meth,meth[0]);
			gd2.addStringField("Start_Delim_Str","-",10);
			gd2.addStringField("End_Delim_Str","-",10);
			gd2.addNumericField("Start_Pos (FW)",0,0);
			gd2.addNumericField("End_Pos (FW)",3,0);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int pcol=gd2.getNextChoiceIndex();
			String colname=gd2.getNextString();
			int parsemeth=gd2.getNextChoiceIndex();
			String sdelim=gd2.getNextString();
			String edelim=gd2.getNextString();
			int startpos=(int)gd2.getNextNumber();
			int endpos=(int)gd2.getNextNumber();
			int nlines=tp.getLineCount();
			List<List<String>> listtable=table_tools.table2listtable(tp);
			for(int i=0;i<nlines;i++){
				List<String> arr=listtable.get(i);
				String val=arr.get(pcol);
				String temp=null;
				if(parsemeth==0) temp=delim(val,sdelim,edelim);
				else temp=fw(val,startpos,endpos);
				arr.add(""+temp);
			}
			TextWindow tw2=new TextWindow(tw.getTitle()+"_1",table_tools.print_string_array(col_labels)+"\t"+colname,table_tools.print_listtable(listtable),200,400);
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public String delim(String orig,String sdelim,String edelim){
		int startindex=orig.indexOf(sdelim);
		if(startindex<0) return null;
		int endindex=orig.indexOf(edelim,startindex+1);
		if(endindex<0) return null;
		return orig.substring(startindex+sdelim.length(),endindex);
	}

	public String fw(String orig,int start,int end){
		if(end>=orig.length()) return orig.substring(start);
		else return orig.substring(start,end);
	}
}
