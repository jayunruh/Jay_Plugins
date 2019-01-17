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
import ij.util.*;
import ij.text.*;
import jguis.*;
import jalgs.*;

public class hist_columns_jru_v1 implements PlugIn {

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
		boolean hasxvals=false;
		gd.addCheckbox("X Vals Column?",hasxvals);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		hasxvals=gd.getNextBoolean();
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String headings=tp.getColumnHeadings();
			String[] col_labels=split_string_tab(headings);
			GenericDialog gd2=new GenericDialog("Choose Columns");
			if(hasxvals){
				gd2.addChoice("X Column",col_labels,col_labels[0]);
			}
			gd2.addChoice("Y Column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int xindex=0;
			if(hasxvals){
				xindex=gd2.getNextChoiceIndex();
			}
			int yindex=gd2.getNextChoiceIndex();
			int nlines=tp.getLineCount();
			float[] xvals=new float[nlines];
			float[] yvals=new float[nlines];
			if(hasxvals){
				int linect=0;
				for(int i=0;i<nlines;i++){
					String temp=tp.getLine(i);
					if(temp.length()<=0) continue;
					String[] data=table_tools.split(temp,"\t",false);
					if(data[yindex].length()<=0) continue;
					try{
						if(xindex<data.length) xvals[linect]=Float.parseFloat(data[xindex]);
						else xvals[linect]=Float.NaN;
						if(yindex<data.length) yvals[linect]=Float.parseFloat(data[yindex]);
						else yvals[linect]=Float.NaN;
						linect++;
					}catch(NumberFormatException e){}
				}
				if(linect!=nlines){
					xvals=(float[])algutils.get_subarray(xvals,0,linect);
					yvals=(float[])algutils.get_subarray(yvals,0,linect);
				}
				new PlotWindow2DHist(titles[index]+" Histogram",col_labels[xindex],col_labels[yindex],xvals,yvals,null).draw();
			} else {
				int linect=0;
				for(int i=0;i<nlines;i++){
					String temp=tp.getLine(i);
					if(temp.length()<=0) continue;
					String[] data=table_tools.split(temp,"\t",false);
					if(data[yindex].length()<=0) continue;
					try{
						if(yindex<data.length) yvals[linect]=Float.parseFloat(data[yindex]);
						else yvals[linect]=Float.NaN;
						linect++;
					}catch(NumberFormatException e){}
				}
				if(linect!=nlines){
					yvals=(float[])algutils.get_subarray(yvals,0,linect);
				}
				new PlotWindowHist(titles[index]+" Histogram",col_labels[yindex],"Occurrences",yvals,3).draw();
			}
		} else {
			IJ.showMessage("wrong window type");
		}
	}

	public String[] split_string_tab(String line){
		String temp;
		if(line.endsWith("\t")){
			temp=line.substring(0,line.length()-1);
		} else {
			temp=line.substring(0);
		}
		return Tools.split(temp,"\t");
	}
}
