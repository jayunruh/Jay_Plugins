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

public class plot_columns_jru_v1 implements PlugIn {

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
		boolean haszvals=false;
		gd.addCheckbox("Z Vals Column?",haszvals);
		boolean haserrs=false;
		gd.addCheckbox("Y Errs Column?",haserrs);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		hasxvals=gd.getNextBoolean();
		haszvals=gd.getNextBoolean();
		haserrs=gd.getNextBoolean();
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
			if(haszvals){
				gd2.addChoice("Z Column",col_labels,col_labels[0]);
			}
			if(haserrs){
				gd2.addChoice("Y Errs Column",col_labels,col_labels[0]);
			}
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int xindex=0;
			if(hasxvals){
				xindex=gd2.getNextChoiceIndex();
			}
			int yindex=gd2.getNextChoiceIndex();
			int zindex=0;
			if(haszvals){
				zindex=gd2.getNextChoiceIndex();
			}
			int errsindex=0;
			if(haserrs) errsindex=gd2.getNextChoiceIndex();
			int nlines=tp.getLineCount();
			float[] xvals=new float[nlines];
			float[] yvals=new float[nlines];
			float[] zvals=null;
			if(haszvals) zvals=new float[nlines];
			float[] errs=null;
			if(haserrs) errs=new float[nlines];
			for(int i=0;i<nlines;i++){
				String[] data=split_string_tab(tp.getLine(i));
				if(hasxvals) xvals[i]=Float.parseFloat(data[xindex]);
				else xvals[i]=(float)i+1;
				yvals[i]=Float.parseFloat(data[yindex]);
				if(haszvals) zvals[i]=Float.parseFloat(data[zindex]);
				if(haserrs) errs[i]=Float.parseFloat(data[errsindex]);
			}
			if(!haszvals){
				PlotWindow4 pw=new PlotWindow4(tw.getTitle()+" Plot",col_labels[xindex],col_labels[yindex],xvals,yvals);
				if(haserrs) pw.addSeriesErrors(0,errs);
				pw.draw();
			}else{
				new PlotWindow3D(tw.getTitle()+" Plot",new Traj3D(col_labels[xindex],col_labels[yindex],col_labels[zindex],xvals,yvals,zvals)).draw();
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
