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
		gd.addCheckbox("X_Vals Column?",hasxvals);
		boolean haszvals=false;
		gd.addCheckbox("Z_Vals Column?",haszvals);
		boolean haserrs=false;
		gd.addCheckbox("Y_Errs Column?",haserrs);
		gd.addCheckbox("Sort by y val?",false);
		gd.addCheckbox("Separate Values?",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		hasxvals=gd.getNextBoolean();
		haszvals=gd.getNextBoolean();
		haserrs=gd.getNextBoolean();
		boolean sorty=gd.getNextBoolean();
		boolean separate=gd.getNextBoolean();
		if(niframes[index] instanceof TextWindow){
			TextWindow tw=(TextWindow)niframes[index];
			TextPanel tp=tw.getTextPanel();
			String headings=tp.getColumnHeadings();
			String[] col_labels=split_string_tab(headings);
			GenericDialog gd2=new GenericDialog("Choose Columns");
			if(hasxvals){
				gd2.addChoice("X_Column",col_labels,col_labels[0]);
			}
			gd2.addChoice("Y_Column",col_labels,col_labels[0]);
			if(haszvals){
				gd2.addChoice("Z_Column",col_labels,col_labels[0]);
			}
			if(haserrs){
				gd2.addChoice("Y_Errs_Column",col_labels,col_labels[0]);
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
			int nlines2=0;
			for(int i=0;i<nlines;i++){
				//String[] data=split_string_tab(tp.getLine(i));
				String[] data=table_tools.split(tp.getLine(i),"\t",false);
				if(data[yindex].length()<=0) continue;
				for(int j=0;j<data.length;j++) if(!table_tools.is_number(data[j])) data[j]="NaN";
				if(hasxvals) xvals[nlines2]=Float.parseFloat(data[xindex]);
				else xvals[nlines2]=(float)i+1;
				yvals[nlines2]=Float.parseFloat(data[yindex]);
				if(haszvals) zvals[nlines2]=Float.parseFloat(data[zindex]);
				if(haserrs) errs[nlines2]=Float.parseFloat(data[errsindex]);
				nlines2++;
			}
			if(nlines2<nlines){
				yvals=(float[])algutils.get_subarray(yvals,0,nlines2);
				xvals=(float[])algutils.get_subarray(xvals,0,nlines2);
				if(haszvals) zvals=(float[])algutils.get_subarray(zvals,0,nlines2);
				if(haserrs) errs=(float[])algutils.get_subarray(errs,0,nlines2);
				nlines=nlines2;
			}
			if(sorty){
				int[] order=jsort.javasort_order(yvals);
				if(haserrs){
					float[] temperrs=new float[errs.length];
					for(int i=0;i<errs.length;i++){
						temperrs[i]=errs[order[i]];
					}
					errs=temperrs;
				}
			}
			String xlab=col_labels[xindex]; String ylab=col_labels[yindex];
			if(!hasxvals){xlab=ylab; ylab="Value";}
			if(separate){
				//in this case, each point becomes its own data set
				float[][] newxvals=new float[xvals.length][1];
				float[][] newyvals=new float[yvals.length][1];
				float[][] newzvals=null; if(haszvals) newzvals=new float[zvals.length][1];
				float[][] newerrs=null; if(haserrs) newerrs=new float[errs.length][1];
				for(int i=0;i<xvals.length;i++){
					newxvals[i][0]=xvals[i];
					newyvals[i][0]=yvals[i];
					if(haszvals) newzvals[i][0]=zvals[i];
					if(haserrs) newerrs[i][0]=errs[i];
				}
				if(!haszvals){
					PlotWindow4 pw=new PlotWindow4(tw.getTitle()+" Plot",xlab,ylab,newxvals,newyvals,null);
					if(haserrs) pw.addErrors(newerrs);
					pw.draw();
				}else{
					new PlotWindow3D(tw.getTitle()+" Plot",new Traj3D(xlab,ylab,col_labels[zindex],newxvals,newyvals,newzvals,null)).draw();
				}
			} else {
				if(!haszvals){
					PlotWindow4 pw=new PlotWindow4(tw.getTitle()+" Plot",xlab,ylab,xvals,yvals);
					if(haserrs) pw.addSeriesErrors(0,errs);
					pw.draw();
				}else{
					new PlotWindow3D(tw.getTitle()+" Plot",new Traj3D(xlab,ylab,col_labels[zindex],xvals,yvals,zvals)).draw();
				}
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
