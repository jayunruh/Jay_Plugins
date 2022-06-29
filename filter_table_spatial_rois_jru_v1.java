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
import jalgs.*;
import ij.plugin.frame.RoiManager;

public class filter_table_spatial_rois_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("X Column",col_labels,col_labels[0]);
		gd.addChoice("Y Column",col_labels,col_labels[1]);
		//gd.addNumericField("Max_Width",512,0);
		//gd.addNumericField("Max_Height",512,0);
		gd.addCheckbox("Add_Roi_Label",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int xcol=gd.getNextChoiceIndex();
		int ycol=gd.getNextChoiceIndex();
		//int width=(int)gd.getNextNumber();
		//int height=(int)gd.getNextNumber();
		float[] xvals=table_tools.get_column_array(listtable,xcol);
		float[] yvals=table_tools.get_column_array(listtable,ycol);
		int width=1+(int)jstatistics.getstatistic("Max",xvals,null);
		int height=1+(int)jstatistics.getstatistic("Max",yvals,null);
		boolean addlabel=gd.getNextBoolean();
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		//boolean[] mask=getMask(rois,width,height);
		float[] mask=getMask2(rois,width,height);
		List<List<String>> selected=new ArrayList<List<String>>();
		for(int i=0;i<listtable.size();i++){
			//float xval=table_tools.get_number(listtable,i,xcol);
			//float yval=table_tools.get_number(listtable,i,ycol);
			float xval=xvals[i];
			float yval=yvals[i];
			float sel=mask[(int)xval+width*(int)yval];
			if(sel>0.0f){
				List<String> row=listtable.get(i);
				if(addlabel) row.add(""+(int)sel);
				selected.add(row);
			}
		}
		if(addlabel){
			String[] col_labels2=new String[col_labels.length+1];
			for(int i=0;i<col_labels.length;i++) col_labels2[i]=col_labels[i];
			col_labels2[col_labels2.length-1]="Roi";
			table_tools.create_table("filtered",selected,col_labels2);
		} else {
			table_tools.create_table("filtered",selected,col_labels);
		}
	}

	public boolean[] getMask(Roi[] rois,int width,int height){
		boolean[] mask=new boolean[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(!mask[j+i*width]){
					for(int k=0;k<rois.length;k++){
						if(rois[k].contains(j,i)){
							mask[j+i*width]=true;
							break;
						}
					}
				}
			}
		}
		return mask;
	}

	public float[] getMask2(Roi[] rois,int width,int height){
		float[] mask=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(mask[j+i*width]==0.0f){
					for(int k=0;k<rois.length;k++){
						if(rois[k].contains(j,i)){
							mask[j+i*width]=(float)(k+1);
							break;
						}
					}
				}
			}
		}
		return mask;
	}

}
