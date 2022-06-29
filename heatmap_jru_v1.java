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
import jguis.*;
import ij.text.*;

public class heatmap_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we create a heat map from a table with row, column, and intensity values
		//use a scaling of 4 pixels per mm
		int cellwidth=4;
		int cellheight=4;
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
			String headings=tp.getColumnHeadings();
			String[] col_labels=table_tools.split_string_tab(headings);
			GenericDialog gd2=new GenericDialog("Choose Columns");
			gd2.addChoice("X Column",col_labels,col_labels[0]);
			gd2.addChoice("Y Column",col_labels,col_labels[0]);
			gd2.addChoice("Intensity Column",col_labels,col_labels[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int xindex=gd2.getNextChoiceIndex();
			int yindex=gd2.getNextChoiceIndex();
			int zindex=gd2.getNextChoiceIndex();
			List<List<String>> list=table_tools.table2listtable(tp);
			int[][] coords=new int[2][list.size()];
			float[] vals=new float[list.size()];
			for(int i=0;i<list.size();i++){
				coords[0][i]=(int)table_tools.get_number(list,i,xindex);
				coords[1][i]=(int)table_tools.get_number(list,i,yindex);
				vals[i]=table_tools.get_number(list,i,zindex);
			}
			int minx=coords[0][0]; int maxx=coords[0][0];
			int miny=coords[1][0]; int maxy=coords[1][0];
			for(int i=0;i<list.size()-1;i++){
				if(coords[0][i]>maxx) maxx=coords[0][i];
				if(coords[0][i]<minx) minx=coords[0][i];
				if(coords[1][i]>maxy) maxy=coords[1][i];
				if(coords[1][i]<miny) miny=coords[1][i];
			}
			int xsize=(maxx-minx+1)*cellwidth;
			int ysize=(maxy-miny+1)*cellheight;
			float[] image=new float[xsize*ysize];
			for(int i=0;i<list.size();i++){
				int xoff=(coords[0][i]-minx)*cellwidth;
				int yoff=(coords[1][i]-miny)*cellheight;
				for(int j=xoff;j<(xoff+cellwidth);j++){
					if(j<xsize){
					for(int k=yoff;k<(yoff+cellheight);k++){
						if(k<ysize) image[j+k*xsize]=vals[i];
					}
					}
				}
			}
			new ImagePlus("Heat Map",new FloatProcessor(xsize,ysize,image,null)).show();
		}
	}

}
