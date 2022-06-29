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
import ij.plugin.*;
import ij.plugin.frame.*;
import jguis.*;
import jalgs.*;
import ij.text.*;
import java.util.*;

public class search_plot_annotations_jru_v1 implements PlugIn {
	//this version searches for substrings, not exact matches

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		GenericDialog gd=new GenericDialog("Name List");
		gd.addCheckbox("Use_Table_Column",false);
		gd.addTextAreas("",null,10,20);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean usetable=gd.getNextBoolean();
		String input=gd.getNextText();
		if(usetable){
			TextWindow[] tw=jutils.selectTables(false,1);
			if(tw!=null && tw[0]!=null){
				TextPanel tp=tw[0].getTextPanel();
				String[] col_labels=table_tools.getcollabels(tp);
				GenericDialog gd2=new GenericDialog("options");
				gd2.addChoice("Search_Column",col_labels,col_labels[0]);
				gd2.showDialog(); if(gd2.wasCanceled()) return;
				int colindex=gd2.getNextChoiceIndex();
				List<List<String>> listtable=table_tools.table2listtable(tp);
				String[] colvals=table_tools.get_listtable_column(listtable,colindex);
				input=table_tools.print_string_array(colvals,3);
			}
		}
		String[] labels=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		if(labels==null){
			labels=new String[npts.length];
			for(int i=0;i<npts.length;i++){
				labels[i]=""+(i+1);
			}
		}
		String xlab=(String)jutils.runPW4VoidMethod(iw,"getxLabel");
		String ylab=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		int[] foundindices=new int[labels.length];
		int nfound=0;
		String[] list=(new delimit_string(' ')).getrows(input);
		for(int i=0;i<list.length;i++){
			list[i]=list[i].trim();
			//IJ.log(""+list[i]);
		}
		//for(int i=0;i<list.length;i++) list[i]=list[i].toLowerCase();
		for(int i=0;i<list.length;i++){
			boolean found=false;
			for(int j=0;j<labels.length;j++){
				if(labels[j].indexOf(list[i])>=0){
					foundindices[nfound]=j;
					nfound++;
					found=true;
					break;
				}
			}
			if(found) IJ.log(list[i]);
		}
		int counter=0;
		float[][] newxvals=new float[nfound][];
		float[][] newyvals=new float[nfound][];
		int[] newnpts=new int[nfound];
		String[] newannot=new String[nfound];
		for(int i=0;i<nfound;i++){
			int ti=foundindices[i];
			newxvals[i]=(float[])algutils.get_subarray(xvals[ti],0,npts[ti]);
			newyvals[i]=(float[])algutils.get_subarray(yvals[ti],0,npts[ti]);
			newnpts[i]=npts[ti];
			newannot[i]=labels[ti];
		}
		new PlotWindow4("Search_Results",xlab,ylab,newxvals,newyvals,newnpts).draw();
	}

}
