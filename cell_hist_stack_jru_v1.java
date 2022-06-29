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

public class cell_hist_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Cell Column",col_labels,col_labels[0]);
		gd.addChoice("Hist Column",col_labels,col_labels[0]);
		gd.addNumericField("Hist_Start",0.0,5,15,null);
		gd.addNumericField("Hist_End",100.0,5,15,null);
		gd.addNumericField("Bin_Size",1.0,5,15,null);
		gd.addCheckbox("Log_x",true);
		gd.addCheckbox("Keep_Separate",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int cellnameindex=gd.getNextChoiceIndex();
		int histcol=gd.getNextChoiceIndex();
		float histstart=(float)gd.getNextNumber();
		float histend=(float)gd.getNextNumber();
		float binsize=(float)gd.getNextNumber();
		boolean logx=gd.getNextBoolean();
		boolean separate=gd.getNextBoolean();
		List<List<String>> listtable=table_tools.table2listtable(tp);
		table_tools.sort_listtable(listtable,cellnameindex);
		List<String> celllist=table_tools.get_cell_list(listtable,cellnameindex);
		ImageStack histstack=null;
		for(int i=0;i<celllist.size();i++){
			String cellname=celllist.get(i);
			//get the cell table
			List<List<String>> subtable=table_tools.get_cell_listtable(listtable,cellnameindex,cellname);
			//now histogram the hist column
			float[] col=table_tools.get_column_array(subtable,histcol);
			PlotHist hist=new PlotHist(col_labels[histcol],"Occurences",col,3);
			hist.setLimits(new float[]{histstart,histend,0.0f,100.0f});
			hist.setBinSizeUnits(binsize);
			hist.setLogAxes(logx,false);
			hist.yautoscale();
			float[] limits=hist.getLimits();
			//limits[3]*=0.5f;
			limits[2]=0.0f;
			hist.setLimits(limits);
			if(separate) new PlotWindowHist(cellname,hist).draw();
			ColorProcessor cp=hist.getProcessor();
			if(histstack==null) histstack=new ImageStack(cp.getWidth(),cp.getHeight());
			histstack.addSlice(cellname,cp);
		}
		new ImagePlus("Cell Hist Stack",histstack).show();
	}

}
