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
import jguis.*;
import jalgs.*;
import ij.text.*;
import java.util.*;

public class set_multi_plot_offsets_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("X Offset List");
		gd.addCheckbox("Use_Table_Column",false);
		gd.addTextAreas("",null,10,20);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean usetable=gd.getNextBoolean();
		String input=gd.getNextText();
		float[] offsets=new float[npts.length];
		if(usetable){
			TextWindow[] tw=jutils.selectTables(false,1);
			if(tw!=null && tw[0]!=null){
				TextPanel tp=tw[0].getTextPanel();
				String[] col_labels=table_tools.getcollabels(tp);
				GenericDialog gd2=new GenericDialog("options");
				gd2.addChoice("Offset_Column",col_labels,col_labels[0]);
				gd2.showDialog(); if(gd2.wasCanceled()) return;
				int colindex=gd2.getNextChoiceIndex();
				List<List<String>> listtable=table_tools.table2listtable(tp);
				String[] colvals=table_tools.get_listtable_column(listtable,colindex);
				//input=table_tools.print_string_array(colvals,3);
				for(int i=0;i<npts.length;i++) offsets[i]=Float.parseFloat(colvals[i]);
			}
		} else {
			String[] list=(new delimit_string(' ')).getrows(input);
			for(int i=0;i<npts.length;i++){
				offsets[i]=Float.parseFloat(list[i]);
			}
		}
		for(int i=0;i<npts.length;i++){
			for(int j=0;j<npts[i];j++){
				xvals[i][j]-=offsets[i];
			}
		}
		jutils.runPW4VoidMethod(iw,"updatePlot");
		jutils.runPW4VoidMethod(iw,"xautoscale");
	}

}
