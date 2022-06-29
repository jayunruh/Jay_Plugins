/*******************************************************************************
 * Copyright (c) 2020 Jay Unruh, Stowers Institute for Medical Research.
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

public class transform_from_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin gets a table with the first two rows of an affine transform and applies it to an image stack
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Transformation Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		//the last 6 columns of the table will be used for the transformation (first two rows of affine matrix)
		int offset=col_labels.length-6;
		float[][][] transformations=new float[listtable.size()][2][3];
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Scale_Factor",1.0,5,15,null);
		gd.addCheckbox("Interpolate?",false);
		gd.addCheckbox("Global Transformations?",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float scalefactor=(float)gd.getNextNumber();
		boolean interp=gd.getNextBoolean();
		boolean global=gd.getNextBoolean();
		for(int i=0;i<listtable.size();i++){
			//List<String> row=listtable.get(i);
			for(int j=0;j<3;j++){
				transformations[i][0][j]=table_tools.get_number(listtable,i,j+offset);
				transformations[i][1][j]=table_tools.get_number(listtable,i,j+offset+3);
			}
			transformations[i][0][2]*=scalefactor;
			transformations[i][1][2]*=scalefactor;
		}
		ImagePlus[] imps=jutils.selectImages(false,1);
		if(imps==null) return;
		ImagePlus aligned=SIFTj.doTransformations(imps[0],transformations,interp,global);
		aligned.show();
	}

}
