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

public class new_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addStringField("Title","");
		int ncols=9;
		gd.addNumericField("Number of Columns",ncols,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String title=gd.getNextString();
		ncols=(int)gd.getNextNumber();
		String[] collabels=new String[ncols];
		GenericDialog gd2=new GenericDialog("Column Labels");
		for(int i=0;i<ncols;i++){
			gd2.addStringField("Label_"+(i+1),""+(i+1));
		}
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		StringBuffer sb=new StringBuffer();
		sb.append(gd2.getNextString());
		for(int i=1;i<ncols;i++){
			sb.append("\t"+gd2.getNextString());
		}
		new TextWindow(title,sb.toString(),"",200,400);
	}
}
