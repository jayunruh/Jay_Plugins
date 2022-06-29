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
import ij.text.*;
import jguis.*;
import java.util.*;

public class table2image_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		List<List<String>> listtable=table_tools.table2listtable(tw[0].getTextPanel());
		int height=listtable.size();
		int width=listtable.get(0).size();
		float[] image=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				image[j+i*width]=table_tools.get_number(listtable,i,j);
			}
		}
		new ImagePlus("Table_Image",new FloatProcessor(width,height,image,null)).show();
	}

}
