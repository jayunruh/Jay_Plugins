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
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import ij.text.*;

public class table_listener_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		if(tw==null) return;
		ImagePlus[] imp=jutils.selectImages(false,1);
		if(imp==null) return;
		float zratio=(float)jutils.get_zratio(imp[0]);
		float psize=(float)jutils.get_psize(imp[0]);
		TextPanel tp=tw[0].getTextPanel();
		String[] collabels=table_tools.getcollabels(tp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("X_Column",collabels,collabels[0]);
		gd.addChoice("Y_Column",collabels,collabels[1]);
		gd.addChoice("Z_Column",collabels,collabels[2]);
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addCheckbox("Calibrated?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int xcol=gd.getNextChoiceIndex();
		int ycol=gd.getNextChoiceIndex();
		int zcol=gd.getNextChoiceIndex();
		zratio=(float)gd.getNextNumber();
		boolean calib=gd.getNextBoolean();
		IJTableListener ijtl=new IJTableListener();
		ijtl.init(tw[0],imp[0],xcol,ycol,zcol,zratio);
		if(calib) ijtl.psize=psize;
		IJTableListener.launch_frame(ijtl);
	}

}
