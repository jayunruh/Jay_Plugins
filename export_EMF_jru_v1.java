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
import ij.plugin.frame.*;
import ij.io.*;
import jguis.*;

public class export_EMF_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		//Plot4 p4=jutils.getPW4PlotCopy(iw);
		SaveDialog sd=new SaveDialog("Save EMF File",arg,".emf");
		String dir=sd.getDirectory();
		String name=sd.getFileName();
		String path=dir+name;
		Object plot=jutils.runReflectionMethod(iw,"getPlot",null,null);
		//IJ.log("got plot");
		//if(name.endsWith(
		if(plot!=null){
			if(name.endsWith(".ps")) jutils.runReflectionMethod(plot,"saveAsPS",new Object[]{path});
			else if(name.endsWith(".eps")) jutils.runReflectionMethod(plot,"saveAsPS",new Object[]{path});
			else if(name.endsWith(".pdf")) jutils.runReflectionMethod(plot,"saveAsPDF",new Object[]{path});
			else jutils.runReflectionMethod(plot,"saveAsEMF",new Object[]{path});	
		}
		//p4.saveAsEMF(dir+name);
	}

}
