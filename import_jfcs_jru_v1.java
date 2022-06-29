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
import jguis.*;
import ij.io.*;
import java.io.*;

public class import_jfcs_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Brightness Correlation?",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean brightcorr=gd.getNextBoolean();
		OpenDialog od=new OpenDialog("Open jfcs File",arg);
		String dir=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		if(AutoCorrFitWindow.is_this(dir+name)){
			AutoCorrFitWindow.launch_from_file(dir+name,brightcorr);
		} else {
			if(CrossCorrFitWindow.is_this(dir+name)){
				CrossCorrFitWindow.launch_from_file(dir+name,brightcorr);
			}
		}
	}

}
