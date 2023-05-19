/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.io.*;

public class set_default_directory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String currdir=OpenDialog.getDefaultDirectory();
		gd.addStringField("Directory",currdir,32);
		gd.addCheckbox("Use Dialog",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		String dir=gd.getNextString();
		if(gd.getNextBoolean()){
			DirectoryChooser dc=new DirectoryChooser("Choose Default Directory");
			dir=dc.getDirectory();
		}
		if(dir!=null) OpenDialog.setDefaultDirectory(dir);
	}

}
