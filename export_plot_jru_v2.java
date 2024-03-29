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
import jguis.*;
import ij.io.*;
import java.io.*;

public class export_plot_jru_v2 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		String title=iw.getTitle();
		if(!title.endsWith(".txt")) title+=".txt";
		SaveDialog sd=new SaveDialog("Save Plot Object File",title,".txt");
		String dir=sd.getDirectory();
		String name=sd.getFileName();
		jutils.savePW4(iw,dir+File.separator+name,1);
	}

}
