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

public class DV_file_reader_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open DV Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		ImagePlus imp=(new DV_file_reader()).open(directory,fname);
		imp.show();
		if(imp.isHyperStack()){imp.setPosition(1,1,1);}
	}

}
