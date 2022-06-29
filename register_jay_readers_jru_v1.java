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
import loci.formats.*;
import jguis.*;

public class register_jay_readers_jru_v1 implements PlugIn {

	public void run(String arg) {
		try{
			Class jayreader=Class.forName("jguis.loci_pw_reader_jru_v1");
			ImageReader.getDefaultReaderClasses().addClass(jayreader);
			//Class jayreader2=Class.forName("loci_sky_reader_jru_v1");
			//ImageReader.getDefaultReaderClasses().addClass(jayreader2);
			Class jayreader3=Class.forName("jguis.loci_kiss_reader_jru_v1");
			ImageReader.getDefaultReaderClasses().addClass(jayreader3);
		} catch(ClassNotFoundException e){
			//if the reader isn't intalled, don't bother
			IJ.log("loci pw or kiss reader wasn't found");
			return;
		}
	}

}
