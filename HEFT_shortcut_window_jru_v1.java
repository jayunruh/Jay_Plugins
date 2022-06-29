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

public class HEFT_shortcut_window_jru_v1 implements PlugIn {

	public void run(String arg) {
		handleextrafiletypes_importer2 heft=new handleextrafiletypes_importer2();
		heft.init();
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Use Custom LOCI Importer",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		heft.usecustomloci=gd.getNextBoolean();
		handleextrafiletypes_importer2.launch_frame(heft);
	}

}
