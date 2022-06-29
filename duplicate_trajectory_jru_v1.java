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

public class duplicate_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		if(!jutils.is3DPlot(iw)){
			if(sel>=0) jutils.getPW4SelCopy(iw);
			else jutils.getPW4Copy(iw);
		} else {
			jutils.getPW3DSelCopy(iw,sel);
		}
	}

}
