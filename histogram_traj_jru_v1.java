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

public class histogram_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		String ylabel=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int selected=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		if(selected<0) selected=0;
		new PlotWindowHist("Histogram",ylabel,"Occurrences",yvals[selected],3).draw();
	}

}
