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

public class dynamic_spectrum_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi roi=imp.getRoi();
		if(roi==null){
			imp.setRoi(imp.getWidth()/4,imp.getHeight()/4,imp.getWidth()/2,imp.getHeight()/2);
		}
		dynamic_profile_panel dp=new dynamic_profile_panel();
		dp.init(imp,"Avg");
		dynamic_profile_panel.launch_frame(dp);
	}

}
