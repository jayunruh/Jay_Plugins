/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class set_imaris_dataset_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		float[] dimensions={imp.getWidth(),imp.getHeight(),imp.getNChannels(),imp.getNSlices(),imp.getNFrames(),(float)jutils.get_psize(imp),(float)jutils.get_pdepth(imp)};
		boolean success=ImarisXT_utils.setDataSet(jutils.stack2array(imp.getStack()),dimensions);
	}

}
