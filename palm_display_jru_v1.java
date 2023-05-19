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
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import jalgs.*;
import jguis.*;

public class palm_display_jru_v1 implements PlugIn{

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();

		final PALMWindow cw = new PALMWindow();
		cw.mw=imp.getWidth();
		cw.init((float[])imp.getProcessor().getPixels());
		PALMWindow.launch_frame(cw);
	}
}

