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

public class roi_2_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we plot an roi
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi roi=imp.getRoi();
		Polygon roipoly=((PolygonRoi)roi).getPolygon();
		int npts=roipoly.npoints;
		float[] x=new float[npts];
		float[] y=new float[npts];
		for(int i=0;i<npts;i++){
			x[i]=(float)roipoly.xpoints[i];
			y[i]=(float)roipoly.ypoints[i];
		}
		new PlotWindow4("Roi Plot","x","y",x,y).draw();
	}

}
