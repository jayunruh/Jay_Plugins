/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;

public class get_2D_hist_roi_coords_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		//float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		//float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		//int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
		RoiManager rman=RoiManager.getInstance();
		Roi roi=rman.getRoisAsArray()[0];
		Polygon poly=roi.getPolygon();
		float[][] coords=new float[poly.npoints][];
		for(int i=0;i<poly.npoints;i++){
			Object[] args1={(Integer)poly.xpoints[i],(Integer)poly.ypoints[i]};
			coords[i]=(float[])jutils.runReflectionMethod(plot,"getPlotCoords",args1);
		}
		table_tools.create_table("Histogram Coordinates",coords,new String[]{"x","y"});
	}

}
