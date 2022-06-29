/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;

public class multi_point_roi_2_roimanager_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		Roi roi=imp.getRoi();
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		if(roi instanceof PolygonRoi){
			Polygon poly=roi.getPolygon();
			for(int i=0;i<poly.npoints;i++){
				rman.addRoi(new PointRoi(poly.xpoints[i],poly.ypoints[i]));
			}
		}
	}

}
