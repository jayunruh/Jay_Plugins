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
import ij.plugin.frame.RoiManager;
import jguis.*;
import ij.io.*;

public class export_emf_image_rois_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		SaveDialog sd=new SaveDialog("Save EMF File",arg,".emf");
		String dir=sd.getDirectory();
		String name=sd.getFileName();
		if(name==null) return;
		String path=dir+name;
		EMFPlotter ep=new EMFPlotter(width,height,path);
		Image img=imp.getProcessor().createImage();
		ep.drawImage(img,0,0);
		for(int i=0;i<rois.length;i++){
			float linewidth=rois[i].getStrokeWidth();
			Color color=rois[i].getStrokeColor(); //getStrokeColor()?
			ep.setColor(color);
			//Polygon poly=((PolygonRoi)rois[i]).getPolygon();
			int[] xpts=null; int[] ypts=null;
			Polygon poly=rois[i].getPolygon();
			if(poly!=null){xpts=poly.xpoints; ypts=poly.ypoints;}
			else {
				Line temp=(Line)rois[i];
				xpts=new int[]{temp.x1,temp.x2};
				ypts=new int[]{temp.y1,temp.y2};
			}
			if(xpts!=null) ep.drawPolyLine(xpts,ypts,true);
		}
		ep.endPlotting();
		IJ.showStatus("finished exporting file");
	}

}
