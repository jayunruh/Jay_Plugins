/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
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

public class roi_2_traj3D_jru_v1 implements PlugIn {
	//this plugin takes a set of points, lines, or polylines from the roi manager and makes a 3D trajectory

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		float zratio=(float)jutils.get_zratio(imp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		zratio=(float)gd.getNextNumber();
		int channels=imp.getNChannels();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		int npts=0;
		for(int i=0;i<rois.length;i++){
			if(rois[i] instanceof Line){
				npts+=2;
			} else if(rois[i] instanceof PolygonRoi){
				npts+=((PolygonRoi)rois[i]).getNCoordinates();
			}
		}
		float[][] coords=new float[3][npts];
		int counter=0;
		for(int i=0;i<rois.length;i++){
			int zpos=rois[i].getZPosition()-1;
			if(rois[i] instanceof Line){
				Polygon pts=((Line)rois[i]).getPoints();
				for(int j=0;j<2;j++){
					coords[0][counter]=(float)pts.xpoints[j];
					coords[1][counter]=(float)pts.ypoints[j];
					coords[2][counter]=zratio*(float)zpos;
					counter++;
				}
			} else {
				Polygon pts=rois[i].getPolygon();
				for(int j=0;j<pts.npoints;j++){
					coords[0][counter]=(float)pts.xpoints[j];
					coords[1][counter]=(float)pts.ypoints[j];
					coords[2][counter]=zratio*(float)zpos;
					counter++;
				}
			}
		}
		Traj3D t3D=new Traj3D("x","y","z",new float[][]{coords[0]},new float[][]{coords[1]},new float[][]{coords[2]},null);
		new PlotWindow3D("3D Trajectory",t3D).draw();
	}

}
