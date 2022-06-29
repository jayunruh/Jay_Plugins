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
import jalgs.*;
import jalgs.jseg.*;
import jguis.*;
import ij.plugin.frame.RoiManager;
import ij.text.*;

public class line_radial_angle_measure_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		//assume there is only one centroid roi
		//need to find it
		int xc=-1;
		int yc=-1;
		int croi=-1;
		for(int k=0;k<rois.length;k++){
			int type=rois[k].getType();
			if(type!=Roi.LINE && type!=Roi.POLYLINE){
				Rectangle bounds=rois[k].getBounds();
				if(type==Roi.POINT){
					xc=bounds.x;
					yc=bounds.y;
				} else {
					int count=0;
					for(int i=bounds.y;i<(bounds.y+bounds.height);i++){
						for(int j=bounds.x;j<(bounds.x+bounds.height);j++){
							if(rois[k].contains(j,i)){
								xc+=j; yc+=i; count++;
							}
						}
					}
					xc=(int)((float)xc/(float)count);
					yc=(int)((float)yc/(float)count);
				}
				croi=k;
				break;
			}
		}
		if(croi<0){
			IJ.log("need a centroid or area marker");
		}
		TextWindow tw=jutils.selectTable("Radial Angle Measures");
		if(tw==null) tw=new TextWindow("Radial Angle Measures","Image\tRoi\tAngle\tdist\tcenterangle","",400,200);
		for(int i=0;i<rois.length;i++){
			if(i!=croi){
				int x1,x2,y1,y2;
				if(rois[i] instanceof Line){
					x1=((Line)rois[i]).x1;
					x2=((Line)rois[i]).x2;
					y1=((Line)rois[i]).y1;
					y2=((Line)rois[i]).y2;
				} else {
					Polygon poly=((PolygonRoi)rois[i]).getPolygon();
					x1=poly.xpoints[0];
					x2=poly.xpoints[1];
					y1=poly.ypoints[0];
					y2=poly.ypoints[1];
				}
				float xc2=0.5f*((float)x1+(float)x2);
				float yc2=0.5f*((float)y1+(float)y2);
				float angle=measure_object.get_inner_angle_points(new float[]{(float)xc,xc2,(float)x2},new float[]{(float)yc,yc2,(float)y2});
				float angledeg=(float)Math.toDegrees(angle);
				if(angledeg>90.0f) angledeg=180.0f-angledeg;
				float dist=(float)Math.sqrt((xc2-xc)*(xc2-xc)+(yc2-yc)*(yc2-yc));
				float centerangle=(float)Math.atan2(yc2-yc,xc2-xc);
				tw.append(imp.getTitle()+"\t"+(i+1)+"\t"+angledeg+"\t"+dist+"\t"+centerangle);
			}
		}
	}

	

}
