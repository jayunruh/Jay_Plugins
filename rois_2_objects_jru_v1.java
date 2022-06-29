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
import jalgs.jseg.*;

public class rois_2_objects_jru_v1 implements PlugIn {

	public void run(String arg) {
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		Polygon[] polys=new Polygon[rois.length];
		int maxx=0;
		int maxy=0;
		for(int i=0;i<rois.length;i++){
			Rectangle r=rois[i].getBounds();
			if((r.x+r.width)>maxx) maxx=r.x+r.width;
			if((r.y+r.height)>maxy) maxy=r.y+r.height;
			polys[i]=rois[i].getPolygon();
		}
		//maxx++;
		//maxy++;
		//if(maxx<maxy) maxx=maxy;
		//else maxy=maxx;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Image_Width",maxx,0);
		gd.addNumericField("Image_Height",maxy,0);
		String[] schemes={"Indexed","Value","Binary"};
		gd.addChoice("Coloring Scheme",schemes,schemes[0]);
		gd.addNumericField("Coloring_Value",1.0f,5,15,null);
		gd.addNumericField("Point Roi Size",3,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		maxx=(int)gd.getNextNumber();
		maxy=(int)gd.getNextNumber();
		int scindex=gd.getNextChoiceIndex();
		float val=(float)gd.getNextNumber();
		int ptsize=(int)gd.getNextNumber();
		for(int i=0;i<rois.length;i++){
			Rectangle r=rois[i].getBounds();
			if(r.width<2 && r.height<2){
				polys[i]=getRectPoly(r.x,r.y,ptsize);
			}
		}
		findblobs3 fb=new findblobs3(maxx,maxy);
		float[] objects=fb.outlines2objects(polys);
		if(scindex==2) val=1.0f;
		if(scindex>0){
			for(int i=0;i<objects.length;i++){
				if(objects[i]>0.0f) objects[i]=val;
			}
		}
		new ImagePlus("Roi Objects",new FloatProcessor(maxx,maxy,objects,null)).show();
	}

	public Polygon getRectPoly(int x,int y,int size){
		int x0=x-size/2;
		int y0=y-size/2;
		int[] xvals={x0,x0+size,x0+size,x0};
		int[] yvals={y0,y0,y0+size,y0+size};
		return new Polygon(xvals,yvals,4);
	}

}
