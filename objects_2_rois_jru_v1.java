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

public class objects_2_rois_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		findblobs3 fb=new findblobs3(width,height);
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Show Squares",true);
		int halfsize=5;
		gd.addNumericField("square_halfsize",halfsize,0);
		int edgebuff=50;
		gd.addNumericField("edge_buffer",edgebuff,0);
		int mindist=10;
		gd.addNumericField("min_separation",mindist,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean squares=gd.getNextBoolean();
		halfsize=(int)gd.getNextNumber();
		edgebuff=(int)gd.getNextNumber();
		mindist=(int)gd.getNextNumber();
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			rman=new RoiManager();
		}
		rman.runCommand("show all");
		//for(int i=0;i<slices;i++){
			float[] objects=null;
			Object pix=imp.getProcessor().getPixels();
			if(pix instanceof float[]){
				objects=(float[])pix;
				fb.set_objects(objects);
			} else {
				objects=fb.dofindblobs(pix,0.5f);
			}
			if(mindist<=0.0f && edgebuff<=0.0f && !squares){
				for(int i=0;i<fb.nobjects;i++){
					Polygon poly=fb.get_object_outline(objects,i+1);
					rman.addRoi(new PolygonRoi(poly,Roi.FREEROI));
				}
			} else {
				float[][] centroids=measure_object.centroids(objects,width,height);
				boolean[] reject=new boolean[centroids.length];
				float mindist2=(float)mindist*(float)mindist;
				for(int i=0;i<(centroids.length-1);i++){
					if(!reject[i]){
						for(int j=i+1;j<centroids.length;j++){
							if(!reject[j]){
								float dist2=(centroids[j][0]-centroids[i][0])*(centroids[j][0]-centroids[i][0])+(centroids[j][1]-centroids[i][1])*(centroids[j][1]-centroids[i][1]);
								if(dist2<mindist2){
									reject[i]=true;
									reject[j]=true;
								}
							}
						}
					}
				}
				for(int i=0;i<centroids.length;i++){
					if(!reject[i]){
						if(squares){
							int x=(int)Math.round(centroids[i][0])-halfsize;
							int y=(int)Math.round(centroids[i][1])-halfsize;
							if(x>edgebuff && y>edgebuff && x<(width-edgebuff) && y<(height-edgebuff)){
								Roi roi=new Roi(x,y,halfsize*2,halfsize*2);
								//roi.setPosition(i+1);
								rman.addRoi(roi);
							}
						} else {
							Polygon poly=fb.get_object_outline(objects,i+1);
							rman.addRoi(new PolygonRoi(poly,Roi.FREEROI));
						}
					}
				}
			}
		//}
	}

}
