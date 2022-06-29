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
import ij.plugin.frame.RoiManager;
import jalgs.*;
import jalgs.jseg.*;
import ij.text.*;

public class object_measure_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length+1];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		titles[wList.length]="Roi_Manager";
		GenericDialog gd2=new GenericDialog("Measurement Options");
		gd2.addChoice("Object Image",titles,titles[0]);
		gd2.addChoice("Measurement Image",titles,titles[0]);
		String[] stats=jstatistics.stats;
		gd2.addChoice("Statistic",stats,stats[0]);
		gd2.addCheckbox("Measure_Circ",false);
		gd2.addNumericField("Circ_Radius",4,0);
		gd2.addCheckbox("Add_To_Roi_Manager",false);
		gd2.addCheckbox("Add_Area",false);
		gd2.addCheckbox("Add_Slice_Label",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		int index1=gd2.getNextChoiceIndex();
		int index2=gd2.getNextChoiceIndex();
		ImagePlus imp1=null;
		if(index1!=wList.length) imp1 = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		String stat=gd2.getNextChoice();
		boolean circ=gd2.getNextBoolean();
		int circrad=(int)gd2.getNextNumber();
		boolean addroi=gd2.getNextBoolean();
		boolean addarea=gd2.getNextBoolean();
		boolean addlabel=gd2.getNextBoolean();

		ImageStack measstack=imp2.getStack();
		int width=imp2.getWidth(); int height=imp2.getHeight(); int slices=measstack.getSize();
		int nchans=imp2.getNChannels();
		ImageStack objstack=null;
		Polygon[] polys=null;
		if(imp1!=null){
			objstack=imp1.getStack();
		} else {
			Roi[] rois=RoiManager.getInstance().getRoisAsArray();
			polys=new Polygon[rois.length];
			for(int i=0;i<rois.length;i++) polys[i]=rois[i].getPolygon();
		}

		StringBuffer sb=new StringBuffer();
		int objslices=objstack.getSize();
		for(int i=0;i<slices;i+=nchans){
			findblobs3 fb=new findblobs3(width,height);
			float[] objects=null;
			if(imp1!=null){
				int tempslice=(int)((float)i/(float)nchans)+1;
				if(tempslice>=objslices) tempslice=objslices;
				Object tempobj=objstack.getPixels(tempslice);
				if(tempobj instanceof float[]){
					objects=(float[])tempobj;
					fb.set_objects(objects);
				} else {
					objects=fb.dofindblobs((byte[])tempobj);
				}
			} else {
				objects=fb.outlines2objects(polys);
			}
			
			int[][] objlims=fb.getallfilllimits(objects);
			Object[] measpixels=new Object[nchans];
			for(int j=0;j<nchans;j++) measpixels[j]=measstack.getPixels(i+j+1);
			int[] areas=fb.get_areas(objects);
			float[][] statvals=new float[nchans][fb.nobjects];
			float[][] circvals=null;
			float[] circobj=null;
			int[][] circlims=null;
			if(circ){
				circobj=fb.get_circ(objects,circrad);
				circlims=fb.getallfilllimits(circobj);
				circvals=new float[nchans][fb.nobjects];
			}

			for(int j=1;j<=fb.nobjects;j++){
				for(int k=0;k<nchans;k++){
					statvals[k][j-1]=fb.get_object_stats(objects,j,measpixels[k],objlims[j-1],stat);
					if(circ) circvals[k][j-1]=fb.get_object_stats(circobj,j,measpixels[k],circlims[j-1],stat);
				}
			}
			String label=measstack.getSliceLabel(i+1);
			for(int j=1;j<=fb.nobjects;j++){
				sb.append(""+(i/nchans+1)+"\t"+j);
				for(int k=0;k<nchans;k++) sb.append("\t"+statvals[k][j-1]);
				if(circ){for(int k=0;k<nchans;k++){sb.append("\t"+circvals[k][j-1]);}}
				if(addarea){sb.append("\t"+areas[j-1]);}
				if(addlabel){sb.append("\t"+label);}
				sb.append("\n");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		String headings="Slice\tObject";
		for(int i=0;i<nchans;i++) headings+="\tStatistic"+(i+1);
		if(circ){
			for(int i=0;i<nchans;i++) headings+="\tCirc"+(i+1);
		}
		if(addarea) headings+="\tarea";
		if(addlabel) headings+="\tlabel";
		new TextWindow("Object Measurements",headings,sb.toString(),400,200);
		if(imp1!=null && addroi){
			Object tempobj=objstack.getPixels(1);
			float[] objects=null;
			findblobs3 fb=new findblobs3(width,height);
			if(tempobj instanceof float[]){
				objects=(float[])tempobj;
				fb.set_objects(objects);
			} else {
				objects=fb.dofindblobs((byte[])tempobj);
			}
			Polygon[] polys2=fb.get_object_outlines(objects);
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			for(int i=0;i<polys2.length;i++){
				rman.addRoi(new PolygonRoi(polys2[i],Roi.FREEROI));
			}
		}
	}

}
