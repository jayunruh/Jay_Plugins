/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;
import jalgs.jseg.*;
import ij.text.*;

public class object_measure_3D_jru_v1 implements PlugIn, gui_interface {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Objects","Measurement"});
		if(imps==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Measurement_Stat",jstatistics.stats,jstatistics.stats[0]);
		gd.addCheckbox("Ouput_Areas_Centroids",true);
		gd.addCheckbox("Is_Already_Indexed",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		boolean outcentarea=gd.getNextBoolean();
		boolean indexed=gd.getNextBoolean();
		ImageStack objstack=imps[0].getStack();
		Object[] objpix=jutils.stack2array(objstack);
		ImageStack measstack=imps[1].getStack();
		int measchans=imps[1].getNChannels();
		int measslices=imps[1].getNSlices();
		int measframes=imps[1].getNFrames();
		Object[][] measpix=new Object[measchans][];
		for(int i=0;i<measchans;i++) measpix[i]=jutils.get3DZSeries(measstack,i,0,measframes,measslices,measchans);
		findblobs3D fb=new findblobs3D(imps[0].getWidth(),imps[0].getHeight(),objstack.getSize(),this);
		float[][] objects=null;
		if(!indexed) objects=fb.dofindblobs(objpix);
		else {
			objects=new float[objpix.length][];
			for(int i=0;i<objpix.length;i++) objects[i]=(float[])objpix[i];
			fb.set_objects(objects);
		}
		float[][] centareas=null;
		float[][] stats=null;
		if(stat.equals("Avg") || stat.equals("Sum")){
			centareas=fb.getCentroidsAreasAvgs(objects,measpix);
			stats=new float[measchans][fb.nobjects];
			for(int i=0;i<stats[0].length;i++){
				for(int j=0;j<stats.length;j++){
					stats[j][i]=centareas[i][j+4];
					if(stat.equals("Sum")) stats[j][i]*=centareas[i][3];
				}
			}
		} else {
			if(outcentarea) centareas=fb.getCentroidsAreas(objects);
			int[][] lims=fb.getallfilllimits(objects);
			stats=fb.get_all_object_stats(objects,measpix,lims,stat);
		}
		String headings="object";
		if(outcentarea) headings+="\tx\ty\tz\tarea";
		for(int i=0;i<measchans;i++) headings+="\t"+stat+(i+1);
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<fb.nobjects;i++){
			sb.append(""+i);
			if(outcentarea) sb.append("\t"+centareas[i][0]+"\t"+centareas[i][1]+"\t"+centareas[i][2]+"\t"+centareas[i][3]);
			for(int j=0;j<measchans;j++) sb.append("\t"+stats[j][i]);
			sb.append("\n");
		}
		new TextWindow("Object Measurements",headings,sb.toString(),400,200);
	}

	public void showMessage(String message){}

	public void showProgress(int currpos,int finalpos){IJ.showProgress(currpos,finalpos);}

}
