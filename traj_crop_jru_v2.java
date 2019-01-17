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

public class traj_crop_jru_v2 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4Copy(iw);
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		int[] npts=pw.getNpts();
		float[] limits=pw.getLimits();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Min",limits[0],5,10,null);
		gd.addNumericField("X_Max",limits[1],5,10,null);
		gd.addNumericField("Y_Min",limits[2],5,10,null);
		gd.addNumericField("Y_Max",limits[3],5,10,null);
		gd.addCheckbox("Delete_Partials",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float xmin=(float)gd.getNextNumber();
		float xmax=(float)gd.getNextNumber();
		float ymin=(float)gd.getNextNumber();
		float ymax=(float)gd.getNextNumber();
		boolean deletepart=gd.getNextBoolean();
		boolean[] fate=new boolean[npts.length];
		if(deletepart){
			for(int i=0;i<npts.length;i++){
				for(int j=0;j<npts[i];j++){
					if(xvals[i][j]<xmin || xvals[i][j]>xmax || yvals[i][j]<ymin || yvals[i][j]>ymax){
						fate[i]=true;
						break;
					}
				}
			}
		} else {
			for(int i=0;i<npts.length;i++){
				if(xvals[i][0]<xmin || xvals[i][0]>xmax || yvals[i][0]<ymin || yvals[i][0]>ymax){
					boolean inside=false;
					for(int j=1;j<npts[i];j++){
						boolean temp=(xvals[i][j]<xmin || xvals[i][j]>xmax || yvals[i][j]<ymin || yvals[i][j]>ymax);
						if(!temp){
							inside=true;
							break;
						}
					}
					if(!inside){
						fate[i]=true;
					}
				}
			}
		}
		int counter=0;
		pw.deleteMultiSeries(fate,false);
	}

}
