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

public class wrap_angle_profile_jru_v1 implements PlugIn {
	//this plugin takes an intensity plot from -360 to 360
	//assume each plot doesn't encompass more than 360 total degrees
	//it takes values <=-180 and >180 and wraps them around

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][] newxvals=new float[npts.length][];
		float[][] newyvals=new float[npts.length][];
		for(int i=0;i<npts.length;i++){
			newxvals[i]=new float[xvals[i].length];
			newyvals[i]=new float[yvals[i].length];
			if(xvals[i][0]<=(-180.0f)){ //here we take the first segment and wrap it around to the end
				//find the position where we get to -180
				int pos=0;
				while(pos<npts[i] && xvals[i][pos]<=(-180.0f)){
					pos++;
				}
				//copy the ending portion
				for(int j=pos;j<npts[i];j++) newxvals[i][j-pos]=xvals[i][j];
				for(int j=pos;j<npts[i];j++) newyvals[i][j-pos]=yvals[i][j];
				//and then copy the beginning
				for(int j=0;j<pos;j++) newxvals[i][j+(npts[i]-pos)]=xvals[i][j]+360.0f;
				for(int j=0;j<pos;j++) newyvals[i][j+(npts[i]-pos)]=yvals[i][j];
			} else { //here we take the ending segment and wrap it around to the beginning
				//find the position where we get to 180
				int pos=0;
				while(pos<npts[i] && xvals[i][pos]<=180.0f){
					pos++;
				}
				//copy the ending portion
				for(int j=pos;j<npts[i];j++) newxvals[i][j-pos]=xvals[i][j]-360.0f;
				for(int j=pos;j<npts[i];j++) newyvals[i][j-pos]=yvals[i][j];
				//and then copy the beginning
				for(int j=0;j<pos;j++) newxvals[i][j+(npts[i]-pos)]=xvals[i][j];
				for(int j=0;j<pos;j++) newyvals[i][j+(npts[i]-pos)]=yvals[i][j];
			}
		}
		new PlotWindow4("Wrapped Angle Plot","x","y",newxvals,newyvals,npts).draw();
	}

}
