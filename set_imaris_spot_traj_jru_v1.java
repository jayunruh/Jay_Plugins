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
import jguis.*;
import jalgs.*;

public class set_imaris_spot_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow[] iw=jutils.selectPlotFamily(false,1,new String[]{"Trajectory Plot"});
		if(iw==null || iw.length<1) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Pixel Units",false);
		gd.addNumericField("Spot Radius",2.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean pixunits=gd.getNextBoolean();
		float rad=(float)gd.getNextNumber();
		int sel=(Integer)jutils.runPW4VoidMethod(iw[0],"getSelected");
		if(sel<0) sel=0;
		int[] npts=((int[][])jutils.runPW4VoidMethod(iw[0],"getNpts"))[0];
		float[] x,y,z;
		if(npts[sel]==1){
			//here we have a bunch of single points in our trajectory window
			float[][] xt=((float[][])jutils.runPW4VoidMethod(iw[0],"getXValues"));
			float[][] yt=((float[][])jutils.runPW4VoidMethod(iw[0],"getYValues"));
			float[][] zt=((float[][][])jutils.runPW4VoidMethod(iw[0],"getZValues"))[0];
			x=new float[npts.length];
			y=new float[npts.length];
			z=new float[npts.length];
			for(int i=0;i<npts.length;i++){
				x[i]=xt[i][0];
				y[i]=yt[i][0];
				z[i]=zt[i][0];
			}
		} else {
			x=((float[][])jutils.runPW4VoidMethod(iw[0],"getXValues"))[sel];
			y=((float[][])jutils.runPW4VoidMethod(iw[0],"getYValues"))[sel];
			if(x.length>npts[sel]){
				x=(float[])algutils.get_subarray(x,0,npts[sel]);
				y=(float[])algutils.get_subarray(y,0,npts[sel]);
			}
			z=new float[npts[sel]];
			if(jutils.is3DPlot(iw[0])){ //this should be 3D plot, right?
				z=((float[][][])jutils.runPW4VoidMethod(iw[0],"getZValues"))[0][sel];
				if(x.length>npts[sel]) z=(float[])algutils.get_subarray(z,0,npts[sel]);
			}
		}
		String[] annot=(String[])jutils.runPW4VoidMethod(iw[0],"getAnnotations");
		int[] timeindices=new int[npts[sel]];
		if(npts[sel]==1) timeindices=new int[npts.length];
		if(annot!=null && npts[sel]>1){
			int start=(int)Float.parseFloat(annot[sel]);
			for(int i=0;i<npts[sel];i++){
				timeindices[i]=start+i;
			}
		}
		boolean success=ImarisXT_utils.setSpotsTraj(pixunits,new float[][]{x,y,z},timeindices,rad);
		if(!success) IJ.log("Imaris Transfer Failed");
	}

}
