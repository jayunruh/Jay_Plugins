/*******************************************************************************
 * Copyright (c) 2020 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.frame.*;
import ij.plugin.*;
import jguis.*;
import jalgs.*;
import ij.io.*;

public class manual_hist_gate_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		//float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		//float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		//int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float mag=(Float)jutils.runPW4VoidMethod(iw,"getmagnification");
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.error("must have roi in RoiManager");
			return;
		}
		Roi roi=rman.getRoisAsArray()[0];
		//we should support polyline or line rois
		int leftborder=	(int)(mag*Plot2DHist.LEFT_MARGIN);
		int rightborder=(int)(mag*(Plot2DHist.LEFT_MARGIN+Plot2DHist.WIDTH));
		int bottomborder=(int)(mag*(Plot2DHist.HEIGHT+Plot2DHist.TOP_MARGIN));
		int width=(int)(mag*Plot2DHist.WIDTH);
		if(!(roi instanceof Line || roi instanceof PolygonRoi)){
			IJ.error("roi must be line or polyline going from left to right");
			return;
		}
		Polygon poly=roi.getPolygon();
		int[] xpts=poly.xpoints;
		int[] ypts=poly.ypoints;
		if(xpts[0]>=xpts[xpts.length-1]){
			IJ.error("roi must be line or polyline going from left to right");
			return;
		}
		//now get the interpolated line profile
		int nlines=xpts.length-1;
		int[][] profile=null;
		if(nlines==1) profile=getYProfile(xpts[0],ypts[0],xpts[1],ypts[1],true);
		else profile=getYProfile(xpts[0],ypts[0],xpts[1],ypts[1],false);
		for(int i=1;i<nlines;i++){
			boolean last=(i==(nlines-1));
			int[][] temp=getYProfile(xpts[i],ypts[i],xpts[i+1],ypts[i+1],last);
			profile[0]=(int[])algutils.combine_arrays(profile[0],temp[0]);
			profile[1]=(int[])algutils.combine_arrays(profile[1],temp[1]);
		}
		//now crop the line profile at the image edges and make it a polygon gate
		int[][] cropped=cropProfile(profile,leftborder,leftborder+width-1);
		cropped[0]=(int[])algutils.combine_arrays(cropped[0],new int[]{rightborder-1,leftborder});
		cropped[1]=(int[])algutils.combine_arrays(cropped[1],new int[]{bottomborder-1,bottomborder-1});
		new PlotWindow4("ROI Profile","xpos","ypos",algutils.convert_arr_float(cropped[0]),algutils.convert_arr_float(cropped[1])).draw();
		translateProfile(cropped,-90,-20);
		PolygonRoi roiout=new PolygonRoi(cropped[0],cropped[1],cropped[0].length,Roi.POLYGON);
		rman.addRoi(roiout);	
		SaveDialog sd=new SaveDialog("Save Roi",arg,".roi");
		String outdir=sd.getDirectory();
		if(outdir==null) return;
		String roiname=sd.getFileName();
		multi_roi_writer.writeRoi(roiout,outdir+roiname);
	}

	public void translateProfile(int[][] profile,int xtrans,int ytrans){
		for(int i=0;i<profile[0].length;i++){
			profile[0][i]+=xtrans;
			profile[1][i]+=ytrans;
		}
		return;
	}

	public int[][] cropProfile(int[][] profile,int xstart,int xend){
		int[][] cropped=algutils.clone_multidim_array(profile);
		int startpos=0;
		int xpos=cropped[0][startpos];
		if(xpos>xstart){
			//starting point is too far right, prepend the left border
			cropped[0]=(int[])algutils.combine_arrays(new int[]{xstart},cropped[0]);
			cropped[1]=(int[])algutils.combine_arrays(new int[]{cropped[1][0]},cropped[1]);
		} else {
			while(xpos<xstart){
				startpos++;
				xpos=cropped[0][startpos];
			}
		}
		int endpos=startpos+1;
		xpos=cropped[0][endpos];
		while(xpos<xend && endpos<(cropped[0].length-1)){
			endpos++;
			xpos=cropped[0][endpos];
		}
		if(xpos<xend){
			//ending point is too far left, append the right border
			cropped[0]=(int[])algutils.combine_arrays(cropped[0],new int[]{xend});
			cropped[1]=(int[])algutils.combine_arrays(cropped[1],new int[]{cropped[1][cropped[1].length-1]});
		}
		//now do the cropping
		cropped[0]=(int[])algutils.get_subarray(cropped[0],startpos,(endpos-startpos+1));
		cropped[1]=(int[])algutils.get_subarray(cropped[1],startpos,(endpos-startpos+1));
		return cropped;
	}

	public int[][] getYProfile(int x1,int y1,int x2,int y2,boolean inclast){
		//here we get a y profile for a line
		float yinc=(float)(y2-y1)/(float)(x2-x1);
		int npts=x2-x1;
		if(inclast) npts++;
		int[] profile=new int[npts];
		int[] xprofile=new int[npts];
		for(int i=0;i<npts;i++){
			profile[i]=y1+(int)(yinc*(float)i);
			xprofile[i]=x1+i;
		}
		return new int[][]{xprofile,profile};
	}

}
