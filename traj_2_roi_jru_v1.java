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
import ij.plugin.frame.RoiManager;

public class traj_2_roi_jru_v1 implements PlugIn {
	float zratio;

	public void run(String arg) {
		//here we create a segmented line roi from a plot
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Image","Trajectory"});
		if(imps==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Roi Options",new String[]{"Roi_Manager","Overlay","Trails","Squares","Circles"},"Roi_Manager");
		gd.addCheckbox("Trails Persist",false);
		gd.addNumericField("Square Size",10,0);
		gd.addNumericField("Z_Ratio (if 3D)",3.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int optindex=gd.getNextChoiceIndex();
		boolean persist=gd.getNextBoolean();
		int swidth=(int)gd.getNextNumber();
		zratio=(float)gd.getNextNumber();
		//ImageWindow iw=WindowManager.getCurrentWindow();
		WindowManager.setCurrentWindow(imps[0].getWindow());
		ImageWindow iw=imps[1].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		boolean is3D=false;
		int[] npts=null;
		float[][] zvals=null;
		if(jutils.is3DPlot(iw)){
			is3D=true;
			npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
			zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];
		} else {
			npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		}
		int[] colors=(int[])jutils.runPW4VoidMethod(iw,"getColors");
		String[] starts=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		Color[] jcolors=new Color[Plot4.java_colors.length];
		for(int i=0;i<jcolors.length;i++) jcolors[i]=Plot4.java_colors[i];
		jcolors[0]=Color.white;
		if(optindex==0){
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			for(int i=0;i<xvals.length;i++){
				Roi roi=null;
				if(is3D) roi=traj2roi(xvals[i],yvals[i],zvals[i],npts[i]); //this is the only place we allow for 3D
				else roi=traj2roi(xvals[i],yvals[i],npts[i]);
				rman.addRoi(roi);
			}
		} else if(optindex==1){
			Overlay overlay=new Overlay();
			for(int i=0;i<xvals.length;i++){
				Roi roi=traj2roi(xvals[i],yvals[i],npts[i]);
				roi.setColor(jcolors[colors[i]%8]);
				overlay.add(roi);
			}
			imps[0].setOverlay(overlay);
			imps[0].updateAndDraw();
		} else if(optindex==2) {
			ImageStack stack=imps[0].getStack();
			int width=imps[0].getWidth(); int height=imps[0].getHeight();
			ImageStack overstack=new ImageStack(width,height);
			int frames=imps[0].getNFrames();
			if(frames==1) frames=imps[0].getNSlices();
			for(int i=0;i<frames;i++) overstack.addSlice("",new int[width*height]);
			int channels=imps[0].getNChannels();
			//IJ.log(""+npts.length);
			for(int i=0;i<npts.length;i++){
				int start=0;
				if(starts!=null){
					start=(int)Float.parseFloat(starts[i]);
				}
				Color color=jcolors[colors[i]%8];
				//Color color=Color.white;
				int[][] coords=traj2int(xvals[i],yvals[i],npts[i]);
				for(int j=(start+1);j<(start+npts[i]);j++){
					ImageProcessor frame=overstack.getProcessor(j+1);
					frame.setColor(color);
					for(int k=0;k<(j-start);k++){
						frame.drawLine(coords[0][k],coords[1][k],coords[0][k+1],coords[1][k+1]);
					}
				}
				if(persist){
					for(int j=(start+npts[i]);j<frames;j++){
						ImageProcessor frame=overstack.getProcessor(j+1);
						frame.setColor(color);
						for(int k=0;k<(npts[i]-1);k++){
							frame.drawLine(coords[0][k],coords[1][k],coords[0][k+1],coords[1][k+1]);
						}
					}
				}
			}
			/*for(int i=0;i<frames;i++){
				ImageRoi roi=new ImageRoi(0,0,overstack.getProcessor(i+1));
				roi.setZeroTransparent(true);
				roi.setPosition(i*channels+1);
				imps[0].setOverlay(new Overlay(roi));
			}
			imps[0].updateAndDraw();*/
			new ImagePlus("Traj Trails",overstack).show();
		} else {
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			int halfwidth=(int)(0.5f*(float)swidth);
			for(int i=0;i<xvals.length;i++){
				if(optindex==3) rman.addRoi(new Roi((int)xvals[i][0]-halfwidth,(int)yvals[i][0]-halfwidth,2*halfwidth,2*halfwidth));
				else rman.addRoi(new OvalRoi((int)xvals[i][0]-halfwidth,(int)yvals[i][0]-halfwidth,2*halfwidth,2*halfwidth));
			}
		}
	}

	public Roi traj2roi(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		if(npts==1) return new PointRoi(xvals2[0],yvals2[0]);
		else return new PolygonRoi(xvals2,yvals2,npts,Roi.POLYLINE);
	}

	public Roi traj2roi(float[] xvals,float[] yvals,float[] zvals,int npts){
		Roi roi=null;
		if(npts==1){
			roi=new PointRoi((int)xvals[0],(int)yvals[0]);
		} else {
			roi=traj2roi(xvals,yvals,npts);
		}
		roi.setPosition((int)(1.5+zvals[0]/zratio));
		IJ.log(""+(int)(1.5+zvals[0]/zratio));
		return roi;
	}

	public int[][] traj2int(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		return new int[][]{xvals2,yvals2};
	}

}
