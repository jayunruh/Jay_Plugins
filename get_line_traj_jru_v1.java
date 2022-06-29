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
import jalgs.*;
import jguis.*;

public class get_line_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		GenericDialog gd2=new GenericDialog("Select Images");
		gd2.addChoice("Binned_Image",titles,titles[0]);
		gd2.addChoice("Original_Image",titles,titles[1]);
		//ltime is 3.85 times 5
		gd2.addNumericField("Time_per_line(ms)",19.25,5,15,null);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		int index1=gd2.getNextChoiceIndex();
		int index2=gd2.getNextChoiceIndex();
		float ltime=(float)gd2.getNextNumber();
		ImagePlus imp1 = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		Roi roi=imp1.getRoi();
		int binnedthickness=(int)roi.getStrokeWidth();
		int binnedwidth=imp1.getWidth();
		int binnedheight=imp1.getHeight();
		int width=imp2.getWidth();
		int height=imp2.getHeight();
		int binbyx=width/binnedwidth;
		int binbyy=height/binnedheight;
		int newthickness=binnedthickness*binbyx;
		if(roi instanceof Line){
			int x1=binbyx*((Line)roi).x1; int x2=binbyx*((Line)roi).x2;
			int y1=binbyy*((Line)roi).y1; int y2=binbyy*((Line)roi).y2;
			ImageStack stack=imp2.getStack();
			float[] carpet=(float[])stack.getProcessor(1).convertToFloat().getPixels();
			float[][] traj=getlinetraj(carpet,width,height,x1,y1,x2,y2,newthickness,ltime,true);
			PlotWindow4 pw=new PlotWindow4("Carpet Trajectory","time(s)","Intensity",traj[0],traj[1]);
			pw.draw();
			int slices=stack.getSize();
			for(int i=1;i<slices;i++){
				float[] carpet2=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
				float[][] traj2=getlinetraj(carpet2,width,height,x1,y1,x2,y2,newthickness,ltime,true);
				pw.addPoints(traj2[0],traj2[1],true);
			}
		} else {
			PolygonRoi polyroi=(PolygonRoi)roi;
			Rectangle r=polyroi.getBounds();
			int[] xvals1=polyroi.getXCoordinates();
			int[] yvals1=polyroi.getYCoordinates();
			int npts=polyroi.getNCoordinates();
			int[] xvals=new int[npts];
			int[] yvals=new int[npts];
			int nlines=npts-1;
			for(int i=0;i<npts;i++){
				xvals[i]=(xvals1[i]+r.x)*binbyx;
				yvals[i]=(yvals1[i]+r.y)*binbyy;
			}
			ImageStack stack=imp2.getStack();
			float[] carpet=(float[])stack.getProcessor(1).convertToFloat().getPixels();
			float[][] traj=getpolylinetraj(carpet,width,height,xvals,yvals,newthickness,ltime);
			PlotWindow4 pw=new PlotWindow4("Carpet Trajectory","time(s)","Intensity",traj[0],traj[1]);
			pw.draw();
			int slices=stack.getSize();
			for(int i=1;i<slices;i++){
				float[] carpet2=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
				float[][] traj2=getpolylinetraj(carpet2,width,height,xvals,yvals,newthickness,ltime);
				pw.addPoints(traj2[0],traj2[1],true);
			}
		}
	}

	private float[][] getpolylinetraj(float[] carpet,int width,int height,int[] xvals,int[] yvals,int thickness,float linetime){
		int nlines=xvals.length-1;
		float[][][] temptraj=new float[nlines][2][];
		int totlength=0;
		for(int i=0;i<(nlines-1);i++){
			float[][] temp=getlinetraj(carpet,width,height,xvals[i],yvals[i],xvals[i+1],yvals[i+1],thickness,linetime,false);
			temptraj[i][0]=temp[0];
			temptraj[i][1]=temp[1];
			totlength+=temp[0].length;
		}
		float[][] temp=getlinetraj(carpet,width,height,xvals[nlines-1],yvals[nlines-1],xvals[nlines],yvals[nlines],thickness,linetime,true);
		temptraj[nlines-1][0]=temp[0];
		temptraj[nlines-1][1]=temp[1];
		totlength+=temp[0].length;
		float[][] outtraj=new float[2][totlength];
		int counter=0;
		for(int i=0;i<nlines;i++){
			System.arraycopy(temptraj[i][0],0,outtraj[0],counter,temptraj[i][0].length);
			System.arraycopy(temptraj[i][1],0,outtraj[1],counter,temptraj[i][1].length);
			counter+=temptraj[i][0].length;
		}
		temptraj=null;
		return outtraj;
	}

	private float[][] getlinetraj(float[] carpet,int width,int height,int x1,int y1,int x2,int y2,int thickness,float linetime,boolean incend){
		int tempy1=y1;
		int tempy2=y2;
		int tempx1=x1;
		int tempx2=x2;
		if(tempy1>tempy2){
			tempy2=y1;
			tempy1=y2;
			tempx1=x2;
			tempx2=x1;
		}
		float xinc=(float)(tempx2-tempx1)/(float)(tempy2-tempy1);
		int start=tempy1;
		if(start<0){start=0;}
		int end=tempy2;
		if(end>=height){end=(height-1);}
		float[][] traj=null;
		if(incend){
			traj=new float[2][end-start+1];
		} else {
			traj=new float[2][end-start];
		}
		for(int i=start;i<end;i++){
			traj[0][i-start]=0.001f*linetime*(float)i;
			float left=xinc*(float)(i-start)+(float)(i*width)+(float)tempx1-0.5f*(float)thickness;
			for(int j=0;j<thickness;j++){
				float position=left+(float)j;
				traj[1][i-start]+=interpolation.interp1D(carpet,carpet.length,position);
			}
		}
		if(incend){
			traj[0][end-start]=0.001f*linetime*(float)end;
			float left=xinc*(float)(end-start)+(float)(end*width)+(float)tempx1-0.5f*(float)thickness;
			for(int j=0;j<thickness;j++){
				float position=left+(float)j;
				traj[1][end-start]+=interpolation.interp1D(carpet,carpet.length,position);
			}
		}
		return traj;
	}

}
