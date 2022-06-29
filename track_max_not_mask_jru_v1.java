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
import java.awt.Frame;
import java.awt.Color;
import java.awt.AWTEvent;
import ij.text.*;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import jalgs.jseg.*;
import jalgs.*;
import jguis.*;

public class track_max_not_mask_jru_v1 implements PlugIn, DialogListener,  track_interface{
	int slices,currframe,width,height,linkdelay,minsize;
	float threshfraction,linkrange;
	boolean com;
	String fracstat;
	findblobs fb;
	ImageStack stack;
	ImagePlus imp;

	public void run(String arg) {
		imp=WindowManager.getCurrentImage();
		width=imp.getWidth(); height=imp.getHeight();
		stack=imp.getStack();
		slices=stack.getSize();
		GenericDialog gd5=new GenericDialog("Options");
		gd5.addChoice("Threshhold Statistic",jstatistics.stats,jstatistics.stats[2]);
		gd5.showDialog(); if(gd5.wasCanceled()) return;
		fracstat=jstatistics.stats[gd5.getNextChoiceIndex()];
		threshfraction=2.0f;
		if(fracstat.equals("Max")) threshfraction=0.5f;
		if(fracstat.equals("Identity")) threshfraction=(float)imp.getProcessor().getStatistics().max;
		fb=new findblobs(width,height,new float[]{0.0f,10000.0f,0.0f,1000.0f,10,5,10});
		float[] criteria=get_criteria(imp,fb);
		if(criteria==null) return;
		imp.setHideOverlay(true);
		fb.setcriteria(criteria);
		fb.usemaxpt=(!com);
		currframe=1;
		tracker tk=new tracker(linkrange,linkdelay);
		List<List<float[]>> trajlist=tk.track2D(this);
		int maxlength=0;
		int ntraj=0;
		for(int i=0;i<trajlist.size();i++){ 
			if(trajlist.get(i).size()>maxlength){maxlength=trajlist.get(i).size();}
			if(trajlist.get(i).size()>=minsize){ntraj++;}
		}
		float[][] xvals=new float[ntraj][maxlength];
		float[][] yvals=new float[ntraj][maxlength];
		int[] npts=new int[ntraj];
		String[] starts=new String[ntraj];
		String collabels="traj #\ttraj frame\tx\ty\tintensity\tarea\tframe1";
		TextWindow tw=new TextWindow("Traj Data",collabels,"",400,200);
		int counter=0;
		for(int i=0;i<trajlist.size();i++){
			List<float[]> temp=trajlist.get(i);
			if(temp.size()>=minsize){
				npts[counter]=temp.size();
				starts[counter]=""+((float[])temp.get(0))[4];
				for(int j=0;j<npts[counter];j++){
					float[] temp2=temp.get(j);
					xvals[counter][j]=temp2[0];
					yvals[counter][j]=temp2[1];
					tw.append(""+(counter+1)+"\t"+j+"\t"+temp2[0]+"\t"+temp2[1]+"\t"+temp2[2]+"\t"+temp2[3]+"\t"+temp2[4]+"\n");
				}
				counter++;
			}
		}
		if(counter>0){
			PlotWindow4 pw=new PlotWindow4("Trajectories","x","y",xvals,yvals,npts);
			pw.p3.setAnnotations(starts);
			pw.draw();
		}
	}

	public float[] get_criteria(ImagePlus imp,findblobs fb){
		GenericDialog gd=new GenericDialog("Adjust Threshold");
		gd.addNumericField("Min_separation (pix)",4,0);
		gd.addNumericField("Thresh_fraction",threshfraction,5,15,null);
		gd.addNumericField("Edge_buffer (pix)",10,0);
		gd.addNumericField("Max_blobs",1000,0);
		gd.addSlider("Display_Frame",1,getNFrames(),imp.getCurrentSlice());
		gd.addNumericField("Link_Range (pix)",10.0,5,15,null);
		gd.addNumericField("Max_Link_Delay (frames)",1,0);
		gd.addNumericField("Min_Traj_Length (frames)",5,0);
		gd.addCheckbox("Track_Center_Of_Mass",true);
		gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		//criteria are 0minarea, 1maxarea, 2searchd(nu), 3maxblobs, 4thresh, 5minsep, 6edgebuf
		float[] criteria={0.0f,10000.0f,0.0f,1000.0f,0.5f,4.0f,10.0f};
		criteria[5]=(float)gd.getNextNumber();
		threshfraction=(float)gd.getNextNumber();
		criteria[6]=(float)gd.getNextNumber();
		criteria[3]=(float)gd.getNextNumber();
		gd.getNextNumber(); //omit the frame slider value
		linkrange=(float)gd.getNextNumber();
		linkdelay=(int)gd.getNextNumber();
		minsize=(int)gd.getNextNumber();
		com=gd.getNextBoolean();
		return criteria;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
		//here we update the current frame's overlay with circles to allow for parameter selection
		if(IJ.escapePressed()){return false;}
		float[] criteria={0.0f,10000.0f,0.0f,10000.0f,0.5f,4.0f,10.0f};
		criteria[5]=(float)gd.getNextNumber();
		threshfraction=(float)gd.getNextNumber();
		criteria[6]=(float)gd.getNextNumber();
		criteria[3]=(float)gd.getNextNumber();
		int frame=(int)gd.getNextNumber();
		linkrange=(float)gd.getNextNumber();
		linkdelay=(int)gd.getNextNumber();
		minsize=(int)gd.getNextNumber();
		com=gd.getNextBoolean();
		fb.setcriteria(criteria);
		fb.usemaxpt=(!com);
		int tempframe=imp.getCurrentSlice();
		if(frame!=tempframe){
			tempframe=frame;
			imp.setPosition(tempframe);
			imp.updateAndDraw();
		}
		List<float[]> currparams=getFrameParams(tempframe);
		if(currparams!=null){
			Overlay overlay=new Overlay();
			int blobwidth=(int)criteria[5];
			for(int i=0;i<currparams.size();i++){
				float[] blob=currparams.get(i);
				Roi roi=new OvalRoi((int)Math.round(blob[0]-0.5f*criteria[5]),(int)Math.round(blob[1]-0.5f*criteria[5]),blobwidth,blobwidth);
				roi.setStrokeColor(Color.red);
				overlay.add(roi);
			}
			imp.setOverlay(overlay);
		}
		return true;
	}

	public int getNFrames(){
		return slices;
	}

	public List<float[]> getFrameParams(int frame){
		float[] pixels=algutils.convert_arr_float2(stack.getPixels(frame));
		float max=jstatistics.getstatistic(fracstat,pixels,null);
		fb.thresh=threshfraction*max;
		if(fb.thresh==0.0f) return null;
		float[][] blobstats=fb.dofindblobs2(pixels,new float[width*height]);
		List<float[]> blobstats2=new ArrayList<float[]>();
		for(int i=0;i<blobstats.length;i++){
			blobstats2.add(blobstats[i]);
		}
		return blobstats2;
	}

	public List<float[]> getNextFrameParams(){
		List<float[]> blobstats=getFrameParams(currframe);
		currframe++;
		return blobstats;
	}

	public void show_progress(int progress,int nframes){
		IJ.showProgress(progress,nframes);
	}

	public void put_assignments(int[] assignments){}

}
