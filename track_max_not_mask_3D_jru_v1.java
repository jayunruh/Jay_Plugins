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

public class track_max_not_mask_3D_jru_v1 implements PlugIn, DialogListener,  track_interface{
	int channels,slices,frames,currframe,width,height,linkdelay,minsize,currchan;
	float threshfraction,linkrange,searchrz,zedgebuf,zratio;
	boolean com;
	String fracstat;
	findblobs fb;
	ImageStack stack;
	ImagePlus imp;
	boolean output;
	float[][] outobjects;

	public void run(String arg) {
		imp=WindowManager.getCurrentImage();
		width=imp.getWidth(); height=imp.getHeight();
		stack=imp.getStack();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		channels=imp.getNChannels();
		currchan=imp.getC()-1;
		GenericDialog gd5=new GenericDialog("Options");
		gd5.addChoice("Threshhold Statistic",jstatistics.stats,jstatistics.stats[2]);
		gd5.addNumericField("Z_Edge_Buffer",4.0,5,15,null);
		gd5.addNumericField("Z_Ratio",3.0,5,15,null);
		gd5.showDialog(); if(gd5.wasCanceled()) return;
		fracstat=jstatistics.stats[gd5.getNextChoiceIndex()];
		searchrz=2.0f;
		zedgebuf=(float)gd5.getNextNumber();
		zratio=(float)gd5.getNextNumber();
		threshfraction=2.0f;
		if(fracstat.equals("Max")) threshfraction=0.5f;
		fb=new findblobs(width,height,new float[]{0.0f,10000.0f,0.0f,1000.0f,10,5,10});
		float[] criteria=get_criteria(imp,fb);
		if(criteria==null) return;
		imp.setHideOverlay(true);
		fb.setcriteria(criteria);
		fb.usemaxpt=(!com);
		currframe=1;
		tracker tk=new tracker(linkrange,linkdelay);
		//output=true;
		//outobjects=new float[slices][width*height];
		List<List<float[]>> trajlist=tk.track3D(this);
		if(outobjects!=null) new ImagePlus("Objects",jutils.array2stack(outobjects,width,height)).show();
		int maxlength=0;
		int ntraj=0;
		for(int i=0;i<trajlist.size();i++){ 
			if(trajlist.get(i).size()>maxlength){maxlength=trajlist.get(i).size();}
			if(trajlist.get(i).size()>=minsize){ntraj++;}
		}
		float[][] xvals=new float[ntraj][maxlength];
		float[][] yvals=new float[ntraj][maxlength];
		float[][] zvals=new float[ntraj][maxlength];
		int[] npts=new int[ntraj];
		String[] starts=new String[ntraj];
		String collabels="traj #\ttraj frame\tx\ty\tz\tintensity\tarea\tframe1";
		TextWindow tw=new TextWindow("Traj Data",collabels,"",400,200);
		int counter=0;
		for(int i=0;i<trajlist.size();i++){
			List<float[]> temp=trajlist.get(i);
			if(temp.size()>=minsize){
				npts[counter]=temp.size();
				starts[counter]=""+((float[])temp.get(0))[5];
				for(int j=0;j<npts[counter];j++){
					float[] temp2=temp.get(j);
					xvals[counter][j]=temp2[0];
					yvals[counter][j]=temp2[1];
					zvals[counter][j]=temp2[2]*zratio;
					tw.append(""+(counter+1)+"\t"+j+"\t"+temp2[0]+"\t"+temp2[1]+"\t"+temp2[2]+"\t"+temp2[3]+"\t"+temp2[4]+"\t"+temp2[5]+"\n");
				}
				counter++;
			}
		}
		if(counter>0){
			Traj3D t3D=new Traj3D("x","y","z",xvals,yvals,zvals,npts);
			PlotWindow3D pw=new PlotWindow3D("Trajectories",t3D);
			pw.p3.setAnnotations(starts);
			pw.draw();
		}
	}

	public float[] get_criteria(ImagePlus imp,findblobs fb){
		GenericDialog gd=new GenericDialog("Adjust Threshold");
		gd.addNumericField("Min_separation (pix)",4,0);
		gd.addNumericField("Thresh_fraction",threshfraction,5,15,null);
		gd.addNumericField("XY_Edge_buffer (pix)",10,0);
		gd.addNumericField("Max_blobs",1000,0);
		gd.addSlider("Display_Frame",1,getNFrames(),imp.getFrame());
		gd.addSlider("Display_Slice",1,slices,imp.getCurrentSlice());
		gd.addNumericField("Link_Range (pix)",10.0,5,15,null);
		gd.addNumericField("Max_Link_Delay (frames)",1,0);
		gd.addNumericField("Min_Traj_Length (frames)",5,0);
		gd.addCheckbox("Track_Center_Of_Mass",true);
		gd.addNumericField("Min_separation_z (slices)",4,0);
		gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		//criteria are 0minarea, 1maxarea, 2searchd(nu), 3maxblobs, 4thresh, 5minsep, 6edgebuf
		float[] criteria={0.0f,10000.0f,0.0f,1000.0f,0.5f,4.0f,10.0f};
		criteria[5]=(float)gd.getNextNumber();
		threshfraction=(float)gd.getNextNumber();
		criteria[6]=(float)gd.getNextNumber();
		criteria[3]=(float)gd.getNextNumber();
		gd.getNextNumber(); //omit the frame slider value
		gd.getNextNumber(); //omit the slice slider value
		linkrange=(float)gd.getNextNumber();
		linkdelay=(int)gd.getNextNumber();
		minsize=(int)gd.getNextNumber();
		com=gd.getNextBoolean();
		searchrz=0.5f*(float)gd.getNextNumber();
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
		int slice=(int)gd.getNextNumber();
		linkrange=(float)gd.getNextNumber();
		linkdelay=(int)gd.getNextNumber();
		minsize=(int)gd.getNextNumber();
		com=gd.getNextBoolean();
		searchrz=0.5f*(float)gd.getNextNumber();
		fb.setcriteria(criteria);
		fb.usemaxpt=(!com);
		int tempframe=imp.getFrame();
		int tempslice=imp.getSlice();
		if(frame!=tempframe || slice!=tempslice){
			tempframe=frame;
			tempslice=slice;
			imp.setPosition(currchan+1,tempslice,tempframe);
			imp.updateAndDraw();
		}
		List<float[]> currparams=getFrameParams(tempframe);
		if(currparams!=null){
			Overlay overlay=new Overlay();
			int blobwidth=(int)criteria[5];
			for(int i=0;i<currparams.size();i++){
				float[] blob=currparams.get(i);
				float zdist=(float)Math.abs(blob[2]-(slice-1));
				float rad=getSpheroidRadius(zdist,criteria[5]);
				Roi roi=new OvalRoi((int)Math.round(blob[0]-rad),(int)Math.round(blob[1]-rad),2*(int)rad,2*(int)rad);
				roi.setStrokeColor(Color.red);
				overlay.add(roi);
			}
			imp.setOverlay(overlay);
		}
		return true;
	}

	public float getSpheroidRadius(float zdist,float width){
		float searchr=0.5f*width;
		float zratio=searchrz/searchr;
		if(zdist/zratio>=searchr) return -1.0f;
		float rad=(float)Math.sqrt(searchr*searchr-zdist*zdist/(zratio*zratio));
		return rad;
	}

	public int getNFrames(){
		return frames;
	}

	public List<float[]> getFrameParams(int frame){
		Object[] stack2=jutils.get3DZSeries(stack,currchan,frame-1,frames,slices,channels);
		float[][] pixels=algutils.convert_arr_float2(stack2);
		float[] spectrum=jstatistics.getspectrum(fracstat,pixels,null);
		//assume that the statistic on the spectrum is the same as on the whole stack
		float max=jstatistics.getstatistic(fracstat,spectrum,null);
		//IJ.log("stat val = "+max);
		fb.thresh=threshfraction*max;
		if(fb.thresh==0.0f) return null;
		float[][] blobstats=null;
		if(!output) blobstats=fb.dofindblobs3D(pixels,null,searchrz,zedgebuf);
		else blobstats=fb.dofindblobs3D(pixels,outobjects,searchrz,zedgebuf);
		output=false;
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
