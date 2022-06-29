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
import ij.plugin.*;
import ij.plugin.frame.*;
import jalgs.*;
import jalgs.jseg.*;
import java.util.*;
import java.awt.Polygon;
import jguis.*;
import ij.text.*;

public class track_objects_jru_v2 implements PlugIn,  track_interface{
	ImageStack stack;
	findblobs3 fb;
	int currframe;
	int nframes;
	ImageStack measurestack;
	int measurech;
	String stat;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Max_Step_Size(pixels)",5.0,5,10,null);
		gd.addNumericField("Min_Traj_Length(frames)",5,0);
		gd.addNumericField("Max_Frames_Off",0,0);
		gd.addCheckbox("Measurement_Image?",false);
		gd.addChoice("Measurement_Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float maxdist=(float)gd.getNextNumber();
		int minsize=(int)gd.getNextNumber();
		int linkdelay=(int)gd.getNextNumber();
		boolean measure=gd.getNextBoolean();
		stat=jstatistics.stats[gd.getNextChoiceIndex()];
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		stack=imp.getStack();
		if(measure){
			ImagePlus measureimp=jutils.selectImages(false,1)[0];
			measurestack=measureimp.getStack();
			measurech=measureimp.getNChannels();
		}
		fb=new findblobs3(width,height);
		nframes=stack.getSize();
		currframe=0;
		tracker tk=new tracker(maxdist,linkdelay);
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
		String collabels="traj #\ttraj frame\tx\ty\tangle\tmajor\tminor";
		for(int i=0;i<measurech;i++) collabels+="\tintensity"+(i+1);
		collabels+="\tframe1";
		TextWindow tw=new TextWindow("Traj Data",collabels,"",400,200);
		int counter=0;
		for(int i=0;i<trajlist.size();i++){
			List<float[]> temp=trajlist.get(i);
			if(temp.size()>=minsize){
				npts[counter]=temp.size();
				for(int j=0;j<npts[counter];j++){
					float[] temp2=temp.get(j);
					xvals[counter][j]=temp2[0];
					yvals[counter][j]=temp2[1];
					String temp3=""+(counter+1)+"\t"+j+"\t"+temp2[0]+"\t"+temp2[1]+"\t"+temp2[2]+"\t"+temp2[3]+"\t"+temp2[4]+"\t"+temp2[5];
					for(int k=0;k<measurech;k++) temp3+="\t"+temp2[6+k];
					temp3+="\n";
					tw.append(temp3);
				}
				counter++;
			}
		}
		new PlotWindow4("Trajectories","x","y",xvals,yvals,npts).draw();
	}

	public int getNFrames(){
		return nframes;
	}

	public List<float[]> getNextFrameParams(){
		byte[] pix=(byte[])stack.getPixels(currframe+1);
		float[] objects=fb.dofindblobs(pix);
		Object[] measurement=null;
		if(measurestack!=null){
			measurement=jutils.get3DCSeries(measurestack,0,currframe,nframes,1,measurech);
		}
		Polygon[] poly=fb.get_object_outlines(objects);
		List<float[]> params=new ArrayList<float[]>();
		for(int i=0;i<poly.length;i++){
			float[] temp=measure_object.get_ellipse_parameters(poly[i]);
			if(measurestack!=null){
				float[] temp2=fb.get_object_stats(objects,i+1,measurement,stat);
				temp=combine_arrays(temp,temp2);
			}
			params.add(temp);
			//params[i]=measure_object.centroid(poly[i]);
		}
		currframe++;
		IJ.showProgress(currframe,nframes);
		return params;
	}

	private float[] combine_arrays(float[] arr1,float[] arr2){
		float[] temp=new float[arr1.length+arr2.length];
		System.arraycopy(arr1,0,temp,0,arr1.length);
		System.arraycopy(arr2,0,temp,arr1.length,arr2.length);
		return temp;
	}

	public void put_assignments(int[] assignments){}

	public void show_progress(int progress,int nframes){
		IJ.showProgress(progress,nframes);
	}

}
