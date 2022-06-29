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

public class channel_alignment_jru_v1 implements PlugIn, DialogListener {
	float[][] shifts;
	Object[] cseries;
	Object[] cclone;
	int nchannels,nslices,nframes,datatype,width,height,currchan,currslice,currframe;
	ImagePlus imp;
	ImageStack stack;
	float psize,cpsize;

	public void run(String arg) {
		//this provides cartesian color alignment with a preview screen
		//note that alignments are in units of the current scaling (not pixels)
		imp=WindowManager.getCurrentImage();
		width=imp.getWidth(); height=imp.getHeight();
		nchannels=imp.getNChannels();
		nslices=imp.getNSlices();
		nframes=imp.getNFrames();
		stack=imp.getStack();
		psize=(float)jutils.get_psize(imp);
		cpsize=psize;
		currchan=imp.getChannel()-1;
		currslice=imp.getSlice()-1;
		currframe=imp.getFrame()-1;
		init_options();
		cseries=jutils.get3DCSeries(stack,currslice,currframe,nframes,nslices,nchannels);
		datatype=algutils.get_array_type(cseries[0]);
		cclone=algutils.clone_obj_array(cseries);
		if(showDialog()){
			for(int i=0;i<nframes;i++){
				for(int j=0;j<nslices;j++){
					if(i==currframe && j==currslice){
						for(int k=0;k<cclone.length;k++){
							float[] shifted=interpolation.shift_image(cclone[k],width,height,shifts[k][0]/cpsize,shifts[k][1]/cpsize);
							cclone[k]=algutils.convert_array(shifted,datatype);
						}
						jutils.set3DCSeries(cclone,stack,j,i,nframes,nslices,nchannels);
					} else {
						Object[] cseries2=jutils.get3DCSeries(stack,j,i,nframes,nslices,nchannels);
						for(int k=0;k<cseries2.length;k++){
							float[] shifted=interpolation.shift_image(cseries2[k],width,height,shifts[k][0]/cpsize,shifts[k][1]/cpsize);
							cseries2[k]=algutils.convert_array(shifted,datatype);
						}
						jutils.set3DCSeries(cseries2,stack,j,i,nframes,nslices,nchannels);
					}
				}
			}
			set_options();
		} else {
			jutils.set3DCSeries(cclone,stack,currslice,currchan,nframes,nslices,nchannels);
		}
		imp.setStack(stack);
		imp.updateAndDraw();
	}

	public boolean showDialog(){
		GenericDialog gd=new GenericDialog("Shift Parameters");
		gd.addCheckbox("Calibrate",true);
		for(int i=0;i<nchannels;i++){
			gd.addNumericField("ch"+(i+1)+"_x_shift",shifts[i][0],5,15,null);
			gd.addNumericField("ch"+(i+1)+"_y_shift",shifts[i][0],5,15,null);
		}
		gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		if(gd.getNextBoolean()) cpsize=psize;
		else cpsize=1.0f;
		for(int i=0;i<nchannels;i++){
			shifts[i][0]=(float)gd.getNextNumber();
			shifts[i][1]=(float)gd.getNextNumber();
		}

		return true;
	}

	public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		if(gd.getNextBoolean()) cpsize=psize;
		else cpsize=1.0f;
		for(int i=0;i<nchannels;i++){
			shifts[i][0]=(float)gd.getNextNumber();
			shifts[i][1]=(float)gd.getNextNumber();
			float[] shifted=interpolation.shift_image(cclone[i],width,height,shifts[i][0]/cpsize,shifts[i][1]/cpsize);
			cseries[i]=algutils.convert_array(shifted,datatype);
		}
		jutils.set3DCSeries(cseries,stack,currslice,currframe,nframes,nslices,nchannels);
		imp.setStack(stack);
		imp.updateAndDraw();
		return true;
	}

	public void init_options(){
		String[] options=jutils.get_plugin_options("channel_alignment_jru_v1",2*nchannels);
		if(options==null) set_options();
		shifts=new float[nchannels][2];
		for(int i=0;i<nchannels;i++){
			shifts[i][0]=Float.parseFloat(options[2*i]);
			shifts[i][1]=Float.parseFloat(options[2*i+1]);
		}
	}

	public void set_options(){
		if(shifts==null){
			shifts=new float[nchannels][2];
		}
		String[] options=new String[100];
		for(int i=0;i<options.length;i++) options[i]="0.0";
		for(int i=0;i<nchannels;i++){
			options[2*i]=""+shifts[i][0];
			options[2*i+1]=""+shifts[i][1];
		}
		jutils.set_plugin_options("channel_alignment_jru_v1",options);
	}

}
