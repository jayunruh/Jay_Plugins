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

public class channel_math_jru_v1 implements PlugIn {
	float[] operators;
	int operator;
	Object[] cseries;
	int nchannels,nslices,nframes,datatype,width,height,currchan,currslice,currframe;
	ImagePlus imp;
	ImageStack stack;
	float psize,cpsize;

	public void run(String arg) {
		//this performs independent mathematical manipulation on each channel
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
		operators=new float[nchannels];
		cseries=jutils.get3DCSeries(stack,currslice,currframe,nframes,nslices,nchannels);
		datatype=algutils.get_array_type(cseries[0]);
		if(!showDialog()) return;
		int counter=0;
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nslices;j++){
				Object[] cseries=jutils.get3DCSeries(stack,j,i,nframes,nslices,nchannels);
				for(int k=0;k<cseries.length;k++){
					float[] temp=algutils.convert_arr_float2(cseries[k]);
					operate(temp,operators[k],operator);
					cseries[k]=algutils.convert_array2(temp,datatype);
					IJ.showProgress(counter,nframes*nslices*nchannels);
					counter++;
				}
				jutils.set3DCSeries(cseries,stack,j,i,nframes,nslices,nchannels);
			}
		}
		imp.setStack(stack);
		imp.updateAndDraw();
	}

	public boolean showDialog(){
		GenericDialog gd=new GenericDialog("Math Parameters");
		String[] operations=new String[11];
		operations[0]="add";
		operations[1]="subtract";
		operations[2]="multiply";
		operations[3]="divide by";
		operations[4]="divide op by";
		operations[5]="set min";
		operations[6]="set max";
		operations[7]="log base";
		operations[8]="to the power";
		operations[9]="power of";
		operations[10]="abs";
		gd.addChoice("operation",operations,operations[0]);
		for(int i=0;i<nchannels;i++){
			gd.addNumericField("ch"+(i+1)+"_operator",operators[i],5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		operator=gd.getNextChoiceIndex();
		for(int i=0;i<nchannels;i++){
			operators[i]=(float)gd.getNextNumber();
		}
		return true;
	}

	public void operate(float[] opvals,float operator,int opindex){
		int length=opvals.length;
		if(opindex==0){for(int i=0;i<length;i++){opvals[i]+=operator;}}
		if(opindex==1){for(int i=0;i<length;i++){opvals[i]-=operator;}}
		if(opindex==2){for(int i=0;i<length;i++){opvals[i]*=operator;}}
		if(opindex==3){for(int i=0;i<length;i++){opvals[i]/=operator;}}
		if(opindex==4){for(int i=0;i<length;i++){opvals[i]=operator/opvals[i];}}
		if(opindex==5){for(int i=0;i<length;i++){if(opvals[i]<operator){opvals[i]=operator;}}}
		if(opindex==6){for(int i=0;i<length;i++){if(opvals[i]>operator){opvals[i]=operator;}}}
		if(opindex==7){for(int i=0;i<length;i++){opvals[i]=(float)Math.log(opvals[i])/(float)Math.log(operator);}}
		if(opindex==8){for(int i=0;i<length;i++){opvals[i]=(float)Math.pow(opvals[i],operator);}}
		if(opindex==9){for(int i=0;i<length;i++){opvals[i]=(float)Math.pow(operator,opvals[i]);}}
		if(opindex==10){for(int i=0;i<length;i++){opvals[i]=Math.abs(opvals[i]);}}
	}

}
