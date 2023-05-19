/*******************************************************************************
 * Copyright (c) 2021 Jay Unruh, Stowers Institute for Medical Research.
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

public class remove_laser_flash_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Max_Flash_Duration",10,0);
		gd.addNumericField("End_Derivative_Threshold (fraction)",0.05f,5,15,null);
		gd.addCheckbox("Output_Derivative",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int maxduration=(int)gd.getNextNumber();
		float enderthresh=(float)gd.getNextNumber();
		boolean outder=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int currchan=imp.getC()-1;
		int currz=imp.getZ()-1;
		int nchan=imp.getNChannels();
		int nframes=imp.getNFrames();
		int nslices=imp.getNSlices();
		if(nframes==1){
			nframes=nslices;
			nslices=1;
		}
		Object[] tstack=jutils.get3DTSeries(stack,currz,currchan,nframes,nslices,nchan);
		//might want to try maximum here
		float[] profile=jstatistics.getspectrum("Avg",tstack,null);
		//calculate the derivative and get it's maximum
		float[] derivative=profile.clone();
		int maxpos=0;
		float maxder=0.0f;
		for(int i=1;i<profile.length;i++){
			derivative[i]=profile[i]-profile[i-1];
			if(derivative[i]>maxder){
				maxder=derivative[i];
				maxpos=i;
			}
		}
		derivative[0]=derivative[1];
		float[] secder=derivative.clone();
		for(int i=1;i<profile.length;i++){
			secder[i]=derivative[i]-derivative[i-1];
		}
		secder[0]=secder[1];
		if(outder) new PlotWindow4("Derivative Plot","Frame","Derivative",derivative).draw();
		//new PlotWindow4("2nd Derivative Plot","Frame","Second Derivative",secder).draw();
		//the start is one before the max position (counted starting at one)
		int bleachstart=maxpos-1;
		IJ.log("bleach start = "+maxpos);
		//now need to see where the derivative goes back to steady state
		//int maxduration=10;
		//start at the maximum flash duration
		float ssder=derivative[maxpos+maxduration-1];
		float thresh=ssder-enderthresh*maxder;
		IJ.log("steady state derivative = "+ssder);
		IJ.log("thresh = "+thresh);
		int bleachend=(maxpos+maxduration-1);
		for(int i=(maxpos+maxduration-1);i>(maxpos-1);i--){
			if(derivative[i]<thresh){
				bleachend=i-1;
				break;
			}
		}
		IJ.log("bleach end = "+(bleachend+1));
		//now perform the deletions
		Object[] stack2=jutils.stack2array(stack);
		int delframes=bleachend-bleachstart+1;
		int deltotframes=delframes*nchan*nslices;
		Object[] stack3=new Object[stack2.length-deltotframes];
		for(int i=0;i<bleachstart*nchan*nslices;i++){
			stack3[i]=stack2[i];
		}
		for(int i=(bleachend+1)*nchan*nslices;i<stack2.length;i++){
			stack3[i-deltotframes]=stack2[i];
		}
		jutils.create_hyperstack("Flash_Removed",jutils.array2stack(stack3,width,height),imp,nframes-delframes,nslices,nchan).show();
	}

}
