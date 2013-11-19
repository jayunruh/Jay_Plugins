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

public class flim_S_G_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin calculates the real and imaginary components of the FFT at harmonic
		//for a FLIM image
		int i,height,width,slices,j,k,shifted,harmonic;
		FloatProcessor fp = null;
		
		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Phase Offset (degrees)",0.0,5,10,null);
		gd.addNumericField("Mod Scaling Factor",1.0,5,10,null);
		gd.addNumericField("Harmonic",1.0,0,10,null);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		float phase_shift=(float)gd.getNextNumber();
		float mod_scale=(float)gd.getNextNumber();
		harmonic=(int)gd.getNextNumber();

		//get the current image and its info
		ImagePlus imp = WindowManager.getCurrentImage();
		height = imp.getHeight();
		width = imp.getWidth();
		ImageStack stack = imp.getStack();
		slices = stack.getSize();
		fp=(FloatProcessor)imp.getProcessor();
		Rectangle r = fp.getRoi();

		float[] S = new float[r.width*r.height];
		float[] G = new float[r.width*r.height];
		float[] intensity = new float[r.width*r.height];
		
		//calculate the fourier transform at each point
		//start by calculating the trig functions
		float[] cosvals = new float[slices];
		float[] sinvals = new float[slices];
		float phase_shift_fraction=phase_shift/360.0f;
		for(i=0;i<slices;i++){
			cosvals[i]=mod_scale*(float)Math.cos((float)harmonic*2.0f*Math.PI*((float)i/(float)slices-phase_shift_fraction));
			sinvals[i]=mod_scale*(float)Math.sin((float)harmonic*2.0f*Math.PI*((float)i/(float)slices-phase_shift_fraction));
		}
		int row=0;
		int column=0;
		for(i=0;i<(r.width*r.height);i++){
			S[i]=0.0f; G[i]=0.0f; intensity[i]=0.0f;
			for(j=0;j<slices;j++){
				float[] temp = (float[])stack.getPixels(j+1);
				float data = temp[column+r.x+width*(row+r.y)];
				intensity[i]+=data;
				S[i]+=sinvals[j]*data;
				G[i]+=cosvals[j]*data;
			}
			if(intensity[i]>0.0f){
				S[i]/=intensity[i];
				G[i]/=intensity[i];
			} else {
				S[i]=0.0f;
				G[i]=0.0f;
			}
			column++;
			if(column>=r.width){column=0; row++;}
		}

		FloatProcessor fp2 = new FloatProcessor(r.width,r.height,intensity,null);
		ImagePlus imp2=new ImagePlus("Intensity",fp2);
		imp2.show();
		FloatProcessor fp3 = new FloatProcessor(r.width,r.height,S,null);
		ImagePlus imp3=new ImagePlus("S",fp3);
		imp3.show();
		FloatProcessor fp4 = new FloatProcessor(r.width,r.height,G,null);
		ImagePlus imp4=new ImagePlus("G",fp4);
		imp4.show();
	}

}
