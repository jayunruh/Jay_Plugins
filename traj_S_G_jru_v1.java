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

public class traj_S_G_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin calculates the real and imaginary components of the FFT at harmonic
		//for a lifetime trajectory
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
		ImageWindow iw=WindowManager.getCurrentWindow();
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		if(sel<0) sel=0;
		float[] data=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[sel];
		slices=data.length;
		
		//calculate the fourier transform at each point
		//start by calculating the trig functions
		float[] cosvals = new float[slices];
		float[] sinvals = new float[slices];
		float phase_shift_fraction=phase_shift/360.0f;
		for(i=0;i<slices;i++){
			cosvals[i]=mod_scale*(float)Math.cos((float)harmonic*2.0f*Math.PI*((float)i/(float)slices-phase_shift_fraction));
			sinvals[i]=mod_scale*(float)Math.sin((float)harmonic*2.0f*Math.PI*((float)i/(float)slices-phase_shift_fraction));
		}

		float S=0.0f; float G=0.0f; float intensity=0.0f;
		for(j=0;j<slices;j++){
			intensity+=data[j];
			S+=sinvals[j]*data[j];
			G+=cosvals[j]*data[j];
		}
		if(intensity>0.0f){
			S/=intensity;
			G/=intensity;
		} else {
			S=0.0f;
			G=0.0f;
		}

		float phase=(float)((Math.atan(S/G)*180.0)/(Math.PI));
		if(G<0.0f){phase+=180.0f;}
		if(S<0.0f && G>0.0f){phase+=360.0f;}
		float mod=(float)Math.sqrt(G*G+S*S);

		IJ.log("Intensity = "+intensity);
		IJ.log("G = "+G);
		IJ.log("S = "+S);
		IJ.log("Phase = "+phase);
		IJ.log("Phase (rad) = "+phase*(float)Math.PI/180.0f);
		IJ.log("Mod = "+mod);
		
	}

}
