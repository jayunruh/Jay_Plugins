/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class get_KISS_spectra_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,8,new String[]{"488laser_488dye","488laser_532dye","561laser_532dye","561laser_550dye","561laser_590dye","561laser_Cy5dye","633laser_590dye","633laser_Cy5dye"});
		if(imps==null || imps.length<8) return;
		ijpluginrunner ipr=new ijpluginrunner("create_spectrum_jru_v1");
		//start with the 488 spectra
		float[][] allspectra=new float[8][];
		float[][] allxvals=new float[8][];
		for(int i=0;i<8;i++){
			allspectra[i]=get_spectrum(imps[i].getStack(),imps[i].getRoi());
			float integral=jstatistics.getstatistic("Sum",allspectra[i],null);
			for(int j=0;j<allspectra[i].length;j++) allspectra[i][j]/=integral;
			allxvals[i]=(float[])ipr.runPluginMethod("get_xvals",new Object[]{imps[i]});
		}
		new PlotWindow4("488ex_spectra","slice","intensity",new float[][]{allxvals[0],allxvals[0]},new float[][]{allspectra[0],allspectra[1]},null).draw();
		new PlotWindow4("561ex_spectra","slice","intensity",new float[][]{allxvals[2],allxvals[2],allxvals[2],allxvals[2]},new float[][]{allspectra[2],allspectra[3],allspectra[4],allspectra[5]},null).draw();
		new PlotWindow4("633ex_spectra","slice","intensity",new float[][]{allxvals[6],allxvals[6]},new float[][]{allspectra[6],allspectra[7]},null).draw();
	}

	public float[] get_spectrum(ImageStack stack,Roi roi){
		Object[] arr=jutils.stack2array(stack);
		float[] spectrum=new float[arr.length];
		if(roi==null) for(int i=0;i<spectrum.length;i++) spectrum[i]=jstatistics.getstatistic("Sum",arr[i],null);
		else{
			boolean[] mask=jutils.roi2mask(roi,stack.getWidth(),stack.getHeight());
			for(int i=0;i<spectrum.length;i++) spectrum[i]=jstatistics.getstatistic("Sum",arr[i],stack.getWidth(),stack.getHeight(),mask,null);
		}
		return spectrum;
	}

}
