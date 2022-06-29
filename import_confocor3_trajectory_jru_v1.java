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
import ij.io.*;
import java.io.*;
import jguis.*;
import jalgs.*;

public class import_confocor3_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=50000.0;
		//sfreq=20000000.0/500.5;
		gd.addNumericField("Sampling Frequency (Hz)?",sfreq,1,10,null);
		gd.addCheckbox("Photon_Mode",false);
		gd.addNumericField("Photon_Mode_Frequency",20000000,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=(double)gd.getNextNumber();
		boolean pmode=gd.getNextBoolean();
		int colfreq=(int)gd.getNextNumber();
		//IJ.log(""+colfreq);
		OpenDialog od = new OpenDialog("Open File","",".raw");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		int[] pmdata=null;
		jdataio jdio=new jdataio();
		//int colfreq=20000000;
		try{
			File infile=new File(directory+name);
			int length=(int)(((float)infile.length()-128.0f)/4.0f);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			jdio.skipstreambytes(instream,12);
			//colfreq=jdio.readintelint(instream);
			jdio.skipstreambytes(instream,128-12);
			pmdata=new int[length];
			if(!jdio.readintelintfile(instream,length,pmdata)){
				instream.close();
				return;
			}
			instream.close();
		}
		catch(IOException e){
			showErrorMessage(e);
			//return;
		}
		if(pmode){
			float[] pmdata2=algutils.convert_arr_float(pmdata);
			(new PlotWindow4(name,"photon","clock spacing",pmdata2)).draw();
		} else {
			float[] data=(new pmodeconvert()).pm2tm(pmdata,sfreq,colfreq);
			float[] xvals=new float[data.length];
			for(int i=0;i<data.length;i++){
				xvals[i]=(float)((double)i/(double)sfreq);
			}
			(new PlotWindow4(name,"x","I",xvals,data)).draw();
		}
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
