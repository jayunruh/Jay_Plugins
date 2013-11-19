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

public class import_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] typeoptions={"Byte","Short","Integer","Float","Double"};
		gd.addChoice("File Type?",typeoptions,typeoptions[1]);
		int offset=0;
		gd.addNumericField("Offset (bytes)?",offset,0);
		double sfreq=10000.0;
		gd.addNumericField("Sampling Frequency (Hz)",sfreq,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int typeindex=gd.getNextChoiceIndex();
		offset=(int)gd.getNextNumber();
		sfreq=gd.getNextNumber();
		OpenDialog od = new OpenDialog("Open File","","");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		int databytes=typeindex+1;
		if(typeindex>1){databytes=4;}
		if(typeindex==4){databytes=8;}
		jdataio jdio=new jdataio();
		try{
			File infile=new File(directory+name);
			int length=(int)((float)(infile.length()-offset)/(float)databytes);
			float[] xvals=new float[length];
			for(int i=0;i<length;i++){
				xvals[i]=(float)i/(float)sfreq;
			}
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			float[] data=new float[length];
			jdio.skipstreambytes(instream,offset);
			if(typeindex==0){
				if(jdio.readintelbytefile(instream,length,data)){(new PlotWindow4(name,"time (s)","Intensity",xvals,data)).draw();}
			} else {
				if(typeindex==1){
					if(jdio.readintelshortfile(instream,length,data)){(new PlotWindow4(name,"time (s)","Intensity",xvals,data)).draw();}
				} else {
					if(typeindex==2){
						if(jdio.readintelintfile(instream,length,data)){(new PlotWindow4(name,"time (s)","Intensity",xvals,data)).draw();}
					} else {
						if(typeindex==3){
							if(jdio.readintelfloatfile(instream,length,data)){(new PlotWindow4(name,"time (s)","Intensity",xvals,data)).draw();}
						} else {
							if(jdio.readinteldoublefile(instream,length,data)){(new PlotWindow4(name,"time (s)","Intensity",xvals,data)).draw();}
						}
					}
				}
			}
			instream.close();
		} catch(IOException e){
			showErrorMessage(e);
			//return;
		}
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
