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
import jalgs.jfft.*;

public class calc_c3_wlswitch_freq_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open File","",".raw");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		int[] pmdata=null;
		jdataio jdio=new jdataio();
		try{
			File infile=new File(directory+name);
			int length=(int)(((float)infile.length()-128.0f)/4.0f);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			jdio.skipstreambytes(instream,128);
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

		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Center Freq (Hz)",10000.0,5,15,null);
		gd.addNumericField("freq resolution",0.001,5,15,null);
		gd.addNumericField("# of frequencies",100,0);
		gd.addNumericField("Photon Mode Freq",20000000,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		double sfreq=gd.getNextNumber();
		double freqres=gd.getNextNumber();
		int nfreq=(int)gd.getNextNumber();
		int pmfreq=(int)gd.getNextNumber();
		double centswfreq=20000000.0/1000.0111; //this is for a switching frequency of 20000.0
		double divider=20000.0/sfreq;
		centswfreq/=divider;
		sfreq=centswfreq-freqres*(double)nfreq/2.0;
		float[] ratio=new float[nfreq];
		float[] freqs=new float[nfreq];
		float[] hists=new float[nfreq*100];
		IJ.log("Switching Frequency (Hz) , Ratio(on/off)");
		double maxfreq=sfreq*1.0000111;
		float maxratio=0.0f;
		for(int i=0;i<nfreq;i++){
			double currswfreq=sfreq+freqres*(double)i;
			int offset=(new pmodeconvert()).calc_wlswitch_offset(pmdata,currswfreq,pmfreq);
			int[] newpmdata=(new pmodeconvert()).subpmoffset(pmdata,offset);
			float[] hist=(new pmodeconvert()).wlswitch_phase_hist(newpmdata,currswfreq,pmfreq);
			int halflength=(int)(0.5*(double)hist.length);
			float temp=0.0f;
			for(int j=0;j<halflength;j++) temp+=hist[j];
			float temp2=0.0f;
			for(int j=halflength;j<hist.length;j++) temp2+=hist[j];
			ratio[i]=temp/temp2;
			if(ratio[i]>maxratio){
				maxratio=ratio[i];
				maxfreq=currswfreq*1.0000111;
			}
			freqs[i]=(float)(currswfreq-centswfreq);
			for(int j=0;j<100;j++){
				hists[j+i*100]=hist[j];
			}
			IJ.log(""+currswfreq*1.0000111+" , "+ratio[i]);
			IJ.showProgress(i,nfreq);
			if(IJ.escapePressed()) break;
		}
		IJ.log("Max Freq = "+maxfreq+" , Max Ratio = "+maxratio);
		IJ.log("Uncorr Max Freq = "+maxfreq/1.0000111);
		if(maxratio<50.0) IJ.log("warning: 50:1 ratio was not exceeded");
		(new PlotWindow4("frequency optimization","freq difference","ratio",freqs,ratio)).draw();
		new ImagePlus("Hist carpet",new FloatProcessor(100,nfreq,hists,null)).show();
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
