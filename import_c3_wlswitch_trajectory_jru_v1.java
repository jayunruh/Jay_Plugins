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

public class import_c3_wlswitch_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Switching_Frequency(Hz)",10000.0,5,15,null);
		gd.addNumericField("Illumination_Delay(us)",2.5,5,15,null);
		gd.addNumericField("Bin_By",1,0);
		gd.addCheckbox("Output_histogram",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		double sfreq=gd.getNextNumber();
		double ill_delay=gd.getNextNumber();
		ill_delay*=1.0e-6;
		int binby=(int)gd.getNextNumber();
		boolean outhist=gd.getNextBoolean();
		OpenDialog od = new OpenDialog("Select_the_green_channel","",".raw");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		String ch=name.substring(name.length()-5,name.length()-4);
		String name2;
		if(ch.equals("1")){
			name2=name.substring(0,name.length()-5)+"2.raw";
		} else {
			name2=name.substring(0,name.length()-5)+"1.raw";
		}
		//IJ.log(name);
		//IJ.log(name2);
		int[] pmdata=null;
		int[] pmdata2=null;
		jdataio jdio=new jdataio();
		try{
			File infile=new File(directory+name);
			File infile2=new File(directory+name2);
			int length=(int)(((double)infile.length()-128.0)/4.0);
			int length2=(int)(((double)infile2.length()-128.0)/4.0);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			InputStream instream2=new BufferedInputStream(new FileInputStream(infile2));
			jdio.skipstreambytes(instream,128);
			jdio.skipstreambytes(instream2,128);
			pmdata=new int[length];
			pmdata2=new int[length2];
			if(!jdio.readintelintfile(instream,length,pmdata)){
				instream.close();
				instream2.close();
				return;
			}
			if(!jdio.readintelintfile(instream2,length2,pmdata2)){
				instream.close();
				instream2.close();
				return;
			}
			instream.close();
			instream2.close();
		}
		catch(IOException e){
			showErrorMessage(e);
			//return;
		}
		double swfreq=20000000.0/1000.0111;
		double divider=20000.0/sfreq;
		swfreq/=divider;
		int offset=(new pmodeconvert()).calc_wlswitch_offset(pmdata,swfreq,20000000);
		float[][] tmdata=(new pmodeconvert()).pm2tm_alex(pmdata,pmdata2,swfreq,20000000,offset,binby,ill_delay);
		float[][] tmdatax=new float[2][tmdata[0].length];
		float xinc=1.0f/(float)swfreq;
		for(int i=0;i<tmdata[0].length;i++){
			tmdatax[0][i]=xinc*(float)i;
			tmdatax[1][i]=xinc*(float)i;
		}
		(new PlotWindow4(name,"time bin","Intensity",tmdatax,tmdata,null)).draw();
		if(outhist){
			float[][] hist=new float[2][];
			int[] newpmdata=(new pmodeconvert()).subpmoffset(pmdata,offset);
			hist[0]=(new pmodeconvert()).wlswitch_phase_hist(newpmdata,swfreq,20000000,ill_delay);
			int[] newpmdata2=(new pmodeconvert()).subpmoffset(pmdata2,offset);
			hist[1]=(new pmodeconvert()).wlswitch_phase_hist(newpmdata2,swfreq,20000000,ill_delay);
			float[][] histx=new float[2][100];
			for(int i=0;i<100;i++){
				histx[0][i]=(float)i/100.0f;
				histx[1][i]=(float)i/100.0f;
			}
			(new PlotWindow4("Phase Histogram","phase bin","Occurences",histx,hist,null)).draw();
		}
	}

	void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured reading the file.\n \n" + msg);
	}

}
