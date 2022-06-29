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
import jalgs.*;
import jguis.*;

public class import_toa_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open File","",".toa");
		String directory=od.getDirectory();
		String name=od.getFileName();
		try{
			File infile=new File(directory+name);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			jdataio jdio=new jdataio();
			int nch=1;
			long pmfreq=(long)80000000;
			double sfreq=100000.0;
			int length=(int)infile.length();
			length/=4;
			float[][] traj=new float[nch][];
			GenericDialog gd=new GenericDialog("Options");
			gd.addNumericField("Sampling_Freq(Hz)",50000.0,5,15,null);
			gd.addCheckbox("Show Photon Mode Trace",false);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			sfreq=gd.getNextNumber();
			boolean show=gd.getNextBoolean();
			long[][] pmdata=new long[1][length];
			jdio.readintelintfile(instream,pmdata[0]);
			long prev=pmdata[0][0];
			for(int i=1;i<pmdata[0].length;i++){
				long temp=pmdata[0][i]-prev;
				//if(prev>0 && pmdata[0][i]<0) IJ.log(""+prev+" , "+pmdata[0][i]+" , "+temp);
				//if(temp<0) IJ.log(""+prev+" , "+pmdata[0][i]+" , "+temp);
				if(temp<=0){temp=1; IJ.log("0 or negative spacing found, truncating at 1");}
				prev=pmdata[0][i];
				pmdata[0][i]=temp;
			}
			if(show) new PlotWindow4("photon mode data","photon","time",algutils.convert_arr_float(pmdata[0])).draw();
			else {
				traj[0]=(new pmodeconvert()).pm2tm(pmdata[0],sfreq,pmfreq);
				float[][] xvals=new float[1][traj[0].length];
				double xinc=1.0/sfreq;
				for(int i=0;i<traj[0].length;i++){
					xvals[0][i]=(float)(xinc*(double)i);
				}

				PlotWindow4 pw=new PlotWindow4(name,"time (sec)","intensity",xvals[0],traj[0]);
				pw.draw();
			}
		}
		catch(IOException e){
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
