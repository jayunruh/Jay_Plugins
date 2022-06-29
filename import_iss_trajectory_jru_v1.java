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

public class import_iss_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open File","",".fcs");
		String directory=od.getDirectory();
		String name=od.getFileName();
		try{
			File infile=new File(directory+name);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			jdataio jdio=new jdataio();
			jdio.skipstreambytes(instream,1);
			String type=jdio.readstring(instream,1);
			int nch=1;
			boolean pmode=Character.isLowerCase(type.charAt(0));
			if(type.toLowerCase().equals("x")) nch=2;
			double sfreq=(double)jdio.readintelint(instream);
			int pmfreq=jdio.readintelint(instream);
			int dtype=jdio.readintelbyte(instream);
			int length=(int)infile.length();
			length-=256;
			if(dtype==0) length/=2;
			else length/=4;
			jdio.skipstreambytes(instream,256-11);
			float[][] traj=new float[nch][];
			if(pmode){
				GenericDialog gd=new GenericDialog("Options");
				gd.addNumericField("Sampling_Freq(Hz)",50000.0,5,15,null);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				sfreq=gd.getNextNumber();
				int[][] pmdata=new int[1][length];
				if(dtype==0) jdio.readintelshortfile(instream,length,pmdata[0]);
				else jdio.readintelintfile(instream,length,pmdata[0]);
				if(nch>1){
					int[][] pmdata2=new int[2][length/2];
					for(int i=0;i<length/2;i++){
						pmdata2[0][i]=pmdata[0][2*i];
						pmdata2[1][i]=pmdata[0][2*i+1];
					}
					pmdata=pmdata2;
				}
				traj[0]=(new pmodeconvert()).pm2tm(pmdata[0],sfreq,pmfreq);
				if(nch>1) traj[1]=(new pmodeconvert()).pm2tm(pmdata[1],sfreq,pmfreq);
			} else {
				traj[0]=new float[length];
				if(dtype==0) jdio.readintelshortfile(instream,length,traj[0]);
				else jdio.readintelintfile(instream,length,traj[0]);
				if(nch>1){
					float[][] tmdata2=new float[2][length/2];
					for(int i=0;i<length/2;i++){
						tmdata2[0][i]=traj[0][2*i];
						tmdata2[1][i]=traj[0][2*i+1];
					}
					traj=tmdata2;
				}
			}
			float[][] xvals=new float[nch][traj[0].length];
			double xinc=1.0/sfreq;
			for(int i=0;i<traj[0].length;i++){
				xvals[0][i]=(float)(xinc*(double)i);
			}
			if(nch>1) xvals[1]=xvals[0];

			PlotWindow4 pw=new PlotWindow4("fcs file","time (sec)","intensity",xvals,traj,null);
			pw.draw();
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
