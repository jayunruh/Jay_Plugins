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

public class import_pt2_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=100000.0;
		gd.addNumericField("Sampling Frequency (Hz)?",sfreq,1,10,null);
		gd.addCheckbox("Output Spacings",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=(double)gd.getNextNumber();
		boolean outspace=gd.getNextBoolean();
		sfreq=2.5e11/(double)((int)(2.5e11/sfreq));
		//sfreq=9.766e8/(double)((int)(9.766e8/sfreq));
		double xinc=1.0/sfreq;
		//IJ.log(""+xinc);
		OpenDialog od = new OpenDialog("Open File","",".pt2");
		String directory=od.getDirectory();
		String name=od.getFileName();
		long tottime=0L;
		long ofl=0L;
		jdataio jdio=new jdataio();
		try{
			File infile=new File(directory+name);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			jdio.skipstreambytes(instream,16+6+18+12+18+2+256);
			jdio.skipstreambytes(instream,2*4);
			int rtchs=jdio.readintelint(instream);
			//IJ.log("router channels = "+rtchs);
			jdio.skipstreambytes(instream,22*4+20+16*4+9*4-3*4); //binhdr
			jdio.skipstreambytes(instream,16+8+33*4); //boardhdr
			jdio.skipstreambytes(instream,7*4); //tttrhdr
			int nrecords=jdio.readintelint(instream);
			//IJ.log("nrecords = "+nrecords);
			int imghdrsz=jdio.readintelint(instream);
			//IJ.log("imghdrsz = "+imghdrsz);
			jdio.skipstreambytes(instream,imghdrsz*4); //imghdr

			byte[] data=new byte[nrecords*4];
			long[][] pmdata=new long[rtchs][nrecords];
			if(!jdio.readintelbytefile(instream,nrecords*4,data)){
				instream.close();
				return;
			}
			instream.close();
			int[] counter=new int[rtchs];
			for(int i=0;i<nrecords;i++){
				byte[] temp={data[4*i],data[4*i+1],data[4*i+2],data[4*i+3]};
				int[] results=read_t3r_record(temp);
				if(results[2]==1){ //marker record
					if(results[3]==1){ //overflow record
						ofl+=210698240L;
					}
				} else {
					if(results[1]>=rtchs){IJ.log("Illegal Channel");}
					else {
						pmdata[results[1]][counter[results[1]]]=(long)results[0]+ofl;
						counter[results[1]]++;
					}
				}
				IJ.showProgress(i,nrecords);
			}
			long[][] pmdata2=new long[rtchs][];
			float[][] pmdataf=null;
			if(outspace) pmdataf=new float[rtchs][];
			for(int i=0;i<rtchs;i++){
				if(counter[i]>0){
					pmdata2[i]=new long[counter[i]];
					//pmdata2[i][0]=(pmdata[i][0]>>8);
					pmdata2[i][0]=pmdata[i][0];
					for(int j=1;j<counter[i];j++){
						//pmdata2[i][j]=((pmdata[i][j]-pmdata[i][j-1])>>8);
						pmdata2[i][j]=(pmdata[i][j]-pmdata[i][j-1]);
					}
					if(outspace){
						pmdataf[i]=new float[counter[i]];
						for(int j=0;j<counter[i];j++){
							pmdataf[i][j]=(float)pmdata2[i][j];
						}
					}
				}
			}
			if(outspace) new PlotWindow4("pmdata","photon","spacing",pmdataf,counter).draw();
			float[][] tmdata=new float[rtchs][];
			for(int i=0;i<rtchs;i++){
				if(counter[i]>0) tmdata[i]=(new pmodeconvert()).pm2tm(pmdata2[i],sfreq,250000000000L);
				//if(counter[i]>0) tmdata[i]=(new pmodeconvert()).pm2tm(pmdata2[i],sfreq,976562500L);
			}
			float[][] xvals=new float[rtchs][];

			for(int i=0;i<rtchs;i++){
				if(counter[i]>0){
					xvals[i]=new float[tmdata[i].length];
					for(int j=0;j<tmdata[i].length;j++) xvals[i][j]=(float)(xinc*(double)j);
				}
			}
			PlotWindow4 pw=new PlotWindow4("pt2 file","time (sec)","intensity",xvals[0],tmdata[0]);
			pw.draw();
			for(int i=1;i<rtchs;i++){
				if(counter[i]>0) pw.addPoints(xvals[i],tmdata[i],true);
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

	int[] read_t3r_record(byte[] record){
		int[] output=new int[4];
		output[0]=((record[3]&0xf)<<24) | ((record[2]&0xff)<<16) | ((record[1]&0xff)<<8) | (record[0]&0xff);
		output[1]=(record[3]&0xf0);
		if(output[1]==0xf0){
			//this is a special record
			output[2]=1;
			int markers=output[0]&0xf0;
			if(markers==0){
				//this is an overflow record
				output[3]=1;
			}
		}
		return output;
	}

}
