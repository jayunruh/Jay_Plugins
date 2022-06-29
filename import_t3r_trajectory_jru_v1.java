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

public class import_t3r_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=20000.0;
		gd.addNumericField("Sampling Frequency (Hz)?",sfreq,1,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=(double)gd.getNextNumber();
		sfreq=10000000.0/(double)((int)(10000000.0/sfreq));
		OpenDialog od = new OpenDialog("Open File","",".t3r");
		String directory=od.getDirectory();
		String name=od.getFileName();
		int tottime=0;
		int ofl=0;
		int counter=0;
		int counter2=0;
		try{
			File infile=new File(directory+name);
			InputStream instream=new BufferedInputStream(new FileInputStream(infile));
			skipstreambytes(instream,540);
			int clockns=readintelintstream(instream);
			//IJ.log("clock ns = "+clockns);
			skipstreambytes(instream,40);
			int nrecords=readintelintstream(instream);
			//IJ.log("n records = "+nrecords);
			int skiplength=4*readintelintstream(instream);
			if(skiplength>0){skipstreambytes(instream,skiplength);}
			byte[] data=new byte[nrecords*4];
			int[] pmdata=new int[nrecords];
			int[] pmdata2=new int[nrecords];
			if(!readintelbytefile(instream,nrecords,data)){
				instream.close();
				return;
			}
			instream.close();
			for(int i=0;i<nrecords;i++){
				byte[] temp={data[4*i],data[4*i+1],data[4*i+2],data[4*i+3]};
				int[] results=read_t3r_record2(temp);
				if(results[1]==0){
					if(results[3]==1){
						ofl+=65536;
					}
				} else {
					if(results[2]==0){
						pmdata[counter]=results[4]+ofl;
						counter++;
					} else{
						pmdata2[counter2]=results[4]+ofl;
						counter2++;
					}
				}
			}
			int[] pmdata3=new int[counter];
			int[] pmdata4=new int[counter2];
			pmdata3[0]=pmdata[0];
			for(int i=1;i<counter;i++){
				pmdata3[i]=pmdata[i]-pmdata[i-1];
			}
			pmdata4[0]=pmdata2[0];
			for(int i=1;i<counter2;i++){
				pmdata4[i]=pmdata2[i]-pmdata2[i-1];
			}
			float[] tmdata=(new pmodeconvert()).pm2tm(pmdata3,sfreq,10000000);
			float[] tmdata2=(new pmodeconvert()).pm2tm(pmdata4,sfreq,10000000);
			float[] xvals=new float[tmdata.length];
			double xinc=1.0/sfreq;
			for(int i=0;i<tmdata.length;i++){
				xvals[i]=(float)(xinc*(double)i);
			}
			float[] xvals2=new float[tmdata2.length];
			for(int i=0;i<tmdata2.length;i++){
				xvals2[i]=(float)(xinc*(double)i);
			}

			PlotWindow4 pw=new PlotWindow4("t3r file","time (sec)","intensity",xvals,tmdata);
			pw.draw();
			if(tmdata2!=null){
				pw.addPoints(xvals2,tmdata2,true);
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

	private boolean readintelbytefile(InputStream instream,int nrecords,byte[] data){
		int dumint;
		//data=new byte[nrecords*4];
		try{dumint=instream.read(data);}
		catch(IOException e){
			showErrorMessage(e);
			return false;
		}
		if(dumint<0){return false;}
                        return true;
	}

	int readintelintstream(InputStream instream){
		byte[] record=new byte[4];
		int dumint;
		try{dumint=instream.read(record);}
		catch(IOException e){
			showErrorMessage(e);
			return -1;
		}
		if(dumint<0){return -1;}
		return (((record[3]&0xff)<<24) | ((record[2]&0xff)<<16) | ((record[1]&0xff)<<8) | (record[0]&0xff));
	}

	void skipstreambytes(InputStream instream, int skip){
		int left;
		try{
			left=skip;
			do{
				int temp=(int)instream.skip(left);
				left-=temp;
			}while(left>0);
		}
		catch(IOException e){
			showErrorMessage(e);
			IJ.showMessage("Error reading data");
			return;
		}
                        return;
	}

	int[] read_t3r_record(byte[] record){
		int reserved=((record[0]&0x80)>>7);
		int valid=((record[0]&0x40)>>6);
		int route=((record[0]&0x30)>>4);
		if(valid==1){
			int channel=((record[0]&0x0f)<<8) | (record[1]&0xff);
			int timetag=((record[2]&0xff)<<8) | (record[3]&0xff);
			int[] temp={reserved,valid,route,channel,timetag};
			return temp;
		} else {
			int overflow=(record[0]&0x08)>>3;
			int timetag=(int)(((record[2]&0xff)<<8) | (record[3]&0xff));
			int[] temp={reserved,valid,route,overflow,timetag};
			return temp;
		}
	}

	int[] read_t3r_record2(byte[] record){
		int reserved=(record[3]&0x80)>>7;
		int valid=(record[3]&0x40)>>6;
		int route=((record[3]&0x30)>>4);
		if(valid==1){
			int channel=((record[3]&0x0f)<<8) | (record[2]&0xff);
			int timetag=((record[1]&0xff)<<8) | (record[0]&0xff);
			int[] temp={reserved,valid,route,channel,timetag};
			return temp;
		} else {
			int overflow=(record[3]&0x08)>>3;
			int timetag=(int)(((record[1]&0xff)<<8) | (record[0]&0xff));
			int[] temp={reserved,valid,route,overflow,timetag};
			return temp;
		}
	}

}
