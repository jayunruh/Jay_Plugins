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

public class read_BH_sdt_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean forcekinetics=false;
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Force Kinetics",forcekinetics);
		gd2.addCheckbox("Output Info",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		forcekinetics=gd2.getNextBoolean();
		boolean output=gd2.getNextBoolean();
		OpenDialog od = new OpenDialog("Open BH sdt File",arg);
		String directory = od.getDirectory();
		String name = od.getFileName();
		if(output) IJ.log(name);
		ImageStack stack,stack2,stack3,stack4;
		int channels=1;
		try{
			InputStream instream = new BufferedInputStream(new FileInputStream(directory+name));
			byte dumbyte=0;
			//find the offset to the data block
			skipstreambytes(instream,14); int current=14;
			int dboffset=readintelintstream(instream); current+=4;
			if(output) IJ.log("Data Block Offset = "+dboffset);
			//find the offset to the measurement description block
			skipstreambytes(instream,2); current+=2;
			int dblength=readintelintstream(instream); current+=4;
			if(output) IJ.log("Data Block Length = "+dblength);
			int mdoffset=readintelintstream(instream); current+=4;
			if(output) IJ.log("Measurement Description Block Offset = "+mdoffset);
			skipstreambytes(instream,2); current+=2;
			int mdlength=(int)(readintelshortstream(instream)&0xffff); current+=2;
			if(output) IJ.log("Measurement Description Block Length = "+mdlength);
			//ffwd to the measurement description block
			skipstreambytes(instream,(mdoffset-current)); current+=(mdoffset-current);
			//now find the adc resolution
			skipstreambytes(instream,82); current+=82;
			int adcres=(int)(readintelshortstream(instream)&0xffff); current+=2;
			if(output) IJ.log("ADC resolution = "+adcres);
			skipstreambytes(instream,89); current+=89;
			if(output) IJ.log("current "+current);
			int width=readintelintstream(instream); current+=4;
			if(output) IJ.log("Width = "+width);
			int height=readintelintstream(instream); current+=4;
			if(output) IJ.log("Height = "+height);
			if((width+height)>0 && !forcekinetics){
				channels=dblength/(2*width*height*adcres);
			} else {
				channels=dblength/(2*adcres);
			}
			if(output) IJ.log("Channels = "+channels);
			skipstreambytes(instream,(dboffset-current));
			skipstreambytes(instream,22);
			if((width+height)>0 && !forcekinetics){
				stack=new ImageStack(width,height);
				stack2=new ImageStack(width,height);
				stack3=new ImageStack(width,height);
				stack4=new ImageStack(width,height);
				for(int i=0;i<adcres;i++){stack.addSlice("",new short[width*height]);}
				if(channels>1){
					for(int i=0;i<adcres;i++){stack2.addSlice("",new short[width*height]);}
					if(channels>2){
						for(int i=0;i<adcres;i++){stack3.addSlice("",new short[width*height]);}
						if(channels>3){for(int i=0;i<adcres;i++){stack4.addSlice("",new short[width*height]);}}
					}
				}
				for(int i=0;i<width*height;i++){
					for(int j=0;j<adcres;j++){
						((short[])stack.getPixels(j+1))[i]=readintelshortstream(instream);
					}
					IJ.showProgress(i,width*height);
				}
				if(channels>1){
					for(int i=0;i<width*height;i++){
						for(int j=0;j<adcres;j++){
								((short[])stack2.getPixels(j+1))[i]=readintelshortstream(instream);
						}
						IJ.showProgress(i,width*height);
					}
					if(channels>2){
						for(int i=0;i<width*height;i++){
							for(int j=0;j<adcres;j++){
									((short[])stack3.getPixels(j+1))[i]=readintelshortstream(instream);
							}
							IJ.showProgress(i,width*height);
						}
						if(channels>3){
							for(int i=0;i<width*height;i++){
								for(int j=0;j<adcres;j++){
										((short[])stack4.getPixels(j+1))[i]=readintelshortstream(instream);
								}
								IJ.showProgress(i,width*height);
							}
						}
					}
				}
				(new ImagePlus(name,stack)).show();
				if(channels>1){
					(new ImagePlus(name+"ch2",stack2)).show();
					if(channels>2){
						(new ImagePlus(name+"ch3",stack3)).show();
						if(channels>3){
							(new ImagePlus(name+"ch4",stack4)).show();
						}
					}
				}
			} else {
				float[][] data=new float[channels][adcres];
				float[][] xvals=new float[channels][adcres];
				float ttimens=12.5f;
				GenericDialog gd=new GenericDialog("Options");
				gd.addNumericField("Total Time Window (ns)",ttimens,5,10,null);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				ttimens=(float)gd.getNextNumber();
				float nschan=(ttimens/(float)adcres);
				for(int i=0;i<channels;i++){
					for(int j=0;j<adcres;j++){
						data[i][j]=(float)(readintelshortstream(instream)&0xffff);
						xvals[i][j]=nschan*(float)j;
					}
				}
				new PlotWindow4(name,"time (ns)","Intensity",xvals,data,null).draw();
			}
			instream.close();
		}
		catch(IOException e){
			showErrorMessage(e);
			return;
		}

	}

	private float readintelfloatstream(InputStream instream){
		byte[] dumbyte=new byte[4];
		try{
			for(int i=0;i<4;i++){dumbyte[i]=(byte)instream.read();}
		}
		catch(IOException e){
			showErrorMessage(e);
			return 0.0f;
		}
		int tmp;
		tmp = (int)(((dumbyte[3]&0xff)<<24) | ((dumbyte[2]&0xff)<<16) | ((dumbyte[1]&0xff)<<8) | (dumbyte[0]&0xff));
		return Float.intBitsToFloat(tmp);
	}

	private short readintelshortstream(InputStream instream){
		byte[] dumbyte=new byte[2];
		try{
			for(int i=0;i<2;i++){dumbyte[i]=(byte)instream.read();}
		}
		catch(IOException e){
			showErrorMessage(e);
			return 0;
		}
                        return (short)(((dumbyte[1]&0xff)<<8) | (dumbyte[0]&0xff));
	}

	private int readintelintstream(InputStream instream){
		byte[] dumbyte=new byte[4];
		try{
			for(int i=0;i<4;i++){dumbyte[i]=(byte)instream.read();}
		}
		catch(IOException e){
			showErrorMessage(e);
			return 0;
		}
                        return (int)(((dumbyte[3]&0xff)<<24) | ((dumbyte[2]&0xff)<<16) | ((dumbyte[1]&0xff)<<8) | (dumbyte[0]&0xff));
	}

	private void skipstreambytes(InputStream instream, int skip){
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

	private void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured writing the file.\n \n" + msg);
	}
}
