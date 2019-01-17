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
import ij.text.*;
import ij.util.*;
import java.io.*;
import jalgs.*;
import jguis.*;

public class create_spectrum_jru_v1 implements PlugIn {
	//this plugin calculates some statistic for each image in a stack (or region of interest) and plots it
	int index,histbins;
	float histstart,histend;

	public void run(String arg) {
		init_options();
		String[] stats=jstatistics.stats;
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Spectrum Statistic?",stats,stats[index]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		index=gd.getNextChoiceIndex();
		String stat=stats[index];
		//get the image and its info
		ImagePlus imp = WindowManager.getCurrentImage();
		int slices=imp.getStack().getSize();
		Roi roi = imp.getRoi();
		float[] histoptions=jutils.getStatsOptions(stat);
		float[] spectral_data=get_spectrum(imp,roi,stat,histoptions);
		float[] slice_array=get_xvals(imp);
		set_options();
		//plot the spectrum
		if(slices==1){IJ.log(stat+" = "+spectral_data[0]);}
		else {
			PlotWindow4 plot = new PlotWindow4(stat+"_Spectrum","slice","intensity",slice_array,spectral_data);
			plot.draw();
		}
	}

	public static float[] get_xvals(ImagePlus imp){
		ImageStack stack = imp.getStack();
		int slices = stack.getSize();
		boolean lambdastack=is_lambda_stack(stack);
		float[] xvals=new float[slices];
		for(int i=0;i<slices;i++){
			if(lambdastack) xvals[i]=(float)get_label_lambda(stack.getSliceLabel(i+1));
			else xvals[i]=(float)(i+1);
		}
		return xvals;
	}

	public static float[] get_spectrum(ImagePlus imp,Roi roi,String stat,float[] extras){
		boolean rect=false;
		Rectangle r=null;
		int width=imp.getWidth();
		int height=imp.getHeight();
		if(roi!=null){
			r=roi.getBounds();
			if(roi.getType()==0) rect=true;
		}
		ImageStack stack = imp.getStack();
		int slices = stack.getSize();
		float[] spectral_data=new float[slices];
		if(roi==null || rect){
			for(int i=0;i<slices;i++){
				Object pixels=stack.getPixels(i+1);
				spectral_data[i]=jstatistics.getstatistic(stat,pixels,width,height,r,extras);
			}
		} else {
			boolean[] mask=jutils.roi2mask(roi,width,height);
			for(int i=0;i<slices;i++){
				Object pixels=stack.getPixels(i+1);
				spectral_data[i]=jstatistics.getstatistic(stat,pixels,width,height,mask,extras);
			}
		}
		return spectral_data;
	}

	public static boolean is_lambda_stack(ImageStack stack){
		if(stack.getSliceLabel(1)==null){return false;}
		for(int i=0;i<stack.getSize();i++){
			String label=stack.getSliceLabel(1);
			if(label.length()!=3){
				return false;
			}
			if(!Character.isDigit(label.charAt(0))){return false;}
			if(!Character.isDigit(label.charAt(1))){return false;}
			if(!Character.isDigit(label.charAt(2))){return false;}
		}
		if(get_label_lambda(stack.getSliceLabel(1))==get_label_lambda(stack.getSliceLabel(2))){return false;}
		return true;
	}

	public static int get_label_lambda(String label){
		return Integer.parseInt(label);
	}

	private float[] gethistoptions(){
		GenericDialog gd=new GenericDialog("Histogram Options");
		gd.addNumericField("Histogram Bins",histbins,0);
		gd.addNumericField("Histogram Start",histstart,5,10,null);
		gd.addNumericField("Histogram End",histend,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		float[] temp=new float[3];
		histbins=(int)gd.getNextNumber();
		histstart=(float)gd.getNextNumber();
		histend=(float)gd.getNextNumber();
		temp[0]=(float)histbins; temp[1]=histstart; temp[2]=histend;
		return temp;
	}

	void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"create_spectrum_jru_v1.jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			index=Integer.parseInt(d.readLine());
			histbins=Integer.parseInt(d.readLine());
			histstart=Float.parseFloat(d.readLine());
			histend=Float.parseFloat(d.readLine());
			d.close();
		}
		catch(IOException e){
			index=0;
			histbins=100;
			histstart=2450.0f;
			histend=2550.0f;
			set_options();
		}
		return;
	}
	
	void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){a.mkdir();}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"create_spectrum_jru_v1.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+index+"\n");
			d.write(""+histbins+"\n");
			d.write(""+histstart+"\n");
			d.write(""+histend+"\n");
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

}
