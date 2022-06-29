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
import jalgs.*;
import java.io.*;
import jguis.*;

public class average_trajectories_jru_v1 implements PlugIn {
	int index,histbins;
	float histstart,histend;

	public void run(String arg) {
		init_options();
		String[] options=jstatistics.stats;
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Statistic?",options,options[index]);
		gd.addCheckbox("Errors?",false);
		gd.addChoice("Error_Stat?",options,options[8]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		index=gd.getNextChoiceIndex();
		boolean showerr=gd.getNextBoolean();
		int index2=gd.getNextChoiceIndex();
		String statistic=options[index];
		String errstat=options[index2];
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[] xvals=((float[][])jutils.runPW4VoidMethod(iw,"getXValues"))[0];
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[] newxvals=new float[npts[0]];
		System.arraycopy(xvals,0,newxvals,0,npts[0]);
		float[] avg=new float[npts[0]];
		float[][] err=new float[1][];
		if(showerr) err[0]=new float[npts[0]];
		int length=yvals.length;
		float[] histoptions=null;
		if(statistic=="Mode"){
			histoptions=gethistoptions();
			if(histoptions==null){return;}
		}
		for(int i=0;i<npts[0];i++){
			float[] temp=new float[length];
			for(int j=0;j<length;j++){
				temp[j]=yvals[j][i];
			}
			avg[i]=jstatistics.getstatistic(statistic,temp,histoptions);
			if(showerr) err[0][i]=jstatistics.getstatistic(errstat,temp,histoptions);
		}
		String[] labels=(String[])jutils.runPW4VoidMethod(iw,"getAllLabels");
		PlotWindow4 pw=new PlotWindow4(labels[0]+"-"+statistic,labels[1],labels[2],newxvals,avg);
		if(showerr) pw.addErrors(err);
		pw.draw();
		set_options();
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
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"average_trajectories_jru_v1.jrn");
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
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"average_trajectories_jru_v1.jrn");
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
