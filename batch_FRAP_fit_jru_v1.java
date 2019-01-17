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
import ij.plugin.frame.*;
import jguis.*;
import jalgs.jfit.*;
import ij.text.*;

public class batch_FRAP_fit_jru_v1 implements PlugIn {
	fit_exp fitfunc;
	boolean twocomp;

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] data=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int maxpts=data[0].length;
		int minnpts=npts[0];
		for(int i=1;i<npts.length;i++) minnpts=Math.min(minnpts,npts[i]);
		GenericDialog gd=new GenericDialog("Options");
		int bf=4;
		gd.addNumericField("Frames_before frap",bf,0);
		int totframes=maxpts;
		gd.addNumericField("Frames_to analyze",totframes,0);
		float mintau=2.0f;
		gd.addNumericField("Minimum tau (frames)",mintau,5,15,null);
		float maxtau=3.0f*(float)(maxpts-bf);
		gd.addNumericField("Maximum tau (frames)",maxtau,5,15,null);
		gd.addCheckbox("Two Component?",false);
		gd.addNumericField("Min_tau_ratio",2.0f,5,15,null);
		gd.addCheckbox("Output_chisquared",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		bf=(int)gd.getNextNumber();
		totframes=(int)gd.getNextNumber();
		mintau=(float)gd.getNextNumber();
		maxtau=(float)gd.getNextNumber();
		twocomp=gd.getNextBoolean();
		float minratio=(float)gd.getNextNumber();
		boolean outchi=gd.getNextBoolean();
		fitfunc=new fit_exp(mintau,maxtau);
		fitfunc.minratio=minratio;
		String labels="Prefrap\tBaseline\tAmplitude\tTau(frames)";
		if(twocomp) labels="Prefrap\tBaseline\tAmplitude1\tTau1(frames)\tAmplitude2\tTau2(frames)";
		if(outchi) labels+="\tc2";
		TextWindow tw=new TextWindow("Fit Results",labels,"",400,200);
		Plot4[] plots=new Plot4[data.length];
		for(int i=0;i<data.length;i++){
			int tempint=totframes;
			if(npts[i]<totframes){tempint=npts[i];}
			float[] temp=fit_frap(data[i],bf,tempint);
			StringBuffer sb=new StringBuffer();
			if(!twocomp){
				sb.append(""+temp[0]+"\t"+temp[1]+"\t"+temp[2]+"\t"+temp[3]);
				if(outchi) sb.append("\t"+temp[4]+"\n");
				else sb.append("\n");
			}else{
				sb.append(""+temp[0]+"\t"+temp[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[4]+"\t"+temp[5]);
				if(outchi) sb.append("\t"+temp[6]+"\n");
				else sb.append("\n");
			}
			tw.append(sb.toString());
			plots[i]=new Plot4("Frame","Intensity",get_data_fit(data[i],bf,tempint,temp),null);
		}
		if(plots.length==1){
			new PlotWindow4("Fits",plots[0]).draw();
		} else {
			new PlotStack4("Fits",plots).draw();
		}
	}

	public float[][] get_data_fit(float[] data,int bf,int totpoints,float[] params){
		float[][] temp=new float[2][totpoints];
		System.arraycopy(data,0,temp[0],0,totpoints);
		for(int i=0;i<bf;i++){
			temp[1][i]=params[0];
		}
		for(int i=0;i<(totpoints-bf);i++){
			temp[1][i+bf]=(float)((params[1]+params[2])-params[2]*Math.exp(-(double)i/params[3]));
			if(twocomp) temp[1][i+bf]-=(float)(params[4]*Math.exp(-(double)i/params[5])-params[4]);
		}
		return temp;
	}

	public float[] fit_frap(float[] data,int bf,int totpoints){
		float[] results=new float[7];
		for(int i=0;i<bf;i++){
			results[0]+=data[i]/(float)bf;
		}
		float[] newdata=new float[totpoints-bf];
		System.arraycopy(data,bf,newdata,0,totpoints-bf);
		float[] params=null;
		if(twocomp) params=fit_rising_2exp(newdata);
		else params=fit_rising_exp(newdata);
		System.arraycopy(params,0,results,1,params.length);
		return results;
	}

	public float[] fit_rising_exp(float[] data){
		double[] params=fitfunc.fitdata(data);
		float[] retvals={(float)(params[0]+params[1]),-(float)params[1],(float)params[2],(float)params[3]};
		return retvals;
	}

	public float[] fit_rising_2exp(float[] data){
		double[] params=fitfunc.fitdata2exp(data);
		float[] retvals={(float)(params[0]+params[1]+params[3]),-(float)params[1],(float)params[2],-(float)params[3],(float)params[4],(float)params[5]};
		return retvals;
	}

}
