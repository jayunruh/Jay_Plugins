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
import java.io.*;
import jalgs.*;
import jalgs.jfft.*;
import jalgs.jfit.*;
import jguis.*;
import ij.io.*;
import java.text.*;
import ij.text.*;

public class analysis_auto_corr_v2 implements PlugIn {
	//this plugin is a gui for fitting correlation curves singly or globally for 
	//3D Gaussian or Gaussian Lorentzian squared point spread functions
	//(solution confocal, solution two photon)
	//copyright 2009 Jay Unruh, Stowers Institute for Medical Research
	float trajkhz,khz;
	int trajlength,binby;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=50000.0;
		gd.addNumericField("Sampling Frequency?",sfreq,1,10,null);
		String[] psfchoice={"3D_Gaussian","2D_Gaussian","2Dxz_Gaussian"};
		gd.addChoice("PSF Type?",psfchoice,psfchoice[0]);
		String[] filetypechoice={"Confocor 3 raw","Short binary trajectory","PlotWindow file","PlotWindow trajectory","Autocorr trajectory"};
		gd.addChoice("File Type?",filetypechoice,filetypechoice[0]);
		boolean showtraj=true;
		gd.addCheckbox("Show Trajectories?",showtraj);
		binby=20;
		gd.addNumericField("Traj Bin?",binby,0,10,null);
		boolean detrend=false;
		gd.addCheckbox("Segmented Linear Detrending?",detrend);
		int segments=2;
		gd.addNumericField("Detrending Segments?",segments,0);
		boolean brightcorr=false;
		gd.addCheckbox("Bright Corr?",brightcorr);
		boolean pad=false;
		gd.addCheckbox("Pad Trajectory?",pad);
		gd.addCheckbox("Simple Analysis?",false);
		int pmfreq=15000000;
		gd.addNumericField("Photon Mode Freq",pmfreq,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=gd.getNextNumber();
		int psfflag=gd.getNextChoiceIndex();
		int fileflag=gd.getNextChoiceIndex();
		showtraj=gd.getNextBoolean();
		binby=(int)gd.getNextNumber();
		detrend=gd.getNextBoolean();
		segments=(int)gd.getNextNumber();
		brightcorr=gd.getNextBoolean();
		pad=gd.getNextBoolean();
		boolean simple=gd.getNextBoolean();
		pmfreq=(int)gd.getNextNumber();
		trajkhz=((float)sfreq/(float)binby)/1000.0f;
		khz=(float)sfreq/1000.0f;
		int nfiles=0;
		Object[] correlations=null;
		float[][] trajectories=null;
		int xmax=0;
		int ymax=0;
		String[] names=null;
		boolean first=true;
		autocorr acclass=null;
		detrend_linear dl=null;
		int size=0;
		int newsize=0;
		float[] xvals=null;
		float[] avg=null;
		float[] var=null;
		trajlength=0;
		binmultilog bml=new binmultilog();
		kstats kstatsfunc=new kstats();
		if(fileflag<3){
			jdataio ioclass=new jdataio();
			File[] filearray=ioclass.openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
			if(filearray.length==0){return;}
			String dir=filearray[0].getAbsolutePath();
			int sepindex=dir.lastIndexOf(File.separator);
			String newdir=dir.substring(0,sepindex+1);
			OpenDialog.setDefaultDirectory(newdir);
			nfiles=filearray.length;
			correlations=new Object[nfiles];
			avg=new float[nfiles];
			var=new float[nfiles];
			names=new String[nfiles+1];
			names[nfiles]="avg";
			for(int i=0;i<nfiles;i++){
				try{
					names[i]=filearray[i].getName();
					int length1=(int)(((double)filearray[i].length()-128.0)/4.0);
					int length2=(int)(((double)filearray[i].length())/2.0);
					InputStream instream=new BufferedInputStream(new FileInputStream(filearray[i]));
					float[] tmdata;
					if(fileflag==0){
						int[] pmdata=new int[length1];
						if(!ioclass.skipstreambytes(instream,128)){showioerror(); instream.close(); return;}
						if(!ioclass.readintelintfile(instream,length1,pmdata)){showioerror(); instream.close(); return;}
						tmdata=(new pmodeconvert()).pm2tm(pmdata,sfreq,pmfreq);
					} else {
						if(fileflag==1){
							tmdata=new float[length2];
							if(!ioclass.readintelshortfile(instream,length2,tmdata)){showioerror(); instream.close(); return;}
						} else {
							Plot4 p4=new Plot4(filearray[i].getPath());
							tmdata=p4.getYValues()[0];
							if(tmdata==null){showioerror(); return;}
							if(first){
								float[][] tempxvals=p4.getXValues();
								sfreq=1.0/((double)tempxvals[0][1]-(double)tempxvals[0][0]);
								khz=(float)sfreq/1000.0f;
							}
						}
					}
					instream.close();
					if(first){
						int shortlength=tmdata.length;
						int p2length=(int)(Math.log((double)shortlength)/Math.log(2.0));
						if(pad){
							p2length++;
						}
						size=(int)Math.pow(2.0,p2length);
						trajlength=(int)(size/binby);
						if(showtraj){trajectories=new float[nfiles][];}
						acclass=new autocorr(size);
						xvals=bml.getxvals(size/2);
						newsize=xvals.length;
						first=false;
					}
					if(detrend){
						tmdata=(new detrend_linear(tmdata.length,segments)).detrend_array(tmdata);
					}
					double[] kstats=kstatsfunc.kstatisticsshort(tmdata);
					avg[i]=(float)kstats[1];
					var[i]=(float)kstats[2];
					if(showtraj){
						trajectories[i]=bintraj(tmdata,kstats[1]);
					}
					float[] temp=acclass.doautocorr_padded(tmdata,brightcorr)[0];
					float[] tempcorr=bml.dobinmultilog(temp,size/2);
					if(brightcorr){
						for(int j=0;j<tempcorr.length;j++){tempcorr[j]*=khz;}
					}
					correlations[i]=tempcorr;
					IJ.showProgress(i,nfiles);
				} catch(IOException e){
					showioerror();
					return;
				}
			}
		} else {
			ImageWindow iw=WindowManager.getCurrentWindow();
			float[][] trajectories2=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			float[][] tempxvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			sfreq=1.0/((double)tempxvals[0][1]-(double)tempxvals[0][0]);
			khz=(float)sfreq/1000.0f;
			nfiles=trajectories2.length;
			names=new String[nfiles+1];
			names[nfiles]="avg";
			correlations=new Object[nfiles];
			avg=new float[nfiles];
			var=new float[nfiles];
			if(fileflag==3){
				int p2length=(int)(Math.log((double)trajectories2[0].length)/Math.log(2.0));
				if(pad){
					p2length++;
				}
				size=(int)Math.pow(2.0,p2length);
				dl=new detrend_linear(trajectories2[0].length,segments);
				acclass=new autocorr(size);
				xvals=bml.getxvals(size/2);
				newsize=xvals.length;
				showtraj=false;
				for(int i=0;i<nfiles;i++){
					names[i]="trajectory "+(i+1);
					float[] temptraj;
					if(detrend){
						temptraj=dl.detrend_array(trajectories2[i]);
					} else {
						temptraj=trajectories2[i];
					}
					float[][] temp=acclass.doautocorr_padded(temptraj,brightcorr);
					//IJ.log(""+temp[1][0]);
					float[] tempcorr=bml.dobinmultilog(temp[0],size/2);
					if(brightcorr){
						for(int j=0;j<tempcorr.length;j++){
							tempcorr[j]*=khz;
						}
					}
					correlations[i]=tempcorr;
					double[] kstats=kstatsfunc.kstatisticsshort(temptraj,size);
					avg[i]=(float)kstats[1];
					var[i]=(float)kstats[2];
					IJ.showProgress(i,nfiles);
				}		
			} else {
				int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
				newsize=1+(int)jstatistics.getstatistic("Min",npts,null);
				//newsize=tempxvals[0].length+1;
				xvals=new float[newsize];
				for(int i=1;i<newsize;i++){
					xvals[i]=tempxvals[0][i-1]*(float)sfreq;
				}
				for(int i=0;i<nfiles;i++){
					correlations[i]=new float[newsize];
					System.arraycopy(trajectories2[i],0,(float[])correlations[i],1,newsize-1);
					avg[i]=1.0f;
					var[i]=1.0f;
				}
			}
		}
		float[][] corr=new float[nfiles][newsize-1];
		for(int i=0;i<nfiles;i++){
			System.arraycopy((float[])correlations[i],1,corr[i],0,newsize-1);
		}
		float[] newxvals=new float[newsize-1];
		for(int i=0;i<newsize-1;i++){
			newxvals[i]=xvals[i+1]/(float)sfreq;
		}
		if(simple){
			double[] dxvals=new double[newsize-1];
			for(int i=0;i<newsize-1;i++) dxvals[i]=(double)newxvals[i];
			//here we do a simple grid search fitting
			GenericDialog gd10=new GenericDialog("Options");
			gd10.addNumericField("Min taud(ms)",1.0f,5,15,null);
			gd10.addNumericField("Max taud(s)",5.0f,5,15,null);
			gd10.addNumericField("Multiplier",1.05,5,15,null);
			gd10.showDialog(); if(gd10.wasCanceled()){return;}
			float mintd=0.001f*(float)gd10.getNextNumber();
			float maxtd=(float)gd10.getNextNumber();
			float mult=(float)gd10.getNextNumber();
			fit_corr fc=new fit_corr(mintd,maxtd,mult,5.0); fc.psftype=psfflag;
			TextWindow tw=jutils.selectTable("AC Results");
			if(tw==null){
				String labels="filename\tbase\tg0\ttd\tc2\tavg(kHz)";
				if(brightcorr) labels="filename\tbase\tB(kHz)\ttd\tc2\tavg(kHz)";
				tw=new TextWindow("AC Results",labels,"",400,400);
			}
			float[][][] corrfit=new float[nfiles][2][];
			for(int i=0;i<nfiles;i++){
				double[] params=fc.fitac(corr[i],newxvals,false,true);
				corrfit[i][0]=corr[i];
				corrfit[i][1]=fc.corfunc_arrayf(params,dxvals);
				tw.append(names[i]+"\t"+table_tools.print_double_array(params)+"\t"+khz*avg[i]);
			}
			if(nfiles==1){
				float[][] tempxvals5=new float[2][]; for(int i=0;i<2;i++) tempxvals5[i]=newxvals;
				PlotWindow4 pwt=null;
				if(brightcorr) pwt=new PlotWindow4("Avg","tau(s)","B(kHz)",tempxvals5,corrfit[0],null);
				else pwt=new PlotWindow4("Avg","tau(s)","G(tau)",tempxvals5,corrfit[0],null);
				pwt.setLogAxes(true,false);
				pwt.draw();
			} else {
				PlotStack4 ps=null;
				if(brightcorr) ps=new PlotStack4("Correlations","tau(s)","G(tau)",newxvals,corrfit);
				else ps=new PlotStack4("Correlations","tau(s)","G(tau)",newxvals,corrfit);
				ps.draw();
				ps.setAllLogAxes(true,false);
			}
		} else {
			final AutoCorrFitWindow cw = new AutoCorrFitWindow();
			cw.init(names,corr,newxvals,trajectories,avg,var,khz,psfflag,brightcorr);
			AutoCorrFitWindow.launch_frame(cw);
		}
	}

	float[] bintraj(float[] data,double avgint){
		int newlength=(int)(data.length/binby);
		if(newlength>trajlength){newlength=trajlength;}
		float[] newtraj=new float[trajlength];
		for(int i=0;i<newlength;i++){
			for(int j=0;j<binby;j++){
				newtraj[i]+=data[j+i*binby];
			}
			newtraj[i]*=trajkhz;
		}
		for(int i=newlength;i<trajlength;i++){
			newtraj[i]=(float)avgint*(float)binby*trajkhz;
		}
		return newtraj;
	}			

	private void showioerror(){
		IJ.showMessage("Error in file io");
	}

}
