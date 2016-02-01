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
import javax.swing.*;
import java.io.*;
import jalgs.*;
import jalgs.jfft.*;
import jalgs.jfit.*;
import jguis.jutils;
import jguis.CrossCorrFitWindow;
import jguis.table_tools;
import jguis.PlotStack4;
import jguis.PlotWindow4;
import java.awt.event.*;
import ij.io.*;
import java.text.*;
import ij.text.*;

public class analysis_cross_corr_v2 implements PlugIn {
	//this plugin is a gui for fitting correlation curves singly or globally for 
	//3D Gaussian or Gaussian Lorentzian squared point spread functions
	//(solution confocal, solution two photon)
	//copyright 2009 Jay Unruh, Stowers Institute for Medical Research
	int binby,trajlength;
	float trajkhz,khz;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=20000.0;
		gd.addNumericField("Sampling Frequency?",sfreq,1,10,null);
		String[] psfchoice={"3D Gaussian","2D Gaussian","2Dxz_Gaussian"};
		gd.addChoice("PSF Type?",psfchoice,psfchoice[0]);
		String[] filetypechoice={"Confocor 3 raw","Short binary trajectory","PlotWindow trajectory","Confocor 3 ALEX"};
		gd.addChoice("File Type?",filetypechoice,filetypechoice[0]);
		boolean ch2green=true;
		gd.addCheckbox("Ch2 is green?",ch2green);
		boolean showtraj=true;
		gd.addCheckbox("Show Trajectories?",showtraj);
		binby=20;
		gd.addNumericField("Traj Bin?",binby,0,10,null);
		boolean detrend=false;
		gd.addCheckbox("Segmented Linear Detrending?",detrend);
		int segments=2;
		gd.addNumericField("Detrending Segments?",segments,0);
		boolean brightcorr=false;
		gd.addCheckbox("Brightcorr?",brightcorr);
		gd.addCheckbox("Pad Data?",false);
		gd.addCheckbox("Simple Analysis?",false);
		gd.addNumericField("ALEX_illumination_delay(us)",2.5,5,15,null);
		int pmfreq=20000000;
		gd.addNumericField("Photon Mode Freq",pmfreq,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=gd.getNextNumber();
		int psfflag=gd.getNextChoiceIndex();
		int fileflag=gd.getNextChoiceIndex();
		ch2green=gd.getNextBoolean();
		showtraj=gd.getNextBoolean();
		binby=(int)gd.getNextNumber();
		detrend=gd.getNextBoolean();
		segments=(int)gd.getNextNumber();
		brightcorr=gd.getNextBoolean();
		boolean pad=gd.getNextBoolean();
		boolean simple=gd.getNextBoolean();
		double ill_delay=gd.getNextNumber();
		ill_delay*=1.0e-6;
		pmfreq=(int)gd.getNextNumber();
		trajkhz=((float)sfreq/(float)binby)/1000.0f;
		khz=(float)sfreq/1000.0f;
		int nfiles=0;
		Object[] correlations=null;
		float[][][] trajectories=null;
		int xmax=0;
		int ymax=0;
		String[] names=null;
		boolean first=true;
		autocorr acclass=null;
		crosscorr ccclass=null;
		detrend_linear dl=null;
		int size=0;
		int newsize=0;
		float[] xvals=null;
		float[][] avg=null;
		float[][] var=null;
		trajlength=0;
		binmultilog bml=new binmultilog();
		kstats kstatsfunc=new kstats();
		if(fileflag!=2){
			jdataio ioclass=new jdataio();
			File[] filearray=ioclass.openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
			if(filearray.length==0){return;}
			String dir=filearray[0].getAbsolutePath();
			int sepindex=dir.lastIndexOf(File.separator);
			String newdir=dir.substring(0,sepindex+1);
			OpenDialog.setDefaultDirectory(newdir);
			nfiles=filearray.length/2;
			correlations=new Object[nfiles];
			avg=new float[3][nfiles];
			var=new float[3][nfiles];
			names=organize_c3_files(filearray);
			for(int i=0;i<nfiles;i++){
				try{
					int length1=(int)(((double)filearray[2*i].length()-128.0)/4.0);
					int length2=(int)(((double)filearray[2*i+1].length()-128.0)/4.0);
					int length3=(int)(((double)filearray[2*i].length())/2.0);
					int length4=(int)(((double)filearray[2*i+1].length())/2.0);
					InputStream instream=new BufferedInputStream(new FileInputStream(filearray[2*i]));
					InputStream instream2=new BufferedInputStream(new FileInputStream(filearray[2*i+1]));
					float[] tmdata,tmdata2;
					if(fileflag==0 || fileflag==3){
						int[] pmdata,pmdata2;
						if(!ioclass.skipstreambytes(instream,128)){showioerror(); instream.close(); return;}
						if(!ioclass.skipstreambytes(instream2,128)){showioerror(); instream2.close(); return;}
						if(ch2green){
							pmdata=new int[length2];
							pmdata2=new int[length1];
							if(!ioclass.readintelintfile(instream,length1,pmdata2)){showioerror(); instream.close(); return;}
							if(!ioclass.readintelintfile(instream2,length2,pmdata)){showioerror(); instream2.close(); return;}
						} else {
							pmdata=new int[length1];
							pmdata2=new int[length2];
							if(!ioclass.readintelintfile(instream,length1,pmdata)){showioerror(); instream.close(); return;}
							if(!ioclass.readintelintfile(instream2,length2,pmdata2)){showioerror(); instream2.close(); return;}
						}
						if(fileflag==3){
							double swfreq=(double)pmfreq/1000.0111;
							double divider=20000.0/sfreq;
							swfreq/=divider;
							/*if(sfreq==20000.0){
							} else {
								if(sfreq==10000.0){
									swfreq/=2.0;
								} else {
									if(sfreq==5000.0){
										swfreq/=4.0;
									} else {
										if(sfreq==2000.0){
											swfreq/=10.0;
										} else {
											IJ.showMessage("Incorrect switching frequency");
											return;
										}
									}
								}
							}*/
							int offset=0;
							offset=(new pmodeconvert()).calc_wlswitch_offset(pmdata,swfreq,pmfreq);
							float[][] tmdatatemp=(new pmodeconvert()).pm2tm_alex(pmdata,pmdata2,swfreq,pmfreq,offset,1,ill_delay);
							tmdata=tmdatatemp[0];
							tmdata2=tmdatatemp[1];
						} else {
							tmdata=(new pmodeconvert()).pm2tm(pmdata,sfreq,pmfreq);
							tmdata2=(new pmodeconvert()).pm2tm(pmdata2,sfreq,pmfreq);
						}
					} else {
						if(ch2green){
							tmdata=new float[length4]; tmdata2=new float[length3];
							if(!ioclass.readintelshortfile(instream,length3,tmdata2)){showioerror(); instream.close(); return;}
							if(!ioclass.readintelshortfile(instream2,length4,tmdata)){showioerror(); instream2.close(); return;}
						} else {
							tmdata=new float[length3]; tmdata2=new float[length4];
							if(!ioclass.readintelshortfile(instream,length3,tmdata)){showioerror(); instream.close(); return;}
							if(!ioclass.readintelshortfile(instream2,length4,tmdata2)){showioerror(); instream2.close(); return;}
						}
					}
					if(first){
						int shortlength=tmdata.length;
						if(tmdata2.length<shortlength){shortlength=tmdata2.length;}
						int p2length=(int)(Math.log((double)shortlength)/Math.log(2.0));
						if(pad){
							p2length++;
						}
						size=(int)Math.pow(2.0,p2length);
						trajlength=(int)(size/binby);
						if(showtraj){trajectories=new float[nfiles][2][trajlength];}
						acclass=new autocorr(size);
						ccclass=new crosscorr(size);
						xvals=bml.getxvals(size/2);
						newsize=xvals.length;
						first=false;
					}
					if(detrend){
						tmdata=(new detrend_linear(tmdata.length,segments)).detrend_array(tmdata);
						tmdata2=(new detrend_linear(tmdata2.length,segments)).detrend_array(tmdata2);
					}
					double[][] kstats=kstatsfunc.kstatisticsshort(tmdata,tmdata2);
					avg[0][i]=(float)kstats[1][0]; avg[1][i]=(float)kstats[0][1]; avg[2][i]=(float)Math.sqrt(kstats[1][0]*kstats[0][1]);
					var[0][i]=(float)kstats[2][0]; var[1][i]=(float)kstats[0][2]; var[2][i]=(float)kstats[1][1];
					if(showtraj){
						trajectories[i][0]=bintraj(tmdata,kstats[1][0]);
						trajectories[i][1]=bintraj(tmdata2,kstats[0][1]);
					}
					float[][] tempcorr=new float[3][];
					float[] temp=acclass.doautocorr_padded(tmdata,brightcorr)[0];
					tempcorr[0]=bml.dobinmultilog(temp,size/2);
					temp=acclass.doautocorr_padded(tmdata2,brightcorr)[0];
					tempcorr[1]=bml.dobinmultilog(temp,size/2);
					temp=ccclass.docrosscorr_padded(tmdata,tmdata2,brightcorr)[0];
					tempcorr[2]=bml.dobinmultilog(temp,size/2);
					if(brightcorr){
						for(int j=0;j<3;j++){
							for(int k=0;k<tempcorr[j].length;k++){
								tempcorr[j][k]*=khz;
							}
						}
					}
					correlations[i]=tempcorr;
					instream.close();
					instream2.close();
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
			int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
			sfreq=1.0/((double)tempxvals[0][1]-(double)tempxvals[0][0]);
			trajkhz=((float)sfreq/(float)binby)/1000.0f;
			khz=(float)sfreq/1000.0f;
			nfiles=trajectories2.length/2;
			names=new String[nfiles+1];
			names[nfiles]="avg";
			correlations=new Object[nfiles];
			avg=new float[3][nfiles];
			var=new float[3][nfiles];
			int shortest=npts[0];
			for(int i=0;i<npts.length;i++){if(shortest<npts[i]){shortest=npts[i];}}
			int p2length=(int)(Math.log((double)shortest)/Math.log(2.0));
			if(pad){
				p2length++;
			}
			size=(int)Math.pow(2.0,p2length);
			acclass=new autocorr(size);
			ccclass=new crosscorr(size);
			xvals=bml.getxvals(size/2);
			newsize=xvals.length;
			first=false;
			showtraj=false;
			for(int i=0;i<nfiles;i++){
				names[i]="trajectory "+(i+1);
				float[] temptraj1=null;
				float[] temptraj2=null;
				if(detrend){
					if(ch2green){
						temptraj2=(new detrend_linear(trajectories2[2*i].length,segments)).detrend_array(trajectories2[2*i]);
						temptraj1=(new detrend_linear(trajectories2[2*i+1].length,segments)).detrend_array(trajectories2[2*i+1]);
					} else {
						temptraj1=(new detrend_linear(trajectories2[2*i].length,segments)).detrend_array(trajectories2[2*i]);
						temptraj2=(new detrend_linear(trajectories2[2*i+1].length,segments)).detrend_array(trajectories2[2*i+1]);
					}
				} else {
					if(ch2green){
						temptraj2=trajectories2[2*i];
						temptraj1=trajectories2[2*i+1];
					} else {
						temptraj1=trajectories2[2*i];
						temptraj2=trajectories2[2*i+1];
					}
				}
				float[][] tempcorr=new float[3][newsize];
				float[] temp=acclass.doautocorr_padded(temptraj1,brightcorr)[0];
				tempcorr[0]=bml.dobinmultilog(temp,size/2);
				temp=acclass.doautocorr_padded(temptraj2,brightcorr)[0];
				tempcorr[1]=bml.dobinmultilog(temp,size/2);
				temp=ccclass.docrosscorr_padded(temptraj1,temptraj2,brightcorr)[0];
				tempcorr[2]=bml.dobinmultilog(temp,size/2);
				if(brightcorr){
					for(int j=0;j<3;j++){
						for(int k=0;k<tempcorr[j].length;k++){
							tempcorr[j][k]*=khz;
						}
					}
				}
				correlations[i]=tempcorr;
				double[][] kstats=kstatsfunc.kstatisticsshort(temptraj1,temptraj2,size);
				avg[0][i]=(float)kstats[1][0]; avg[1][i]=(float)kstats[0][1]; avg[2][i]=(float)Math.sqrt(kstats[1][0]*kstats[0][1]);
				var[0][i]=(float)kstats[2][0]; var[1][i]=(float)kstats[0][2]; var[2][i]=(float)kstats[1][1];
				IJ.showProgress(i,nfiles);
			}
		}
		float[][][] corr=new float[nfiles][3][newsize-1];
		for(int i=0;i<nfiles;i++){
			for(int j=0;j<3;j++){
				System.arraycopy(((float[][])correlations[i])[j],1,corr[i][j],0,newsize-1);
			}
		}
		float[] newxvals=new float[newsize-1];
		for(int i=0;i<newsize-1;i++){
			newxvals[i]=xvals[i+1]/(float)sfreq;
		}
		if(simple){
			double[] dxvals=new double[newsize-1];
			for(int i=0;i<newsize-1;i++) dxvals[i]=(double)newxvals[i];
			//here we do a simple grid search fitting
			//the cross corr taud can be fixed to the longest taud, the avg taud, or fit in the analysis
			GenericDialog gd10=new GenericDialog("Options");
			gd10.addNumericField("Min_green taud(ms)",1.0f,5,15,null);
			gd10.addNumericField("Max_green taud(s)",5.0f,5,15,null);
			gd10.addNumericField("Min_red taud(ms)",1.0f,5,15,null);
			gd10.addNumericField("Max_red taud(s)",5.0f,5,15,null);
			gd10.addNumericField("Min_cc taud(ms)",1.0f,5,15,null);
			gd10.addNumericField("Max_cc taud(s)",5.0f,5,15,null);
			String[] ccoptions={"Avg_green_red","Max_green_red","fit"};
			gd10.addChoice("CC_taud_fitting",ccoptions,ccoptions[0]);
			gd10.showDialog(); if(gd10.wasCanceled()){return;}
			float mingtd=0.001f*(float)gd10.getNextNumber();
			float maxgtd=(float)gd10.getNextNumber();
			float minrtd=0.001f*(float)gd10.getNextNumber();
			float maxrtd=(float)gd10.getNextNumber();
			float mincctd=0.001f*(float)gd10.getNextNumber();
			float maxcctd=(float)gd10.getNextNumber();
			int ccoptindex=gd10.getNextChoiceIndex();
			fit_corr fcg=new fit_corr(mingtd,maxgtd,1.05,5.0); fcg.psftype=psfflag;
			fit_corr fcr=new fit_corr(minrtd,maxrtd,1.05,5.0); fcr.psftype=psfflag;
			fit_corr fccc=new fit_corr(mincctd,maxcctd,1.05,5.0); fccc.psftype=psfflag;
			TextWindow tw=jutils.selectTable("Results");
			if(tw==null){
				String labels="filename\tbaseg\tg0g\ttdg\tc2g\tbaser\tg0r\ttdr\tc2r\tbasecc\tg0cc\ttdcc\tc2cc\tIg\tIr";
				tw=new TextWindow("Results",labels,"",400,400);
			}
			float[][][] corrfit=new float[nfiles][6][];
			for(int i=0;i<nfiles;i++){
				double[] gparams=fcg.fitac(corr[i][0],newxvals,false,true);
				double[] rparams=fcr.fitac(corr[i][1],newxvals,false,true);
				double[] ccparams=new double[gparams.length];
				if(ccoptindex==0){
					double td=0.5*(gparams[2]+rparams[2]);
					double[] temp=fccc.fit_linear_ac(td,corr[i][2],newxvals,false,true);
					ccparams[0]=temp[0]; ccparams[1]=temp[1]; ccparams[2]=td;
					ccparams[3]=fccc.c2(ccparams,corr[i][2],newxvals,false);
				} else {
					if(ccoptindex==1){
						double td=Math.max(gparams[2],rparams[2]);
						double[] temp=fccc.fit_linear_ac(td,corr[i][2],newxvals,false,true);
						ccparams[0]=temp[0]; ccparams[1]=temp[1]; ccparams[2]=td;
						ccparams[3]=fccc.c2(ccparams,corr[i][2],newxvals,false);
					} else {
						ccparams=fccc.fitac(corr[i][2],newxvals,false,true);
					}
				}
				corrfit[i][0]=corr[i][0]; corrfit[i][1]=corr[i][1]; corrfit[i][2]=corr[i][2];
				corrfit[i][3]=fcg.corfunc_arrayf(gparams,dxvals);
				corrfit[i][4]=fcr.corfunc_arrayf(rparams,dxvals);
				corrfit[i][5]=fccc.corfunc_arrayf(ccparams,dxvals);
				tw.append(names[i]+"\t"+table_tools.print_double_array(gparams)+"\t"+table_tools.print_double_array(rparams)+"\t"+table_tools.print_double_array(ccparams)+"\t"+khz*avg[0][i]+"\t"+khz*avg[1][i]);
			}
			if(nfiles==1){
				float[][] tempxvals5=new float[6][]; for(int i=0;i<6;i++) tempxvals5[i]=newxvals;
				PlotWindow4 pwt=new PlotWindow4("Avg","tau(s)","G(tau)",tempxvals5,corrfit[0],null);
				pwt.setLogAxes(true,false);
				pwt.draw();
			} else {
				PlotStack4 ps=new PlotStack4("Correlations","tau(s)","G(tau)",newxvals,corrfit);
				ps.draw();
				ps.setAllLogAxes(true,false);
			}
		} else {
			//here we do more advanced averaging or global analysis
			final CrossCorrFitWindow cw = new CrossCorrFitWindow();
			cw.init(names,corr,newxvals,trajectories,avg,var,khz,psfflag,brightcorr);
			CrossCorrFitWindow.launch_frame(cw);
		}
	}

	private String[] organize_c3_files(File[] filearray){
		int length=filearray.length;
		int[] assign=new int[length];
		for(int i=0;i<length;i++){
			assign[i]=-1;
		}
		int counter=0;
		File[] temp=new File[length];
		String[] outnames=new String[length/2+1];
		outnames[length/2]="avg";
		for(int i=0;i<length;i++){
			if(assign[i]<0){
				String bait=filearray[i].getName();
				String baittrunc=bait.substring(0,bait.length()-5);
				for(int j=i+1;j<length;j++){
					if(assign[j]<0){
						String target=filearray[j].getName();
						if(target.substring(0,target.length()-5).equals(baittrunc)){
							assign[j]=counter;
							assign[i]=counter;
							if(bait.charAt(bait.length()-5)>target.charAt(target.length()-5)){
								temp[2*counter]=filearray[i];
								temp[2*counter+1]=filearray[j];
								outnames[counter]=target.substring(0,target.length()-8);
							} else {
								temp[2*counter]=filearray[j];
								temp[2*counter+1]=filearray[i];
								outnames[counter]=bait.substring(0,bait.length()-8);
							}
							counter++;
							break;
						}
					}
				}
			}
		}
		filearray=temp;
		return outnames;
	}
			

	private void showioerror(){
		IJ.showMessage("Error in file io");
	}

	float[] bintraj(float[] data,double avgint){
		int newlength=(int)(data.length/binby);
		float[] newtraj=new float[trajlength];
		if(newlength>trajlength){newlength=trajlength;}
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

	

}
