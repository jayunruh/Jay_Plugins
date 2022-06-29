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
import jalgs.jfft.*;
import jguis.*;

public class traj_fft_2corr_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here are all incarnations of dual color fcs
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		int length=yvals[0].length;
		float xinc=xvals[0][1]-xvals[0][0];

		int p2length=(int)(Math.log((double)length)/Math.log(2.0));
		p2length++;
		String[] p2lengths=new String[p2length-1];
		for(int i=0;i<p2length-1;i++){
			p2lengths[i]=""+(int)(Math.pow(2.0,p2length-i));
		}
		GenericDialog gd2 = new GenericDialog("Options");
		gd2.addChoice("Analysis_length",p2lengths,p2lengths[0]);
		gd2.addCheckbox("Do_Autocorr1?",true);
		gd2.addCheckbox("Do_Crosscorr?",true);
		gd2.addCheckbox("Do_Autocorr2?",true);
		gd2.addCheckbox("Do_Tricorr1?",false);
		gd2.addCheckbox("Do_Tricorr2?",false);
		gd2.addCheckbox("Do_Tricorr21?",false);
		gd2.addCheckbox("Do_Tricorr12?",false);
		boolean binlog=true;
		gd2.addCheckbox("log_bin?",binlog);
		boolean brightcorr=false;
		gd2.addCheckbox("Brightcorr?",brightcorr);
		gd2.showDialog();
		if(gd2.wasCanceled()){return;}
		p2length = p2length-gd2.getNextChoiceIndex();
		int ncurves=0;
		boolean[] selections=new boolean[7];
		for(int i=0;i<7;i++){
			selections[i]=gd2.getNextBoolean();
			if(selections[i]){ncurves++;}
		}
		binlog=gd2.getNextBoolean();
		brightcorr=gd2.getNextBoolean();
		int size = (int)Math.pow(2.0,p2length);
		int segments=(int)(length/size);
		if(segments==0){segments=1;}		

		if(ncurves>0){
			float[][] cc = new float[ncurves][size/2];
			crosscorr ccclass=new crosscorr(size);
			autocorr acclass=new autocorr(size);
			triautocorr tcclass=new triautocorr(size);
			corr21 tc21class=new corr21(size);
			for(int k=0;k<segments;k++){
				float[] tempfloat,tempfloat2;
				if(size<=length){
					tempfloat=new float[size];
					tempfloat2=new float[size];
					System.arraycopy(yvals[0],k*size,tempfloat,0,size);
					if(yvals.length>1) System.arraycopy(yvals[1],k*size,tempfloat2,0,size);
				} else {
					tempfloat=new float[length];
					tempfloat2=new float[length];
					System.arraycopy(yvals[0],k*size,tempfloat,0,length);
					if(yvals.length>1) System.arraycopy(yvals[1],k*size,tempfloat2,0,length);
				}
				float[] tempfloat3;
				int counter=0;
				if(selections[0]){
					tempfloat3=acclass.doautocorr_padded(tempfloat,brightcorr)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[1]){
					tempfloat3=ccclass.docrosscorr_padded(tempfloat,tempfloat2,brightcorr)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[2]){
					tempfloat3=acclass.doautocorr_padded(tempfloat2,brightcorr)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[3]){
					tempfloat3=tcclass.doautocorr_padded(tempfloat,brightcorr)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[4]){
					tempfloat3=tcclass.doautocorr_padded(tempfloat2,brightcorr)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[5]){
					tempfloat3=tc21class.docrosscorr_padded(tempfloat,tempfloat2,brightcorr?1:0)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
				if(selections[6]){
					tempfloat3=tc21class.docrosscorr_padded(tempfloat2,tempfloat,brightcorr?1:0)[0];
					for(int j=0;j<(size/2);j++){cc[counter][j]+=tempfloat3[j]/(float)segments;}
					counter++;
				}
			}
			if(binlog){
				binmultilog bml=new binmultilog();
				float[] xvals5=bml.getxvals(size/2);
				int newsize=xvals5.length;
				float[][] newresult=new float[ncurves][newsize];
				for(int i=0;i<ncurves;i++){
					float[] temp=new float[size/2];
					System.arraycopy(cc[i],0,temp,0,size/2);
					temp=bml.dobinmultilog(temp,size/2);
					System.arraycopy(temp,0,newresult[i],0,newsize);
				}
				for(int j=0;j<newsize;j++){
					xvals5[j]*=xinc;
				}
				PlotWindow4 pw5=new PlotWindow4("Correlation","tau","G(tau)",xvals5,newresult[0]);
				pw5.draw();
				for(int i=1;i<ncurves;i++){
					pw5.addPoints(xvals5,newresult[i],true);
				}
				pw5.setLogAxes(true,false);	
			}
			else {
				float[] xvals5=new float[size/2];
				for(int i=0;i<(size/2);i++){xvals5[i]=(float)i*xinc;}
				PlotWindow4 pw5=new PlotWindow4("Correlation","tau","G(tau)",xvals5,cc[0]);
				pw5.draw();
				for(int i=1;i<ncurves;i++){
					pw5.addPoints(xvals5,cc[i],true);
				}
				pw5.setLogAxes(true,false);
			}
		}
	}
}
