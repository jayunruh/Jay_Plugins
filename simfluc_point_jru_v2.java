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
import jalgs.*;
import jalgs.jsim.*;
import jguis.*;
import java.io.*;

public class simfluc_point_jru_v2 implements PlugIn, simflucinterface  {
	boolean output;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("2D",true);
		gd.addNumericField("Time_step(us)",50.0,5,15,null);
		gd.addNumericField("#_of_steps",1048576,0,15,null);
		gd.addNumericField("w0(um)",0.17f,5,15,null);
		gd.addNumericField("z0/w0",5.0f,5,15,null);
		gd.addNumericField("#_of_species1",100,0);
		gd.addNumericField("D1(um^2/sec)",1.0,5,15,null);
		gd.addNumericField("B1_ch1(cpsm)",100000.0,5,15,null);
		gd.addNumericField("B1_ch2(cpsm)",0.0,5,15,null);
		gd.addNumericField("Bsteps1_ch1",1,0);
		gd.addNumericField("Bsteps1_ch2",1,0);
		gd.addNumericField("#_of_species2",0,0);
		gd.addNumericField("D2(um^2/sec)",0.01,5,15,null);
		gd.addNumericField("B2_ch1(cpsm)",100000.0,5,15,null);
		gd.addNumericField("B2_ch2(cpsm)",0.0,5,15,null);
		gd.addNumericField("Bsteps2_ch1",1,0);
		gd.addNumericField("Bsteps2_ch2",1,0);
		gd.addNumericField("Bleach_rate_ch1(1/s_or_1/photons)",0.0,5,15,null);
		gd.addNumericField("Bleach_rate_ch2(1/s_or_1/photons)",0.0,5,15,null);
		gd.addCheckbox("Photon_Budget_Bleaching",false);
		gd.addNumericField("Triplet_fraction",0.2,5,15,null);
		gd.addNumericField("Triplet_tau(us)",100.0,5,15,null);
		gd.addNumericField("Back1(cps)",0.0,5,15,null);
		gd.addNumericField("Back2(cps)",0.0,5,15,null);
		gd.addCheckbox("Save_Output",false);
		String userdir=System.getProperty("user.home");
		userdir=userdir.replace('\\','/');
		gd.addStringField("Save_path(use/)",userdir+"/sim1.pw");
		gd.addCheckbox("Output_progress",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean conf=gd.getNextBoolean();
		float exptime=(float)gd.getNextNumber();
		int frames=(int)gd.getNextNumber();
		float w0=(float)gd.getNextNumber();
		float zratio=(float)gd.getNextNumber();
		int N1=(int)gd.getNextNumber();
		float D1=(float)gd.getNextNumber();
		float B11=(float)gd.getNextNumber();
		float B12=(float)gd.getNextNumber();
		int bsteps11=(int)gd.getNextNumber();
		int bsteps12=(int)gd.getNextNumber();
		int N2=(int)gd.getNextNumber();
		float D2=(float)gd.getNextNumber();
		float B21=(float)gd.getNextNumber();
		float B22=(float)gd.getNextNumber();
		int bsteps21=(int)gd.getNextNumber();
		int bsteps22=(int)gd.getNextNumber();
		float brate=(float)gd.getNextNumber();
		float brate2=(float)gd.getNextNumber();
		boolean bleachbyphoton=gd.getNextBoolean();
		float Tf=(float)gd.getNextNumber();
		float Ttau=(float)gd.getNextNumber();
		float back1=(float)gd.getNextNumber();
		float back2=(float)gd.getNextNumber();
		boolean save=gd.getNextBoolean();
		String savepath=gd.getNextString();
		output=gd.getNextBoolean();
		if(save){
			savepath=savepath.replace('/',File.separatorChar);
		}
		
		//need to set up four species, every other one dark (triplet)
		float[] D={D1,D1,D2,D2};
		float[][] bright=new float[4][2];
		bright[0][0]=B11; bright[0][1]=B12; bright[2][0]=B21; bright[2][1]=B22;
		int[] num=new int[4];
		num[1]=(int)(Tf*(float)N1); num[0]=N1-num[1];
		num[3]=(int)(Tf*(float)N2); num[2]=N2-num[3];
		float[][] trate=new float[4][4]; //indices are [from][to];
		float ktot=1000000.0f/Ttau;
		trate[0][1]=Tf*ktot; trate[1][0]=ktot-trate[0][1];
		trate[2][3]=Tf*ktot; trate[3][2]=ktot-trate[2][3];
		float[][] kbleach={{brate,brate2},{brate,brate2},{brate,brate2},{brate,brate2}};
		float[] back={back1,back2};
		int[][] bsteps={{bsteps11,bsteps12},{bsteps11,bsteps12},{bsteps21,bsteps22},{bsteps21,bsteps22}};
		
		simfluc_scanning ss=new simfluc_scanning(this);
		ss.bleachbyphoton=bleachbyphoton;
		int confineindex=conf?1:0;
		float[][] data=ss.do_simfluc_point(3.2f,64,w0,zratio*w0,frames,exptime,back,0.0f,confineindex,false,false,false,
			0,0,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,num,bright,D,trate,kbleach,bsteps);
		float[] xvals=new float[frames];
		float xinc=exptime*(float)(1.0e-6);
		for(int i=0;i<frames;i++){
			xvals[i]=xinc*(float)i;
		}
		PlotWindow4 pw=new PlotWindow4("Point Simulation","time","Intensity",xvals,data[0]);
		pw.draw();
		if(B12>0.0f || B22>0.0f){
			pw.addPoints(xvals,data[1],true);
		}
		if(save){
			pw.saveAsObject(savepath);
			pw.close();
		}
	}

	public void showprogress(int i,int end){
		if(output) IJ.showProgress(i,end);
	}

	public void showmessage(String message){
		if(output) IJ.log(message);
	}

}
