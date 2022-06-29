/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jsim.*;
import jguis.*;

public class sim_kt_tracking_jru_v1 implements PlugIn {
	double[] kouttable,kintable;

	public void run(String arg) {
		//this plugin simulates microtubule end tracking according to joglekar and hunt, 2002, biophys j.
		//double kappa=1800; //units s^-1
		double kappa=7200;
		double s=33.0/340.0; //unitless
		double r=0.96; //unitless
		double beta=340.0; //units s^-1
		//double beta=85.0; //version for stu2 which has slower depolymerization
		double tension=20.0; //units are piconewton, 20 works for Boris' equation, 4 works for joglekar version

		int M=8; //total binding sites
		int N=M; //initial binding position
		double tau=0.00001; //time step in seconds/frame, used 0.00001 for Boris' model
		double tottime=1000.0; //sim time in seconds
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("kappa (1/s)",kappa,5,15,null);
		gd.addNumericField("s",s,5,15,null);
		gd.addNumericField("r",r,5,15,null);
		gd.addNumericField("beta (1/s)",beta,5,15,null);
		gd.addNumericField("tension (pN)",tension,5,15,null);
		gd.addNumericField("binding_sites (M)",M,5,15,null);
		gd.addNumericField("sim_time (s)",tottime,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		kappa=gd.getNextNumber();
		s=gd.getNextNumber();
		r=gd.getNextNumber();
		beta=gd.getNextNumber();
		tension=gd.getNextNumber();
		M=(int)gd.getNextNumber();
		tottime=gd.getNextNumber();

		double f=Math.exp(-tension*0.615/(2.0*4.114)); //kT is 4.114 pN*nm at rt
		IJ.log("f = "+f);
		//double f=0.3f; //unitless--load effect, smaller value means we are resisting pulling
		int frames=(int)(tottime/tau);
		int skip=100;
		float[] traj=new float[frames/skip];
		float[] xvals=new float[frames/skip];
		for(int i=0;i<frames/skip;i++) xvals[i]=(float)(tau*(double)i*(double)skip);
		//double kappa2=kappa*tau;
		//double beta2=beta*tau;
		buildtables(kappa,s,r,f,beta,M,tau);
		rngs random=new rngs();
		double len=run_sim(N,M,frames,true,random,skip,traj,xvals,tau);
		new PlotWindow4("kt_tracking_sim","time (s)","rel. position",xvals,traj).draw();
		int nsims=1000;
		float[] lens=new float[nsims];
		for(int i=0;i<nsims;i++){
			lens[i]=(float)run_sim(N,M,frames,false,random,skip,traj,xvals,tau);
			IJ.showStatus("run "+i+" of "+nsims);
		}
		new PlotWindow4("kt_loss_time","run","loss time (s)",lens).draw();
	}

	public double run_sim(int N,int M,int frames,boolean plot,rngs random,int skip,float[] traj,float[] xvals,double tau){
		if(plot){
			traj[0]=N;
			int counter=1;
			for(int i=1;i<frames;i++){
				double kin1=kintable[N];
				double kout1=kouttable[N];
				double random1=random.unidev(1.0,0.0);
				if(random1<kin1){
					if(random1>=0.25) IJ.log("time step exceeded");
					N++;
					if(N>M) N=M;
				} else {
					if(random1>(1.0-kout1)){
						if(random1<0.75) IJ.log("time step exceeded");
						N--;
						if(N<=0){
							N=0;
							double ltime=tau*(double)i;
							IJ.log("lost microtubule at t = \t"+ltime);
							//break;
							return ltime;
						}
					}
				}
				if((i%skip)==0 && counter<traj.length){
					traj[counter]=(float)N;
					counter++;
					IJ.showProgress(i,frames);
				}
				if(IJ.escapePressed()) return -2.0f;
			}
			return xvals[xvals.length-1];
		} else {
			int counter=1;
			for(int i=1;i<frames;i++){
				double kin1=kintable[N];
				double kout1=kouttable[N];
				double random1=random.unidev(1.0,0.0);
				if(random1<kin1){
					if(random1>=0.25) IJ.log("time step exceeded");
					N++;
					if(N>M) N=M;
				} else {
					if(random1>(1.0-kout1)){
						if(random1<0.75) IJ.log("time step exceeded");
						N--;
						if(N<=0){
							N=0;
							double ltime=tau*(double)i;
							IJ.log("lost microtubule at t =\t "+ltime);
							//break;
							return xvals[counter];
						}
					}
				}
				if((i%skip)==0 && counter<traj.length){
					//traj[counter]=(float)N;
					counter++;
					IJ.showProgress(i,frames);
				}
				if(IJ.escapePressed()) return -2.0f;
			}
			return xvals[xvals.length-1];
		}
	}

	public void buildtables(double kappa,double s,double r,double f,double beta,int M,double tau){
		kouttable=new double[M+1];
		kintable=new double[M+1];
		for(int i=0;i<=M;i++){
			kouttable[i]=kout(kappa,s,r,f,beta,M,i);
			kintable[i]=kin(kappa,s,r,f,beta,M,i);
			kouttable[i]*=tau;
			kintable[i]*=tau;
			IJ.log("N = \t"+i+"\t kin = \t"+kintable[i]+"\t kout = \t"+kouttable[i]);
		}
	}

	public double kout(double kappa,double s,double r,double f,double beta,int M,int N){
		return kappa*s*intpow(r,M-N)/f+beta*s;
	}

	public double kin(double kappa,double s,double r,double f,double beta,int M,int N){
		return kappa*intpow(r,M-N)*f/s; //this is Boris's version
		//return kappa*intpow(r,M-N)*f; //this is the joglekar version
	}

	public double intpow(double val,int pow){
		if(pow<=0) return 1.0;
		double retval=val;
		for(int i=1;i<pow;i++){
			retval=retval*val;
		}
		return retval;
	}

}
