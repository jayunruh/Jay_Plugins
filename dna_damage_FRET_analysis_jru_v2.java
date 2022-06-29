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
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;
import ij.text.*;

public class dna_damage_FRET_analysis_jru_v2 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Acceptor_First",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean a1=gd.getNextBoolean();
		int asp=0; int dsp=1; int aop=2; int dop=3; int atp=4; int dtp=5;
		if(!a1){
			asp=1; dsp=0; aop=3; dop=2; atp=5; dtp=4;
		}
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int maxpts=(int)jstatistics.getstatistic("Max",npts,null);
		int nsets=npts.length/6;
		float[][] rtxvals=new float[nsets][maxpts];
		float[][] rtavals=new float[nsets][maxpts];
		float[][] rtdvals=new float[nsets][maxpts];
		int[] rtnpts=new int[nsets];
		int fretlength=10;
		float[][] stfretvals=new float[nsets][fretlength];
		float[][] offfretvals=new float[nsets][fretlength];
		TextWindow tw=jutils.selectTable("DNA Damage FRET");
		if(tw==null) tw=new TextWindow("DNA Damage FRET","title\tacceptor\tdonor\testripe\teoff\tmaxart\tmaxdrt","",400,200);
		//note that st stands for stripe and off is for areas off the damage stripe
		for(int i=0;i<npts.length/6;i++){
			int len=npts[i*6];
			int damageindex=findbleach(yvals[i*6+dsp],len);
			int fretindex=findbleach(yvals[i*6+atp],len);
			IJ.log("set "+i+" damage pos = "+damageindex+" , fret pos = "+fretindex);
			int predamagestart=damageindex-4;
			int prefretstart=fretindex-5;
			float staccpredam=getavg(yvals[i*6+asp],len,predamagestart,damageindex-1);
			float stdonpredam=getavg(yvals[i*6+dsp],len,predamagestart,damageindex-1);
			float nucaccpredam=getavg(yvals[i*6+atp],len,predamagestart,damageindex-1);
			float nucdonpredam=getavg(yvals[i*6+dtp],len,predamagestart,damageindex-1);
			float stdonprefret=getavg(yvals[i*6+dsp],len,prefretstart,fretindex-2);
			float stdonafret=getavg(yvals[i*6+dsp],len,fretindex,fretindex+3);
			float offdonprefret=getavg(yvals[i*6+dop],len,prefretstart,fretindex-2);
			float offdonafret=getavg(yvals[i*6+dop],len,fretindex,fretindex+3);
			float estripe=1.0f-stdonprefret/stdonafret;
			float eoff=1.0f-offdonprefret/offdonafret;
			rtnpts[i]=len;
			for(int j=0;j<len;j++){
				rtxvals[i][j]=j-damageindex-1;
				rtavals[i][j]=(yvals[i*6+asp][j]/staccpredam)/(yvals[i*6+atp][j]/nucaccpredam);
				rtdvals[i][j]=(yvals[i*6+dsp][j]/stdonpredam)/(yvals[i*6+dtp][j]/nucdonpredam);
			}
			float[] smart=(float[])algutils.get_subarray(rtavals[i],0,fretindex);
			float[] smdrt=(float[])algutils.get_subarray(rtdvals[i],0,fretindex);
			jsmooth.blur1D(smart,2.0f);
			jsmooth.blur1D(smdrt,2.0f);
			float maxart=0.0f;
			float maxdrt=0.0f;
			for(int j=0;j<fretindex-2;j++){
				if(smart[j]>maxart) maxart=smart[j];
				if(smdrt[j]>maxdrt) maxdrt=smdrt[j];
			}
					
			stfretvals[i]=getregion(yvals[i*6+dsp],len,prefretstart,fretlength);
			offfretvals[i]=getregion(yvals[i*6+dop],len,prefretstart,fretlength);
			tw.append(iw.getTitle()+"-"+(i+1)+"\t"+staccpredam+"\t"+stdonpredam+"\t"+estripe+"\t"+eoff+"\t"+maxart+"\t"+maxdrt);
		}
		new PlotWindow4("Stripe_FRET_profiles","time","intensity",stfretvals,null).draw();
		new PlotWindow4("OffStripe_FRET_profiles","time","intensity",offfretvals,null).draw();
		new PlotWindow4("Acc_Rt_profiles","time","intensity",rtxvals,rtavals,rtnpts).draw();
		new PlotWindow4("Don_Rt_profiles","time","intensity",rtxvals,rtdvals,rtnpts).draw();
	}

	public float[] getregion(float[] yvals,int npts,int start,int length){
		float[] temp=new float[length];
		for(int i=start;i<(start+length);i++){
			float val=yvals[0];
			if(i>=0 && i<npts) val=yvals[i];
			if(i>=npts) val=yvals[npts-1];
			temp[i-start]=val;
		}
		return temp;
	}

	public float getavg(float[] yvals,int npts,int start,int end){
		int start2=start;
		if(start2<0) start2=0;
		int end2=end;
		if(end2>=npts) end2=(npts-1);
		float sum=0.0f;
		int count=0;
		for(int i=start2;i<=end2;i++){
			sum+=yvals[i];
			count++;
		}
		return sum/(float)count;
	}

	public int findbleach(float[] yvals,int npts){
		float minder=0.0f;
		int minindex=0;
		float prevval=yvals[0];
		for(int i=1;i<npts;i++){
			float derivative=yvals[i]-prevval;
			if(derivative<minder){
				minder=derivative;
				minindex=i;
			}
			prevval=yvals[i];
		}
		return minindex;
	}

}
