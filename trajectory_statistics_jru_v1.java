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
import jguis.*;

public class trajectory_statistics_jru_v1 implements PlugIn {
	//this plugin calculates statistics for a trajectory
	//if pts to avg is less than the total number of pts, a trajectory of the statistic is created
	//this can be used for N and B analysis
	int avgframes,index;
	boolean avgall,incthird;

	public void run(String arg) {
		GenericDialog gd = new GenericDialog("Options");
		avgframes=500;
		gd.addNumericField("Points to Avg",avgframes,0);
		avgall=false;
		gd.addCheckbox("Avg all points",avgall);
		String[] vartypes={"Variance","Covariance","Temporal Covar"};
		gd.addChoice("Variance Type",vartypes,vartypes[index]);
		gd.addCheckbox("Include third moment",incthird);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		avgframes=(int)gd.getNextNumber();
		avgall=gd.getNextBoolean();
		index=gd.getNextChoiceIndex();
		incthird=gd.getNextBoolean();
		if(index!=0){incthird=false;}

		ImageWindow iw=WindowManager.getCurrentWindow();
		String title=iw.getTitle();
		int selected=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		if(selected<0){selected=0;} if(selected>=yvals.length){selected=0;}
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int length=npts[selected];

		int newlength=(int)((float)length/(float)avgframes);
		if(avgall){newlength=1; avgframes=length;}
		if(index==0){
			double[] second = new double[newlength];
			float[] var=new float[newlength];
			double[] first = new double[newlength];
			float[] avg=new float[newlength];
			double[] third=new double[newlength];
			float[] dthird=new float[newlength];
			float[] K2=new float[newlength];
			float[] B1=new float[newlength];
			float[] K3=new float[newlength];
			float[] B2=new float[newlength];
			for(int i=0;i<newlength;i++){
				for(int j=0;j<avgframes;j++){
					double temp=(double)yvals[selected][j+i*avgframes];
					first[i]+=temp/(double)avgframes;
					second[i]+=(temp*temp)/(double)avgframes;
					if(incthird){third[i]+=(temp*temp*temp)/(double)avgframes;}
				}
				IJ.showProgress(i,newlength);
			}
			for(int i=0;i<newlength;i++){
				avg[i]=(float)first[i];
				var[i]=(float)(second[i]-first[i]*first[i]);
				K2[i]=(float)(second[i]-first[i]*first[i]-first[i]);
				B1[i]=K2[i]/avg[i];
				if(incthird){
					dthird[i]=(float)(third[i]+2.0*avg[i]*avg[i]*avg[i]-3.0*second[i]*avg[i]);
					K3[i]=(float)(third[i]+2.0*avg[i]*avg[i]*avg[i]-3.0*second[i]*avg[i]-3.0*second[i]+3.0*first[i]*first[i]+2.0*first[i]);
					B2[i]=K3[i]/K2[i];
				}
			}
			if(!avgall){
				float[] xvals=new float[newlength];
				for(int i=0;i<newlength;i++){
					xvals[i]=i+1;
				}
				(new PlotWindow4(title+"_avg","Bin","avg",xvals,avg)).draw();
				(new PlotWindow4(title+"_var","Bin","var",xvals,var)).draw();
				(new PlotWindow4(title+"_K2","Bin","K2",xvals,K2)).draw();
				(new PlotWindow4(title+"_B1","Bin","K2/K1",xvals,B1)).draw();
				if(incthird){
					(new PlotWindow4(title+"_third","Bin","third",xvals,dthird)).draw();
					(new PlotWindow4(title+"_K3","Bin","K3",xvals,K3)).draw();
					(new PlotWindow4(title+"_B2","Bin","K3/K2",xvals,B2)).draw();
				}
			} else {
				IJ.log("Statistics for "+title);
				IJ.log("Avg = "+avg[0]);
				IJ.log("Var = "+var[0]);
				IJ.log("K2 = "+K2[0]);
				IJ.log("K2/K1 = "+B1[0]);
				if(incthird){
					IJ.log("Third Central Moment = "+dthird[0]);
					IJ.log("K3 = "+K3[0]);
					IJ.log("K3/K2 = "+B2[0]);
				}
			}
		}
		if(index==1){
			double[] first1 = new double[newlength];
			float[] avg1=new float[newlength];
			double[] first2 = new double[newlength];
			float[] avg2=new float[newlength];
			double[] second1=new double[newlength];
			float[] var1=new float[newlength];
			double[] second2=new double[newlength];
			float[] var2=new float[newlength];
			double[] second12=new double[newlength];
			float[] covar=new float[newlength];
			float[] g0cc=new float[newlength];
			for(int i=0;i<newlength;i++){
				for(int j=0;j<avgframes;j++){
					double temp1=(double)yvals[selected][j+i*avgframes];
					double temp2=(double)yvals[selected+1][j+i*avgframes];
					first1[i]+=temp1/(double)avgframes;
					second1[i]+=(temp1*temp1)/(double)avgframes;
					first2[i]+=temp2/(double)avgframes;
					second1[i]+=(temp1*temp1)/(double)avgframes;
					second12[i]+=(temp1*temp2)/(double)avgframes;
				}
				IJ.showProgress(i,newlength);
			}
			for(int i=0;i<newlength;i++){
				avg1[i]=(float)first1[i];
				avg2[i]=(float)first2[i];
				var1[i]=(float)(second1[i]-first1[i]*first1[i]);
				var2[i]=(float)(second2[i]-first2[i]*first2[i]);
				covar[i]=(float)(second12[i]-first1[i]*first2[i]);
				g0cc[i]=covar[i]/(avg1[i]*avg2[i]);
			}
			if(!avgall){
				float[] xvals=new float[newlength];
				for(int i=0;i<newlength;i++){
					xvals[i]=i+1;
				}
				(new PlotWindow4(title+"_avg1","Bin","avg1",xvals,avg1)).draw();
				(new PlotWindow4(title+"_avg2","Bin","avg2",xvals,avg2)).draw();
				(new PlotWindow4(title+"_var1","Bin","var1",xvals,var1)).draw();
				(new PlotWindow4(title+"_var2","Bin","var2",xvals,var2)).draw();
				(new PlotWindow4(title+"_covar","Bin","covar",xvals,covar)).draw();
				(new PlotWindow4(title+"_G0cc","Bin","G(0)cc",xvals,g0cc)).draw();
			} else {
				IJ.log("Statistics for "+title+" "+selected+" and "+(selected+1));
				IJ.log("Avg1 = "+avg1[0]);
				IJ.log("Var1 = "+var1[0]);
				IJ.log("Avg2 = "+avg2[0]);
				IJ.log("Var2 = "+var2[0]);
				IJ.log("covar = "+covar[0]);
				IJ.log("G(0)cc = "+g0cc[0]);
			}
		}
		/*if(index==2){
			float[] var = new float[width*newheight];
			float[] avg = new float[width*newheight];
			if(ip instanceof FloatProcessor){
				float[] pixels = (float[])ip.getPixels();
				for(int i=0;i<newheight;i++){
					for(int j=0;j<(avgframes-1);j++){
						for(int k=0;k<width;k++){
							float temp1=pixels[(i*avgframes+j)*width+k];
							float temp2=pixels[(i*avgframes+j+1)*width+k];
							avg[k+i*width]+=(temp1+temp2)/(float)(2*avgframes);
							var[k+i*width]+=(temp1*temp2)/(float)avgframes;
						}
					}
					IJ.showProgress(i,newheight);
				}
			}
			if(ip instanceof ShortProcessor){
				short[] pixels = (short[])ip.getPixels();
				for(int i=0;i<newheight;i++){
					for(int j=0;j<(avgframes-1);j++){
						for(int k=0;k<width;k++){
							float temp1=pixels[(i*avgframes+j)*width+k]&0xffff;
							float temp2=pixels[(i*avgframes+j+1)*width+k]&0xffff;
							avg[k+i*width]+=(temp1+temp2)/(float)(2*avgframes);
							var[k+i*width]+=(temp1*temp2)/(float)avgframes;
					}
					}
					IJ.showProgress(i,newheight);
				}
			}
			for(int i=0;i<newheight*width;i++){
				var[i]-=avg[i]*avg[i];
			}
			if(!avgall){
				(new ImagePlus(trunctitle+"_avg",new FloatProcessor(width,newheight,avg,null))).show();
				(new ImagePlus(trunctitle+"_g1var",new FloatProcessor(width,newheight,var,null))).show();
			} else {
				float[] xvals=new float[width];
				for(int i=0;i<width;i++){
					xvals[i]=i+1;
				}
				(new PlotWindow(trunctitle+"_avg","Pixel","avg",xvals,avg)).draw();
				(new PlotWindow(trunctitle+"_g1var","Pixel","g1var",xvals,var)).draw();
			}
		}*/
	}


}
