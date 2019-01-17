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
import jguis.*;

public class traj_calculator_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] operations=new String[7];
		operations[0]="add";
		operations[1]="subtract";
		operations[2]="multiply";
		operations[3]="divide_by";
		operations[4]="to_the_power";
		operations[5]="sub_from";
		operations[6]="divided_from";

		Object[] windowList=jutils.getPlotWindowList(false);
		String[] titles=(String[])windowList[1];
		int[] ids=(int[])windowList[0];
		//for(int i=0;i<ids.length;i++) IJ.log(WindowManager.getImage(ids[i]).getWindow().getClass().getName());
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Trajectory_1",titles,titles[0]);
		gd.addChoice("Operation",operations,operations[0]);
		gd.addChoice("Trajectory_2",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int opindex=gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		ImagePlus imp = WindowManager.getImage(ids[index1]);
		ImagePlus imp2 = WindowManager.getImage(ids[index2]);

		ImageWindow iw=imp.getWindow();
		ImageWindow iw2=imp2.getWindow();
		float[] xvals1,xvals2,yvals1,yvals2;
		String title1,title2;
		int length1,length2;
		if(iw!=iw2){
			title1=imp.getTitle();
			title2=imp2.getTitle();
			xvals1=((float[][])jutils.runPW4VoidMethod(iw,"getXValues"))[0];
			yvals1=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
			xvals2=((float[][])jutils.runPW4VoidMethod(iw2,"getXValues"))[0];
			yvals2=((float[][])jutils.runPW4VoidMethod(iw2,"getYValues"))[0];
			length1=xvals1.length;
			length2=xvals2.length;
			if(length1!=length2){IJ.error("trajectories are not the same length"); return;}
		} else {
			float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			xvals1=xvals[0]; yvals1=yvals[0]; xvals2=xvals[1]; yvals2=yvals[1];
			length1=yvals1.length; length2=yvals2.length;
			title1=imp.getTitle()+1;
			title2=imp.getTitle()+2;
		}
		float[] retyvals=new float[length1];
		float[] retxvals=new float[length1];
		for(int i=0;i<length1;i++){retxvals[i]=xvals1[i];}
		if(opindex==0){for(int i=0;i<length1;i++){retyvals[i]=yvals1[i]+yvals2[i];}}
		if(opindex==1){for(int i=0;i<length1;i++){retyvals[i]=yvals1[i]-yvals2[i];}}
		if(opindex==2){for(int i=0;i<length1;i++){retyvals[i]=yvals1[i]*yvals2[i];}}
		if(opindex==3){for(int i=0;i<length1;i++){retyvals[i]=yvals1[i]/yvals2[i];}}
		if(opindex==4){for(int i=0;i<length1;i++){retyvals[i]=(float)Math.pow(yvals1[i],yvals2[i]);}}
		if(opindex==5){for(int i=0;i<length1;i++){retyvals[i]=yvals2[i]-yvals1[i];}}
		if(opindex==6){for(int i=0;i<length1;i++){retyvals[i]=yvals2[i]/yvals1[i];}}
		String operator="";
		if(opindex==0){operator="+";}
		if(opindex==1 || opindex==5){operator="-";}
		if(opindex==2){operator="*";}
		if(opindex==3 || opindex==6){operator="/";}
		if(opindex==4){operator="^";}
		PlotWindow4 newpw=null;
		if(opindex<5) newpw=new PlotWindow4(title1+operator+title2,"x","y",retxvals,retyvals);
		else newpw=new PlotWindow4(title2+operator+title1,"x","y",retxvals,retyvals);
		newpw.draw();
	}

}
