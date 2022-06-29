/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jsim.*;
import jguis.*;

public class random_number_generator_jru_v1 implements PlugIn {

	public void run(String arg) {
		int nnumbers=100;
		String[] types={"Uniform","Gaussian","Poisson","Binary","Exponential","Rand_Order"};
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Number_Type",types,types[0]);
		gd2.addNumericField("Number_of_values",nnumbers,0);
		gd2.showDialog();  if(gd2.wasCanceled()) return;
		int type=gd2.getNextChoiceIndex();
		nnumbers=(int)gd2.getNextNumber();
		rngs random=new rngs();
		float[] nums=new float[nnumbers];
		if(type==5){
			int[] temp=random.random_order(nnumbers);
			nums=algutils.convert_arr_float(temp);
			new PlotWindow4("Random Numbers","sample","number",nums).draw();
			return;
		}
		float[] extras=new float[3];
		extras[0]=10.0f;
		GenericDialog gd=new GenericDialog("Options");
		if(type==0){
			gd.addNumericField("Lower Limit",0.0,5,15,null);
			gd.addNumericField("Upper Limit",1.0,5,15,null);
		}
		if(type==1){
			gd.addNumericField("Avg",0.0,5,15,null);
			gd.addNumericField("StDev",1.0,5,15,null);
		}
		if(type==2 || type==4){
			gd.addNumericField("Avg",0.0,5,15,null);
		}
		if(type==3){
			gd.addNumericField("Prob_of_success",0.5,5,15,null);
		}
		gd.showDialog(); if(gd.wasCanceled()) return;
		extras[0]=(float)gd.getNextNumber();
		if(type<2) extras[1]=(float)gd.getNextNumber();
		for(int i=0;i<100;i++){
			if(type==0) nums[i]=(float)random.unidev((double)extras[1],(double)extras[0]);
			if(type==1) nums[i]=(float)random.gasdev((double)extras[1],(double)extras[0]);
			if(type==2) nums[i]=(float)random.poidev((double)extras[0]);
			if(type==3) nums[i]=((float)random.unidev(1.0,0.0)>extras[0])?0.0f:1.0f;
			if(type==4) nums[i]=(float)random.expdev((double)extras[0]);
		}
		new PlotWindow4("Random Numbers","sample","number",nums).draw();
	}

}
