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

public class traj_math_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] operations=new String[11];
		operations[0]="add";
		operations[1]="subtract";
		operations[2]="multiply";
		operations[3]="divide by";
		operations[4]="divide op by";
		operations[5]="set min";
		operations[6]="set max";
		operations[7]="log base";
		operations[8]="to the power";
		operations[9]="power of";
		operations[10]="abs";
		GenericDialog gd = new GenericDialog("Parameters");
		boolean opx=false;
		gd.addCheckbox("Operate_on_X",opx);
		gd.addChoice("operation",operations,operations[0]);
		float operator=0.0f;
		gd.addNumericField("operator",operator,5,20,"");
		boolean opall=false;
		gd.addCheckbox("Analyze_All_Series",opall);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		opx=gd.getNextBoolean();
		int opindex=gd.getNextChoiceIndex();
		operator=(float)gd.getNextNumber();
		opall=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=null;
		if(opall){
			pw=jutils.getPW4Copy(iw);
		} else {
			pw=jutils.getPW4SelCopy(iw);
		}
		String title=pw.getTitle();
		pw.getImagePlus().setTitle("Result of "+title);
		float[][] xvals=pw.getXValues();
		float[][] yvals=pw.getYValues();
		if(opall){
			for(int i=0;i<yvals.length;i++){
				if(opx) operate(xvals[i],operator,opindex);
				else operate(yvals[i],operator,opindex);
			}
		} else {
			if(opx) operate(xvals[0],operator,opindex);
			else operate(yvals[0],operator,opindex);
		}
		pw.autoscale();
	}

	public void operate(float[] opvals,float operator,int opindex){
		int length=opvals.length;
		if(opindex==0){for(int i=0;i<length;i++){opvals[i]+=operator;}}
		if(opindex==1){for(int i=0;i<length;i++){opvals[i]-=operator;}}
		if(opindex==2){for(int i=0;i<length;i++){opvals[i]*=operator;}}
		if(opindex==3){for(int i=0;i<length;i++){opvals[i]/=operator;}}
		if(opindex==4){for(int i=0;i<length;i++){opvals[i]=operator/opvals[i];}}
		if(opindex==5){for(int i=0;i<length;i++){if(opvals[i]<operator){opvals[i]=operator;}}}
		if(opindex==6){for(int i=0;i<length;i++){if(opvals[i]>operator){opvals[i]=operator;}}}
		if(opindex==7){for(int i=0;i<length;i++){opvals[i]=(float)Math.log(opvals[i])/(float)Math.log(operator);}}
		if(opindex==8){for(int i=0;i<length;i++){opvals[i]=(float)Math.pow(opvals[i],operator);}}
		if(opindex==9){for(int i=0;i<length;i++){opvals[i]=(float)Math.pow(operator,opvals[i]);}}
		if(opindex==10){for(int i=0;i<length;i++){opvals[i]=Math.abs(opvals[i]);}}
	}

}
