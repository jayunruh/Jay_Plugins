/*******************************************************************************
 * Copyright (c) 2020 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import jguis.*;
import java.util.*;

public class SIFTj_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Initial_Blur (pix)",1.6,5,15,null);
		gd.addNumericField("Steps_per_octave",3,0);
		gd.addNumericField("Min_Octave_Size",64,0);
		gd.addNumericField("Max_Octave_Size",1024,0);
		gd.addNumericField("Feature_Descriptor_Size",4,0);
		gd.addNumericField("Feature_Orientation_Bins",8,0);
		gd.addChoice("Transformation",SIFTj.modelNames,SIFTj.modelNames[0]);
		gd.addCheckbox("Interpolate",false);
		gd.addCheckbox("Align_Secondary",false);
		gd.addNumericField("N_Threads",1,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float[] params={1.6f,3.0f,64.0f,1024.0f,4.0f,8.0f,0.0f};
		params[0]=(float)gd.getNextNumber();
		params[1]=(float)gd.getNextNumber();
		params[2]=(float)gd.getNextNumber();
		params[3]=(float)gd.getNextNumber();
		params[4]=(float)gd.getNextNumber();
		params[5]=(float)gd.getNextNumber();
		params[6]=(float)gd.getNextChoiceIndex();
		boolean interp=gd.getNextBoolean();
		boolean secondary=gd.getNextBoolean();
		int nthreads=(int)gd.getNextNumber();
		ImagePlus simp=null;
		if(secondary){
			ImagePlus[] temp=jutils.selectImages(false,1,new String[]{"Secondary_Image"});
			if(temp==null) return;
			simp=temp[0];
		}
		SIFTj siftj=new SIFTj(params,interp);
		Object[] output=null;
		if(nthreads<2) output=siftj.runSIFT(imp,simp);
		else output=siftj.runSIFT(imp,simp,nthreads);
		ImagePlus outimp=(ImagePlus)output[0];
		outimp.show();
		//String[] outtrans=(String[])output[1];
		double[][] outtrans=(double[][])output[1];
		List<List<String>> listtable=new ArrayList<List<String>>();
		for(int i=0;i<outtrans.length;i++){
			//IJ.log(outtrans[i]);
			//String tempr=outtrans[i].replace('[',',');
			//tempr=tempr.replace(']',',');
			//String[] split1=table_tools.split(tempr,",",true);
			List<String> row=new ArrayList<String>();
			row.add(""+(i+1));
			//for(int j=1;j<split1.length;j++){
			//	if(!split1[j].trim().isEmpty()) row.add(split1[j]);
			//}
			double[] temp={outtrans[i][0],outtrans[i][2],outtrans[i][4],outtrans[i][1],outtrans[i][3],outtrans[i][5]};
			for(int j=0;j<6;j++) row.add(""+temp[j]);
			listtable.add(row);
		}
		List<String> col_labels=new ArrayList<String>();
		col_labels.add("frame");
		for(int i=1;i<listtable.get(1).size();i++) col_labels.add("col"+i);
		table_tools.create_table("Transformations",listtable,col_labels);
		float[][] translation=(float[][])output[2];
		new PlotWindow4("SIFT Translation","x","y",translation[0],translation[1]).draw();
		if(simp!=null){
			((ImagePlus)output[3]).show();
		}
	}

}
