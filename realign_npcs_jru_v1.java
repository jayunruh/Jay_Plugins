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
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class realign_npcs_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin realigns all Npc's in an (approximately spherical) nucleus so that we are viewing down the pore
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Image","Trajectory"});
		if(imps==null) return;
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		Object[] stack=jutils.stack2array(imps[0].getStack());
		int nchans=imps[0].getNChannels();
		int nslices=imps[0].getNSlices();
		if(nslices==1) nslices=stack.length/nchans;
		Object[][] stack2=new Object[nchans][nslices];
		for(int i=0;i<nchans;i++){
			for(int j=0;j<nslices;j++){
				stack2[i][j]=stack[j*nchans+i];
			}
		}
		ImageWindow iw=imps[1].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
		float[][] zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];
		float zratio=(float)jutils.get_zratio(imps[0]);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_ratio",zratio,5,15,null);
		gd.addNumericField("Reconstruction_Size (cube)",10,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		zratio=(float)gd.getNextNumber();
		int recsize=(int)gd.getNextNumber();
		//start by getting the trajectory center of mass (assumed to be nuclear com)
		//expect each point to be its own trajectory
		float[] com=new float[3];
		for(int i=0;i<xvals.length;i++){
			com[0]+=xvals[i][0];
			com[1]+=yvals[i][0];
			com[2]+=zvals[i][0];
		}
		com[0]/=(float)xvals.length;
		com[1]/=(float)xvals.length;
		com[2]/=(float)xvals.length;
		IJ.log("com = "+com[0]+" , "+com[1]+" , "+com[2]);
		//now get the unit vector for each npc through the com
		float[][] realigned=new float[nchans*recsize*xvals.length][];
		float[][] angles=new float[xvals.length][3];
		String[] labels=new String[xvals.length*recsize*nchans];
		String[][] tabledata=new String[xvals.length][3]; 
		for(int i=0;i<xvals.length;i++){
			float[] vec={xvals[i][0]-com[0],yvals[i][0]-com[1],zvals[i][0]-com[2]};
			float[] center={xvals[i][0],yvals[i][0],zvals[i][0]};
			vec=measure_object.norm_vector(vec);
			//get the inner angle to the z axis
			float angle=measure_object.get_inner_angle(vec,new float[]{0.0f,0.0f,1.0f});
			angles[i][1]=angle;
			angles[i][1]-=0.5f*(float)Math.PI; //shift the angle so that equatorial is 0, up is -90 and down is 90
			angles[i][1]*=180.0f/(float)Math.PI; //convert to degrees
			angles[i][0]=i+1;
			//get the xy projection angle to the x axis
			//float[] projvec={vec[0],vec[1],0.0f};
			//projvec=measure_object.norm_vector(projvec);
			//angles[i][2]=measure_object.get_inner_angle(projvec,new float[]{1.0f,0.0f,0.0f});
			angles[i][2]=(float)Math.atan2(vec[1],vec[0]);
			angles[i][2]*=180.0f/(float)Math.PI; //convert to degrees
			tabledata[i][0]=""+(int)angles[i][0];
			tabledata[i][1]=""+angles[i][1];
			tabledata[i][2]=""+angles[i][2];
			for(int j=0;j<nchans;j++){
				float[][] temp=profiler.getRotated3DImage(stack2[j],width,height,center,zratio,vec,recsize,recsize/2,recsize);
				for(int k=0;k<recsize;k++){
					realigned[i*nchans*recsize+k*nchans+j]=temp[k];
					labels[i*nchans*recsize+k*nchans+j]=""+(int)angles[i][0];
				}
			}
		}
		ImagePlus recimp=jutils.create_hyperstack("Realigned_Npcs",jutils.array2stack(realigned,recsize,recsize),xvals.length,recsize,nchans,true,null);
		ImageStack recstack=recimp.getStack();
		for(int i=0;i<recstack.getSize();i++) recstack.setSliceLabel(labels[i],i+1);
		recimp.show();
		table_tools.create_table("Npc_Angles",tabledata,new String[]{"Npc","Inclination_Angle","Azimuthal_Angle"});
	}

}
