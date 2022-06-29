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
import jguis.*;
import jalgs.jseg.*;
import jalgs.*;

public class KISS_analysis_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] labels={"Masked_Chromosomes","Unmixed_Image","Spectral_Image(optional)","Spectra(optional)"};
		ImagePlus[] imps=jutils.selectImages(true,4,labels);
		if(imps==null){return;}
		if(imps[0]==null){return;}
		float[] mask=(float[])imps[0].getStack().getPixels(2);
		findblobs3 fb=new findblobs3(imps[0].getWidth(),imps[0].getHeight());
		float[] objects=fb.dofindblobs(mask,0.5f);
		WaitForUserDialog dg=new WaitForUserDialog("Optional Input","Place RoiManager Points on Chromosome Segments (if desired)");
		dg.show();
		if(!dg.escPressed()){
			RoiManager rman=RoiManager.getInstance();
			while(rman!=null && rman.getCount()>1){
				Roi[] rois=rman.getRoisAsArray();
				int[] ids=new int[rois.length];
				for(int i=0;i<rois.length;i++){
					Rectangle r=rois[i].getBounds();
					ids[i]=(int)objects[r.x+fb.width*r.y];
				}
				objects=fb.link_objects(objects,ids);
				rman.reset();
				dg=new WaitForUserDialog("Optional Input","Place More RoiManager Points on Chromosome Segments (if desired)");
				dg.show();
				if(dg.escPressed()) break;
			}
		}
		int[] areas=fb.get_areas(objects);
		int[] temprank=jsort.get_javasort_order(areas);
		int[] arearank=jsort.get_javasort_order(temprank);
		for(int i=0;i<fb.nobjects;i++){
			arearank[i]=fb.nobjects-arearank[i]-1;
		}
		//if the spectra are available, get them
		float[][][] spectra=null;
		Object[] data=null;
		if(imps[1]!=null && imps[2]!=null && imps[3]!=null){
			ImageWindow iw=imps[3].getWindow();
			if(iw.getClass().getName().equals("jguis.PlotWindow4")){
				float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
				data=jutils.stack2array(imps[2].getStack());
				Object[] coef=jutils.stack2array(imps[1].getStack());
				spectra=new float[fb.nobjects][2][];
				for(int i=0;i<fb.nobjects;i++){
					spectra[i][0]=fb.get_object_spectrum(objects,(i+1),data,"Sum");
					spectra[i][1]=new float[yvals[0].length];
					float[] tempcoef=fb.get_object_spectrum(objects,(i+1),coef,"Sum");
					for(int j=0;j<yvals[0].length;j++){
						for(int k=0;k<5;k++){
							spectra[i][1][j]+=tempcoef[k]*yvals[k][j];
						}
					}
				}
			}
		}
		CompositeImage imp=(CompositeImage)imps[0];
		imp.setPosition(1,1,1);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		imp.setChannelColorModel(graylut);
		imp.setPosition(2,1,1);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		imp.setChannelColorModel(redlut);
		imp.setPosition(1,1,1);
		imp.updateAndRepaintWindow();
		KissPanel sp=new KissPanel();
		int skychan=6;
		if(imps[1]!=null) skychan=imps[1].getNChannels();
		//assume that the sky image has 6 channels and that the second is the unknown green
		//shift the unknown green to the end
		ImagePlus skyimp=null;
		if(imps[1]!=null){
			Object[] skystack=jutils.stack2array(imps[1].getStack());
			//Object[] skystack2={skystack[0],skystack[2],skystack[3],skystack[4],skystack[5],skystack[1]};
			Object[] skystack2=null;
			if(skychan==6) skystack2=new Object[]{skystack[0],skystack[2],skystack[3],skystack[4],skystack[5]};
			else skystack2=new Object[]{skystack[0],skystack[1],skystack[2],skystack[3],skystack[4]};
			skyimp=new ImagePlus("rearranged",jutils.array2stack(skystack2,imps[1].getWidth(),imps[1].getHeight()));
		}
		int nch=5;
		if(skyimp!=null) nch=skyimp.getStack().getSize();
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addNumericField("Area Accuracy (percent)",30,0);
		for(int i=0;i<nch;i++){
			gd2.addNumericField("Ch_"+(i+1)+"_Contr_Thresh",0.35,5,15,null);
		}
		//gd2.addNumericField("Contribution Threshold",0.35,5,15,null);
		gd2.addCheckbox("Mouse?",false);
		Object[] codes=KissPanel.getCustomCodes();
		String[] codenames=new String[]{"none"};
		if(codes!=null){
			String[] temp=(String[])codes[0];
			codenames=new String[temp.length+1]; codenames[0]="none";
			for(int i=0;i<temp.length;i++) codenames[i+1]=temp[i];
		}
		gd2.addChoice("Custom_Code",codenames,codenames[0]);
		gd2.addNumericField("Box_Width",150,0);
		gd2.addNumericField("Box_Height",100,0);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		sp.areathresh=(float)gd2.getNextNumber();
		sp.objthresh2=new float[nch];
		for(int i=0;i<nch;i++) sp.objthresh2[i]=(float)gd2.getNextNumber();
		//sp.objthresh=(float)gd2.getNextNumber();
		boolean mouse=gd2.getNextBoolean();
		int codeindex=gd2.getNextChoiceIndex();
		int bwidth=(int)gd2.getNextNumber();
		int bheight=(int)gd2.getNextNumber();
		int[] colorindices={4,1,2,6,3};
		GenericDialog gd3=new GenericDialog("Color Options");
		for(int i=0;i<5;i++) gd3.addChoice("Ch"+(i+1)+" Color",KissPanel.colornames,KissPanel.colornames[colorindices[i]]);
		gd3.showDialog(); if(gd3.wasCanceled()) return;
		for(int i=0;i<5;i++) colorindices[i]=gd3.getNextChoiceIndex();
		sp.colorindices=colorindices;
		sp.nch=5;
		sp.dapilast=false;
		sp.cellwidth=bwidth;
		sp.cellheight=bheight;
		int[][] custcode=null;
		if(codeindex>0) custcode=(int[][])codes[codeindex+1];
		sp.init(imps[0],skyimp,objects,areas,arearank,fb,true,spectra,data,mouse,custcode);
		KissPanel.launch_frame(sp);
	}

}
