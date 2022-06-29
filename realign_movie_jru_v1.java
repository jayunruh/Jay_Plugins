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
import jalgs.*;

public class realign_movie_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length+1];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){
				titles[i]=imp.getTitle();
			} else {
				titles[i]="NA";
			}
		}
		titles[wList.length]="null";
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addChoice("Image",titles,titles[0]);
		gd2.addChoice("Track",titles,titles[0]);
		gd2.addChoice("Angles (interpolation only, degrees)",titles,titles[wList.length]);
		gd2.addCheckbox("Interpolate?",true);
		gd2.addCheckbox("Difference_Track?",false);
		gd2.showDialog();
		if(gd2.wasCanceled()){return;}
		int index1 = gd2.getNextChoiceIndex();
		int index2 = gd2.getNextChoiceIndex();
		int index3=gd2.getNextChoiceIndex();
		boolean interp=gd2.getNextBoolean();
		boolean difftrack=gd2.getNextBoolean();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		int frames=imp.getNFrames();
		if(frames==1){
			frames=slices;
			slices=1;
		}
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		ImageWindow iw=imp2.getWindow();
		float[] xvals1=((float[][])jutils.runPW4VoidMethod(iw,"getXValues"))[0];
		float[] yvals1=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
		float[] angles1=new float[xvals1.length];
		if(index3<wList.length) angles1=((float[][])jutils.runPW4VoidMethod(WindowManager.getImage(wList[index3]).getWindow(),"getYValues"))[0];
		float[] xvals=new float[xvals1.length];
		float[] yvals=new float[xvals1.length];
		float[] angles=new float[xvals1.length];
		if(!difftrack){
			for(int i=0;i<xvals.length;i++){
				xvals[i]=-xvals1[i];
				yvals[i]=-yvals1[i];
				angles[i]=-angles1[i];
			}
		} else {
			xvals[0]=0.0f;
			yvals[0]=0.0f;
			angles[0]=0.0f;
			for(int i=1;i<xvals.length;i++){
				xvals[i]=xvals[i-1]+(xvals1[i]-xvals1[0]);
				yvals[i]=yvals[i-1]+(yvals1[i]-yvals1[0]);
				angles[i]=angles[i-1]+(angles1[i]-angles1[0]);
			}
		}
		float xmin=xvals[0]; float xmax=xvals[0];
		float ymin=yvals[0]; float ymax=yvals[0];
		for(int i=1;i<xvals.length;i++){
			if(xvals[i]<xmin){xmin=xvals[i];}
			if(xvals[i]>xmax){xmax=xvals[i];}
			if(yvals[i]<ymin){ymin=yvals[i];}
			if(yvals[i]>ymax){ymax=yvals[i];}
		}
		int intxmax=1+(int)xmax;
		int intxmin=(int)xmin-1;
		int intymax=1+(int)ymax;
		int intymin=(int)ymin-1;
		int newwidth=width+intxmax-intxmin+1;
		int newheight=height+intymax-intymin+1;
		float newcenterx=xvals1[0]+xvals[0]-(float)intxmin;
		float newcentery=yvals1[0]+yvals[0]-(float)intymin;
		//IJ.log(""+newcenterx+" , "+newcentery);
		ImageStack retstack=new ImageStack(newwidth,newheight);
		for(int k=0;k<frames;k++){
			for(int l=0;l<(slices*channels);l++){
				Object currimage=stack.getPixels(1+l+k*slices*channels);
				Object tempimage;
				int typeindex;
				if(currimage instanceof float[]){
					tempimage=new float[newwidth*newheight];
					typeindex=0;
				} else {
					if(currimage instanceof short[]){
						tempimage=new short[newwidth*newheight];
						typeindex=1;
					} else {
						tempimage=new byte[newwidth*newheight];
						typeindex=2;
					}
				}
				if(interp){
					float[] tempimage2=interpolation.shift_image(currimage,width,height,xvals[k]-(float)intxmin,yvals[k]-(float)intymin);
					if(angles[k]!=0.0f) tempimage2=interpolation.rotate_image(tempimage2,width,height,-(float)Math.toRadians(angles[k]),newcenterx,newcentery);
					float[] tempimage3=new float[newwidth*newheight];
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							tempimage3[j+i*newwidth]=tempimage2[j+i*width];
						}
					}
					tempimage=algutils.convert_array(tempimage3,2-typeindex);
				} else {
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							//need to add rotation here
							int index5=j+(int)xvals[k]-intxmin+1+(i+(int)yvals[k]-intymin+1)*newwidth;
							int index6=j+i*width;
							if(typeindex==0){((float[])tempimage)[index5]=((float[])currimage)[index6];}
							if(typeindex==1){((short[])tempimage)[index5]=((short[])currimage)[index6];}
							if(typeindex==2){((byte[])tempimage)[index5]=((byte[])currimage)[index6];}
						}
					}
				}
				retstack.addSlice("",tempimage);
			}
			IJ.showProgress(k,frames);
		}
		ImagePlus imp3=new ImagePlus("Realigned",retstack);
		imp3.copyScale(imp);
		if(slices*channels>1){
			imp3.setOpenAsHyperStack(true);
			imp3.setDimensions(channels,slices,frames);
			if(channels>1 && imp.isComposite()){
				CompositeImage ci=new CompositeImage(imp3,((CompositeImage)imp).getMode());
				ci.copyLuts(imp);
				ci.show();
			} else {
				imp3.show();
			}
		} else {
			imp3.show();
		}
	}

}
