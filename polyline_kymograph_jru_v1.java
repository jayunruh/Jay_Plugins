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
import ij.measure.*;
import jalgs.*;
import jguis.*;

public class polyline_kymograph_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi thisroi=imp.getRoi();
		if(!(thisroi instanceof PolygonRoi || thisroi instanceof Line)){IJ.showMessage("Linear Selection Required"); return;}
		int linewidth=(int)thisroi.getStrokeWidth();
		if(linewidth<1) linewidth=1;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Profile_Width",linewidth,0);
		boolean connected=false;
		gd.addCheckbox("Circular?",connected);
		boolean cumulative=false;
		gd.addCheckbox("Cumulative?",cumulative);
		String[] orientation={"Centered","Outside","Inside"};
		gd.addChoice("Positioning",orientation,orientation[0]);
		gd.addCheckbox("Output_Straightened",false);
		gd.addCheckbox("Single_Frame_Profile",false);
		gd.addCheckbox("All_Channels (for single frame)",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		linewidth=(int)gd.getNextNumber();
		connected=gd.getNextBoolean();
		cumulative=gd.getNextBoolean();
		int or_index=gd.getNextChoiceIndex();
		boolean straighten=gd.getNextBoolean();
		boolean single=gd.getNextBoolean();
		boolean both=gd.getNextBoolean();

		int width=imp.getWidth();
		int height=imp.getHeight();
		int nchannels=imp.getNChannels();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		if(frames==1){
			frames=slices;
			slices=1;
		}
		if(single){
			frames=1; 
			if(!both) nchannels=1;
		}

		int[] xvals,yvals;
		int nlines,npts;
		if(thisroi instanceof Line){
			xvals=new int[2]; xvals[0]=((Line)thisroi).x1; xvals[1]=((Line)thisroi).x2;
			yvals=new int[2]; yvals[0]=((Line)thisroi).y1; yvals[1]=((Line)thisroi).y2;
			npts=2;
			nlines=1;
		} else {
			PolygonRoi roi=(PolygonRoi)thisroi;
			Rectangle r = roi.getBounds();
			int[] xvals1=roi.getXCoordinates();
			int[] yvals1=roi.getYCoordinates();
			npts=roi.getNCoordinates();
			xvals=new int[npts];
			yvals=new int[npts];
			nlines=npts-1;
			for(int i=0;i<npts;i++){
				xvals[i]=xvals1[i]+r.x;
				yvals[i]=yvals1[i]+r.y;
			}
		}
		Polygon polyroi=new Polygon(xvals,yvals,npts);
		//for(int i=0;i<xvals.length;i++) IJ.log(""+xvals[i]+" , "+yvals[i]);
		//IJ.log(""+profiler.getPolygonLength(polyroi,false));
		
		int length=0;
		boolean first=true;

		//ImageStack kstack=new ImageStack(length,frames);
		ImageStack kstack=null;
		if(!single){
			for(int s=0;s<slices;s++){
				for(int ch=0;ch<nchannels;ch++){
					float[] kymograph=null;
					if(!first) kymograph=new float[length*frames];
					for(int t=0;t<frames;t++){
						Object pixels=stack.getPixels(ch+s*nchannels+t*nchannels*slices+1);
						float[] temp=profiler.getProfile(pixels,width,height,polyroi,connected,linewidth,or_index);
						if(temp==null) return;
						if(first){
							length=temp.length;
							kstack=new ImageStack(length,frames);
							kymograph=new float[length*frames];
							first=false;
						}
						System.arraycopy(temp,0,kymograph,t*length,length);
					}
					if(cumulative){
						for(int i=0;i<length*frames;i++){
							kymograph[i]*=(float)linewidth;
						}
					}
					kstack.addSlice("",(Object)kymograph);
				}
			}
			if(straighten){
				ImageStack sstack=new ImageStack(linewidth,length);
				for(int t=0;t<frames;t++){
					for(int s=0;s<slices;s++){
						for(int ch=0;ch<nchannels;ch++){
							Object pixels=stack.getPixels(ch+s*nchannels+t*nchannels*slices+1);
							float[] temp=profiler.getStraightened(pixels,width,height,polyroi,connected,linewidth,or_index);
							sstack.addSlice("",temp);
						}
					}
				}
				ImagePlus simp=new ImagePlus("Straightened",sstack);
				simp.copyScale(imp);
				simp.setOpenAsHyperStack(true);
				simp.setDimensions(nchannels,slices,frames);
				if(nchannels>1 && imp.isComposite()){
					CompositeImage ci=new CompositeImage(simp,CompositeImage.COMPOSITE);
					ci.copyLuts(imp);
					ci.show();
				} else {
					simp.show();
				}
			}
		} else {
			if(!both || nchannels==1){
				Object pixels=imp.getProcessor().getPixels();
				float[] temp=profiler.getProfile(pixels,width,height,polyroi,connected,linewidth,or_index);
				length=temp.length;
				kstack=new ImageStack(length,1);
				if(cumulative){
					for(int i=0;i<length;i++) temp[i]*=(float)linewidth;
				}
				kstack.addSlice("",temp);
				if(straighten){
					float[] temp2=profiler.getStraightened(pixels,width,height,polyroi,connected,linewidth,or_index);
					ImagePlus simp=new ImagePlus("Straightened",new FloatProcessor(linewidth,length,temp2,null));
					simp.copyScale(imp);
					simp.show();
				}
			} else {
				int frame=imp.getT();
				int slice=imp.getZ();
				Object[] pixels=jutils.get3DCSeries(stack,slice-1,frame-1,frames,slices,nchannels);
				Object[] temp2=new Object[pixels.length];
				for(int i=0;i<pixels.length;i++){
					float[] temp=profiler.getProfile(pixels[i],width,height,polyroi,connected,linewidth,or_index);
					length=temp.length;
					if(i==0) kstack=new ImageStack(length,1);
					if(cumulative){
						for(int j=0;j<length;j++) temp[j]*=(float)linewidth;
					}
					kstack.addSlice("",temp);
					if(straighten){
						temp2[i]=profiler.getStraightened(pixels,width,height,polyroi,connected,linewidth,or_index);
					}
				}
				if(straighten){
					ImagePlus simp=new ImagePlus("Straightened",jutils.array2stack(temp2,linewidth,length));
					simp.copyScale(imp);
					simp.setOpenAsHyperStack(true);
					simp.setDimensions(pixels.length,1,1);
					CompositeImage ci=new CompositeImage(simp,CompositeImage.COMPOSITE);
					ci.copyLuts(imp);
					ci.show();
				}
			}
		}
		if(frames==1){
			Calibration cal=imp.getCalibration();
			String unit=cal.getUnit();
			float[] xvals2=new float[length];
			for(int i=0;i<length;i++){
				xvals2[i]=(float)cal.pixelWidth*(float)i;
			}
			PlotWindow4 pw=new PlotWindow4("Polyline Profile",unit,"Avg. Intensity",xvals2,(float[])kstack.getPixels(1));
			for(int i=1;i<nchannels;i++){
				pw.addPoints(xvals2,(float[])kstack.getPixels(i+1),true);
			}
			pw.draw();
		} else {
			ImagePlus imp2=new ImagePlus("Kymograph",kstack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(nchannels,slices,1);
			if(nchannels>1 && imp.isComposite()){
				CompositeImage ci=new CompositeImage(imp2,CompositeImage.COMPOSITE);
				ci.copyLuts(imp);
				ci.show();
			} else {
				imp2.show();
			}
		}
	}

}
