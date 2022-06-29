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

public class carpet_sub_avg_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("subtract_temporal_avg",false);
		boolean maintainavg=true;
		gd.addCheckbox("maintain_average",maintainavg);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		boolean temporal=gd.getNextBoolean();
		maintainavg=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		ImageStack retstack=new ImageStack(width,height);
		float[][] colavg=null;
		if(!temporal){
			colavg=new float[slices][width];
			for(int k=0;k<slices;k++){
				float[] pix=(float[])stack.getProcessor(k+1).convertToFloat().getPixels();
				float[] result=new float[width*height];
				for(int i=0;i<width;i++){
					for(int j=0;j<height;j++){
						colavg[k][i]+=pix[j*width+i]/(float)height;
					}
				}
				if(maintainavg){
					float[] rowavg=new float[height];
					for(int i=0;i<width;i++){
						for(int j=0;j<height;j++){
							rowavg[j]+=pix[j*width+i]/(float)width;
						}
					}
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							result[i*width+j]=pix[i*width+j]-colavg[k][j]+rowavg[i];
						}
					}
				} else {
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							result[i*width+j]=pix[i*width+j]-colavg[k][j];
						}
					}
				}
				retstack.addSlice("",result);
			}
		} else {
			colavg=new float[slices][height];
			for(int k=0;k<slices;k++){
				float[] pix=(float[])stack.getProcessor(k+1).convertToFloat().getPixels();
				float[] result=new float[width*height];
				for(int i=0;i<width;i++){
					for(int j=0;j<height;j++){
						colavg[k][j]+=pix[j*width+i]/(float)width;
					}
				}
				if(maintainavg){
					float[] rowavg=new float[width];
					for(int i=0;i<width;i++){
						for(int j=0;j<height;j++){
							rowavg[i]+=pix[j*width+i]/(float)height;
						}
					}
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							result[i*width+j]=pix[i*width+j]-colavg[k][i]+rowavg[j];
						}
					}
				} else {
					for(int i=0;i<height;i++){
						for(int j=0;j<width;j++){
							result[i*width+j]=pix[i*width+j]-colavg[k][i];
						}
					}
				}
				retstack.addSlice("",result);
			}
		}
		ImagePlus imp2=new ImagePlus("Subtracted",retstack);
		imp2.copyScale(imp);
		if(slices>1){
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(slices,1,1);
			(new CompositeImage(imp2,CompositeImage.COMPOSITE)).show();
		} else {
			imp2.show();
		}
		(new PlotWindow4("Avg Profiles","position","Intensity",colavg,null)).draw();
	}

}
