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

public class carpet_sub_mov_avg_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean addspatialconst,addtemporalconst;
		int period,i,j,k;
		GenericDialog gd = new GenericDialog("subtract moving average");
		String[] sitems=new String[3];
		sitems[0]="No Const Add";
		sitems[1]="Maintain Spatial Average";
		sitems[2]="Maintain Temporal Average";
		gd.addChoice("Avg Correction",sitems,sitems[0]);
		gd.addNumericField("period (must be odd)",5,0);
		gd.showDialog();
		if(gd.wasCanceled()) return;
		int index=gd.getNextChoiceIndex();
		addspatialconst=false;
		addtemporalconst=false;
		if(index==1){addspatialconst=true;}
		if(index==2){addtemporalconst=true;}
		period = (int)gd.getNextNumber();
		ImagePlus imp = WindowManager.getCurrentImage();
		int width = imp.getWidth();
		int height = imp.getHeight();
		int periodshift = (period-1)/2;
		int startslice = periodshift;
		int endslice = height-startslice-1;
		ImageStack stack=imp.getStack();
		int frames=stack.getSize();
		ImageStack resstack=new ImageStack(width,endslice-startslice+1);
		ImageStack movavgstack=new ImageStack(width,endslice-startslice+1);
		for(int s=0;s<frames;s++){
			float[] pixels=(float[])stack.getProcessor(s+1).convertToFloat().getPixels();
			float[] oldcurravg = new float[width];
			float[] totavg=new float[width];
			float[] result=new float[width*(endslice-startslice+1)];
			float[] movavg=new float[width*(endslice-startslice+1)];
			for(i=startslice;i<=endslice;i++){
				//calculate the moving average
				float[] newcurravg = new float[width];
				if(i==startslice){
					//at the beginning calculate it directly
					for(j=i-periodshift;j<=i+periodshift;j++){
						for(k=0;k<width;k++){
							newcurravg[k]+=pixels[k+j*width]/((float)period);
						}
					}
					if(addtemporalconst){
						for(k=0;k<width;k++){
							totavg[k]=(newcurravg[k]*(float)period)/(float)height;
						}
					}
				} else{
					//later, just subtract the first image from the average and add a new image
					for(k=0;k<width;k++){newcurravg[k]=oldcurravg[k]-pixels[k+(i-periodshift-1)*width]/((float)period)+pixels[k+(i+periodshift)*width]/((float)period);}
					if(addtemporalconst){
						for(k=0;k<width;k++){
							totavg[k]+=pixels[k+(i+periodshift)*width]/(float)height;
						}
					}
				}
				//now subtract it from the current image
				for(k=0;k<width;k++){result[(i-startslice)*width+k]=pixels[i*width+k]-newcurravg[k];}
				//add the average intensity back in if called for
				if(addspatialconst){
					float avgint=0.0f;
					for(k=0;k<width;k++){avgint+=pixels[i*width+k]/(float)width;}
					for(k=0;k<width;k++){result[(i-startslice)*width+k]+=avgint;}
				}
				for(k=0;k<width;k++){movavg[(i-startslice)*width+k]=newcurravg[k];}
				for(k=0;k<width;k++){oldcurravg[k]=newcurravg[k];}
				IJ.showProgress(i-startslice,endslice-startslice);
			}
			if(addtemporalconst){
				int outslices=(endslice-startslice+1);
				for(i=0;i<outslices;i++){
					for(j=0;j<width;j++){
						result[j+i*width]+=totavg[j];
					}
				}
			}
			resstack.addSlice("",result);
			movavgstack.addSlice("",movavg);
		}
		//output the results
		ImagePlus resimp=new ImagePlus("Subtracted Carpet",resstack);
		resimp.copyScale(imp);
		resimp.setOpenAsHyperStack(true);
		resimp.setDimensions(frames,1,1);
		new CompositeImage(resimp,CompositeImage.COMPOSITE).show();
		ImagePlus movavgimp=new ImagePlus("Mov Avg Carpet",movavgstack);
		movavgimp.copyScale(imp);
		movavgimp.setOpenAsHyperStack(true);
		movavgimp.setDimensions(frames,1,1);
		new CompositeImage(movavgimp,CompositeImage.COMPOSITE).show();
	}

}
