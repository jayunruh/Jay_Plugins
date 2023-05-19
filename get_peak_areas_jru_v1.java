/*******************************************************************************
 * Copyright (c) 2023 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;

public class get_peak_areas_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin takes a line plot and a series of line rois (approximating the baselines) and calculates areas
		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4SelCopy(iw);
		//float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		//float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		//int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		//int selected=(int)jutils.runPW4VoidMethod(iw,"getSelected");
		RoiManager rman=RoiManager.getInstance();
		if(rman==null || rman.getCount()<1){
			IJ.log("Please put baseline selections in Roi Manager");
			return;
		}
		//get all of the peak positions on the x axis
		Roi[] rois=rman.getRoisAsArray();
		int npeaks=rois.length;
		//these are start and end points
		int[][] coords=new int[npeaks][4];
		Plot4 plot=pw.getPlot();
		float[] limits=pw.getLimits();
		float xpixsize=(limits[1]-limits[0])/(float)Plot4.WIDTH;
		float ypixsize=(limits[3]-limits[2])/(float)Plot4.HEIGHT;
		float pixvol=xpixsize*ypixsize;
		IJ.log("Pixel Volume = "+pixvol);
		for(int i=0;i<npeaks;i++){
			//Rectangle r=rois[i].getBounds();
			Line r=(Line)rois[i];
			coords[i]=new int[]{r.x1-Plot4.LEFT_MARGIN,r.y1-Plot4.TOP_MARGIN,r.x2-Plot4.LEFT_MARGIN,r.y2-Plot4.RIGHT_MARGIN};
			float xc=0.5f*(float)(r.x1+r.x2);
			float yc=0.5f*(float)(r.y1+r.y2);
			//IJ.log("peak "+(i+1)+": "+xc+", "+yc);
			//coords[i]=plot.getPlotCoords(r.x+r.width/2,r.y+r.height/2)[0];
			//IJ.log("peak "+(i+1)+": "+coords[i]);
		}
		int[] cidxs=plot.getColors();
		cidxs[0]=3;
		int[] cpix=(int[])plot.getProcessor().getPixelsCopy();
		byte[] bwpix=new byte[Plot4.WIDTH*Plot4.HEIGHT];
		int totwidth=Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN;
		pw.close();
		for(int i=0;i<Plot4.HEIGHT;i++){
			int ypos=i+Plot4.TOP_MARGIN;
			for(int j=0;j<Plot4.WIDTH;j++){
				int xpos=j+Plot4.LEFT_MARGIN;
				int idx=xpos+ypos*totwidth;
				if(cpix[idx]!=0xFFFF0000){
					bwpix[j+i*Plot4.WIDTH]=(byte)255;
				} else {
					bwpix[j+i*Plot4.WIDTH]=(byte)0;
				}
			}
		}
		//now "rain" fill the top half of the image
		for(int j=0;j<Plot4.WIDTH;j++){
			for(int i=0;i<Plot4.HEIGHT;i++){
				if(bwpix[j+i*Plot4.WIDTH]==(byte)0){
					break;
				} else {
					bwpix[j+i*Plot4.WIDTH]=(byte)0;
				}
			}
		}
		//now "rain" fill below the roi lines
		for(int i=0;i<npeaks;i++){
			fillBelowLine(bwpix,Plot4.WIDTH,Plot4.HEIGHT,coords[i]);
		}
		//finally fill all columns that aren't occupied by a line and sum the others
		int[] areas=new int[npeaks];
		for(int i=0;i<Plot4.WIDTH;i++){
			boolean inline=false;
			for(int j=0;j<npeaks;j++){
				if(i>=coords[j][0] && i<=coords[j][2]){
					inline=true;
					for(int k=0;k<Plot4.HEIGHT;k++){
						if(bwpix[i+k*Plot4.WIDTH]!=(byte)0) areas[j]++;
					}
				}
			}
			if(!inline){
				for(int j=0;j<Plot4.HEIGHT;j++){
					bwpix[i+j*Plot4.WIDTH]=(byte)0;
				}
			}
		}
		for(int j=0;j<npeaks;j++){
			IJ.log("peak "+(j+1)+" area (pixels) = "+areas[j]);
		}
		new ImagePlus("Peak_Areas",new ByteProcessor(Plot4.WIDTH,Plot4.HEIGHT,bwpix)).show();
	}

	public void fillBelowLine(byte[] img,int width,int height,int[] coords){
		//scan the line in x and drop fill from there to the bottom
		float yinc=(float)(coords[3]-coords[1])/(float)(coords[2]-coords[0]);
		for(int i=coords[0];i<=coords[2];i++){
			int start=1+(int)(coords[1]+yinc*(i-coords[0]));
			for(int j=start;j<height;j++){
				img[i+j*width]=(byte)0;
			}
		}
		return;
	}

}
