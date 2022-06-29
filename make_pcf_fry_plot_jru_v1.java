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
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import jalgs.jsim.*;
import jalgs.*;

public class make_pcf_fry_plot_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin simulates particles of a particular size according to a pair correlation function
		//correlation function is normalized to shell volume so have to unnormalize
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float dx=xvals[0][1]-xvals[0][0];
		float maxx=xvals[0][npts[0]-1];
		double[] pdist=algutils.convert_arr_double(yvals[0].clone());
		double sumval=pdist[0];
		for(int i=1;i<pdist.length;i++){
			pdist[i]*=(double)(dx*(float)i);
			sumval+=pdist[i];
		}
		for(int i=0;i<pdist.length;i++){
			pdist[i]/=sumval;
		}
		//now make and image with random points placed according to the probability distribution
		//make the pixel size equal to dx
		rngs random=new rngs();
		int maxdist=(int)(maxx/dx)+1;
		int imgsize=2*maxdist;
		int nsims=1000;
		float pointrad=0.05f;
		int scaling=4;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Number of Points",nsims,0);
		gd.addNumericField("Point Radius",pointrad,5,15,null);
		gd.addNumericField("Image Scaling",scaling,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		nsims=(int)gd.getNextNumber();
		pointrad=(float)gd.getNextNumber();
		scaling=(int)gd.getNextNumber();
		float pointradpix=pointrad/dx;
		float[][] coords=new float[2][nsims];
		byte[] simimg=new byte[imgsize*imgsize*scaling*scaling];
		ByteProcessor bp=new ByteProcessor(imgsize*scaling,imgsize*scaling,simimg);
		//bp.setValue(255);
		bp.setColor(Color.white);
		for(int i=0;i<nsims;i++){
			double randdist=random.arbdev(pdist);
			if(randdist>=(double)maxdist) randdist=(double)maxdist;
			double randtheta=random.unidev(2.0*Math.PI,0.0);
			double xpos=randdist*Math.cos(randtheta);
			double ypos=randdist*Math.sin(randtheta);
			//coords[0][i]=(float)(xpos+(double)maxdist);
			//coords[1][i]=(float)(ypos+(double)maxdist);
			coords[0][i]=(float)xpos;
			coords[1][i]=(float)ypos;
			drawCircle(bp,imgsize*scaling,imgsize*scaling,scaling*(int)(xpos+(double)maxdist),scaling*(int)(ypos+(double)maxdist),scaling*pointradpix);
		}
		new PlotWindow4("Random PCF Simulation","x (pixels)","y",coords[0],coords[1]).draw();
		new ImagePlus("Fry Image",bp).show();
		//now generate the radial PCF image
		float[] radialimg=new float[imgsize*imgsize];
		for(int i=0;i<imgsize;i++){
			for(int j=0;j<imgsize;j++){
				float dist=(float)Math.sqrt((float)((j-maxdist)*(j-maxdist)+(i-maxdist)*(i-maxdist)));
				radialimg[j+i*imgsize]=interpolation.interp1D(yvals[0],npts[0],dist);
			}
		}
		new ImagePlus("Probability Map",new FloatProcessor(imgsize,imgsize,radialimg)).show();
	}

	public void drawCircle(ByteProcessor bp,int width,int height,int xc,int yc,float radius){
		if(radius<=0.5f){
			bp.drawPixel(xc,yc);
		} else {
			int xstart=xc-(int)radius;
			int ystart=yc+(int)radius;
			int xend=xstart+(int)(2.0f*radius);
			int yend=ystart-(int)(2.0f*radius);
			jutils.fill_ellipse(bp,xstart,ystart,xend,yend,1.0f);
		}
	}

}
