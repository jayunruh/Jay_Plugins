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
import jalgs.*;
import jalgs.jseg.*;
import jguis.*;

public class cell_quadrant_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_divisions",4,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int ndiv=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		Roi roi=imp.getRoi();
		int x1=((Line)roi).x1;
		int x2=((Line)roi).x2;
		int y1=((Line)roi).y1;
		int y2=((Line)roi).y2;
		float length=(float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		float rad=length/2.0f;
		float xc=(float)Math.abs(0.5f*(double)(x2+x1));
		float yc=(float)Math.abs(0.5f*(double)(y2+y1));
		int xstart=(int)(xc-rad);
		if(xstart<0) xstart=0;
		int xend=(int)(xc+rad);
		if(xend>=width) xend=width-1;
		int ystart=(int)(yc-rad);
		if(ystart<0) ystart=0;
		int yend=(int)(yc+rad);
		if(yend>=height) yend=height-1;
		float[][] profiles=new float[ndiv][slices];
		double dtheta=2.0*Math.PI/(double)ndiv;
		double startangle=(double)measure_object.get_angle((float)x2,(float)y2,xc,yc);
		//int[] counts=new int[4];
		for(int k=0;k<slices;k++){
			float[] image=algutils.convert_arr_float2(stack.getPixels(k+1));
			for(int i=ystart;i<=yend;i++){
				float ydist2=((float)i-yc)*((float)i-yc);
				for(int j=xstart;j<=xend;j++){
					float xdist2=((float)j-xc)*((float)j-xc);
					if((xdist2+ydist2)<(rad*rad)){
						float value=image[j+i*width];
						if(value!=0.0f){
							double angle=(double)measure_object.get_angle((float)j,(float)i,xc,yc);
							angle-=startangle;
							if(angle<0) angle+=2.0*Math.PI;
							boolean found=false;
							if(angle<(0.5*dtheta) || angle>(2.0*Math.PI-0.5*dtheta)){
								profiles[0][k]+=value;
								found=true;
							}
							if(!found){
								for(int l=1;l<ndiv;l++){
									if(angle<((double)l*dtheta+0.5*dtheta)){
										profiles[l][k]+=value;
										found=true;
										break;
									}
								}
							}
						}
					}
				}
			}
		}
		if(slices==1){
			float[] profiles2=new float[ndiv];
			float[] angles=new float[ndiv];
			for(int i=0;i<ndiv;i++){
				profiles2[i]=profiles[i][0];
				angles[i]=360.0f*(float)i/(float)ndiv;
			}
			new PlotWindow4("Cell Quadrant Profile","angle","sum intensity",angles,profiles2).draw();
		} else {
			//IJ.log(""+counts[0]+" , "+counts[1]+" , "+counts[2]+" , "+counts[3]);
			new PlotWindow4("Cell Quadrant Profiles","frame","sum intensity",profiles,null).draw();
		}
	}

}
