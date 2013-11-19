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

public class circular_avg_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean mirror=false;
		gd.addCheckbox("Show Mirror Image?",mirror);
		gd.showDialog();  if(gd.wasCanceled()){return;}
		mirror=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int rsize=width;
		if(height<rsize){rsize=height;}
		rsize/=2;
		float twopi=2.0f*(float)Math.PI;
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		float[] avgvals=new float[rsize*slices];
		float psize=(float)jutils.get_psize(imp);
		float xc=(float)(width/2);
		float yc=(float)(height/2);
		for(int i=0;i<slices;i++){
			float[] data=(float[])stack.getPixels(i+1);
			avgvals[i*rsize]=data[width/2+width*height/2];
			for(int j=1;j<rsize;j++){
				float angleincrement=1.0f/(float)j;
				int angles=(int)(twopi/angleincrement);
				angleincrement=twopi/(float)angles;
				for(int k=0;k<angles;k++){
					float angle=angleincrement*(float)k;
					float x=(float)Math.cos(angle)*(float)j+xc;
					float y=(float)Math.sin(angle)*(float)j+yc;
					avgvals[j+i*rsize]+=interpolation.interp2D(data,width,height,x,y);
				}
				avgvals[j+i*rsize]/=(float)angles;
			}
		}
		if(mirror){
			float[] newavgvals=new float[(2*rsize-1)*slices];
			int halfsize=rsize-1;
			int newwidth=2*rsize-1;
			//here make sure we don't replicate the center point twice
			for(int i=0;i<slices;i++){
				newavgvals[halfsize+i*newwidth]=avgvals[i*rsize];
				for(int j=1;j<rsize;j++){
					newavgvals[j+halfsize+i*newwidth]=avgvals[j+i*rsize];
					newavgvals[halfsize-j+i*newwidth]=avgvals[j+i*rsize];
				}
			}
			avgvals=newavgvals;
		}
		if(slices==1){
			float[] xvals=new float[avgvals.length];
			if(!mirror){
				for(int i=0;i<avgvals.length;i++){xvals[i]=psize*(float)i;}
			} else {
				for(int i=0;i<avgvals.length;i++){xvals[i]=psize*(float)(i-(rsize-1));}
			}
			new PlotWindow4("Circular Avg","r","Intensity",xvals,avgvals).draw();
		} else {
			ImagePlus imp5=null;
			if(!mirror) imp5=new ImagePlus("Circular Avg",new FloatProcessor(rsize,slices,avgvals,null));
			else imp5=new ImagePlus("Circular Avg",new FloatProcessor(2*rsize-1,slices,avgvals,null));
			imp5.copyScale(imp);
			imp5.show();
		}
	}

}
