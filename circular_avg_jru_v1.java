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
	public float twopi;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean mirror=false;
		gd.addCheckbox("Show Mirror Image?",mirror);
		gd.addCheckbox("Use_Center_of_Mass",false);
		gd.showDialog();  if(gd.wasCanceled()){return;}
		mirror=gd.getNextBoolean();
		boolean com=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		float xc=(float)(width/2);
		float yc=(float)(height/2);
		Roi roi=imp.getRoi();
		if(roi!=null){
			Rectangle r=roi.getBounds();
			xc=(float)r.x+0.5f*(float)r.width;
			yc=(float)r.y+0.5f*(float)r.height;
		} else if(com){
			float[] com2=interpolation.center_of_mass_2D(stack.getPixels(1),width,height);
			xc=com2[0]; yc=com2[1];
			IJ.log("Center of Mass = "+xc+" , "+yc);
		}
		Object retvals=exec(imp,xc,yc,mirror);
		if(retvals instanceof Plot4) new PlotWindow4("Circular Avg",(Plot4)retvals).draw();
		else ((ImagePlus)retvals).show();
	}

	public Object exec(ImagePlus imp,float xc,float yc,boolean mirror){
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		//int slices=stack.getSize();
		int nslices=imp.getNSlices();
		int nchan=imp.getNChannels();
		int nframes=imp.getNFrames();
		int rsize=get_closest_edge_dist(xc,yc,width,height);
		//IJ.log(""+rsize);
		twopi=2.0f*(float)Math.PI;
		float[][] avgvals=new float[nchan*nframes][];
		float psize=(float)jutils.get_psize(imp);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nchan;j++){
				Object[] data=jutils.get3DZSeries(stack,j,i,nframes,nslices,nchan);
				avgvals[j+i*nchan]=analyzeStack(data,width,height,rsize,xc,yc);
				if(mirror) avgvals[j+i*nchan]=makeMirror(avgvals[j+i*nchan],rsize,nslices);
			}
		}
		if(nslices==1){
			float[] xvals=new float[avgvals[0].length];
			if(!mirror){
				for(int i=0;i<avgvals[0].length;i++){xvals[i]=psize*(float)i;}
			} else {
				for(int i=0;i<avgvals[0].length;i++){xvals[i]=psize*(float)(i-(rsize-1));}
			}
			float[][] xvals2=new float[nchan*nframes][];
			for(int i=0;i<nchan*nframes;i++) xvals2[i]=xvals;
			//IJ.log(""+xvals2[0].length);
			return new Plot4("r","Intensity",xvals2,avgvals,null);
		} else {
			ImagePlus imp5=null;
			if(!mirror) imp5=new ImagePlus("Circular Avg",jutils.array2stack(avgvals,rsize,nslices));
			else imp5=new ImagePlus("Circular Avg",jutils.array2stack(avgvals,2*rsize-1,nslices));
			imp5.copyScale(imp);
			imp5.setOpenAsHyperStack(true);	
			imp5.setDimensions(nchan,1,nframes);
			return imp5;
		}
	}

	public float[] makeMirror(float[] avgvals,int rsize,int slices){
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
		return newavgvals;
	}

	public float[] analyzeStack(Object[] data,int width,int height,int rsize,float xc,float yc){
		float[] avgvals=new float[data.length*rsize];
		for(int i=0;i<data.length;i++){
			float[] profile=circavg(data[i],width,height,rsize,xc,yc);
			System.arraycopy(profile,0,avgvals,i*rsize,rsize);
		}
		return avgvals;
	}

	public float[] circavg(Object data,int width,int height,int rsize,float xc,float yc){
		float[] avgvals=new float[rsize];
		avgvals[0]=interpolation.interp2D(data,width,height,xc,yc);
		for(int j=1;j<rsize;j++){
			float angleincrement=1.0f/(float)j;
			int angles=(int)(twopi/angleincrement);
			angleincrement=twopi/(float)angles;
			for(int k=0;k<angles;k++){
				float angle=angleincrement*(float)k;
				float x=(float)Math.cos(angle)*(float)j+xc;
				float y=(float)Math.sin(angle)*(float)j+yc;
				avgvals[j]+=interpolation.interp2D(data,width,height,x,y);
			}
			avgvals[j]/=(float)angles;
		}
		return avgvals;
	}

	public int get_closest_edge_dist(float x,float y,int width,int height){
		int dist=(int)((float)(width-1)-x);
		int dist2=(int)x; if(dist2<dist) dist=dist2;
		dist2=(int)y; if(dist2<dist) dist=dist2;
		dist2=(int)((float)(height-1)-y); if(dist2<dist) dist=dist2;
		return dist;
	}

}
