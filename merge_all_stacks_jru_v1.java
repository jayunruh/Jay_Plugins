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
import jalgs.*;
import jguis.*;

public class merge_all_stacks_jru_v1 implements PlugIn {

	public void run(String arg) {
		//assume all images have same number of channels and slices
		int[] wList=WindowManager.getIDList();
		ImagePlus imp=WindowManager.getImage(wList[0]);
		ImageStack stack=imp.getStack();
		int tchannels=imp.getNChannels();
		//int tslices=imp.getNSlices();
		int totsize=0;
		int[] widths=new int[wList.length];
		int[] heights=new int[wList.length];
		//int[] channels=new int[wList.length];
		int[] slices=new int[wList.length];
		int maxwidth=0;
		int maxheight=0;
		float psize=(float)jutils.get_psize(imp);
		//int maxchans=0;
		int maxslices=0;
		for(int i=0;i<wList.length;i++){
			ImagePlus imp2=WindowManager.getImage(wList[i]);
			if(imp2!=null){
				widths[i]=imp2.getWidth();
				heights[i]=imp2.getHeight();
				//channels[i]=imp2.getNChannels();
				slices[i]=imp2.getNSlices();
				if(widths[i]>maxwidth) maxwidth=widths[i];
				if(heights[i]>maxheight) maxheight=heights[i];
				//if(channels[i]>maxchans) maxchans=channels[i];
				if(slices[i]>maxslices) maxslices=slices[i];
			}
		}
		ImageStack newstack=new ImageStack(maxwidth,maxheight);
		for(int i=0;i<wList.length;i++){
			ImagePlus imp2=WindowManager.getImage(wList[i]);
			if(imp2!=null){
				ImageStack stack2=imp2.getStack();
				//int size=stack2.getSize();
				int frames=imp2.getNFrames();
				int channels=imp2.getNChannels();
				float xshift=0.5f*(float)(maxwidth-widths[i]);
				float yshift=0.5f*(float)(maxheight-heights[i]);
				float zshift=0.5f*(float)(maxslices-slices[i]);
				//IJ.log(""+xshift+" , "+yshift+" , "+zshift);
				int size=frames*channels*maxslices;
				totsize+=size;
				int type=algutils.get_array_type(stack2.getPixels(1));
				for(int j=0;j<frames;j++){
					Object[] copied2=new Object[channels*maxslices];
					for(int k=0;k<channels;k++){
						Object[] pix=jutils.get3DZSeries(stack2,k,j,frames,slices[i],channels);
						for(int l=0;l<maxslices;l++){
							//Object copied=algutils.create_array(maxwidth*maxheight,type);
							float[] copied=new float[maxwidth*maxheight];
							for(int m=0;m<maxheight;m++){
								for(int n=0;n<maxwidth;n++){
									copied[n+m*maxwidth]=interpolation.interp3D(pix,widths[i],heights[i],n-xshift,m-yshift,l-zshift);
								}
							}
							copied2[k+l*channels]=copied;
						}
						//interpolation.shift_copy_image(pixels,widths[i],heights[i],copied,maxwidth,maxheight,xshift,yshift);
						//stack2.deleteSlice(1);
					}
					for(int k=0;k<copied2.length;k++) newstack.addSlice(imp2.getTitle(),algutils.convert_array2(copied2[k],type));
				}
				imp2.changes=false;
				imp2.close();
			}
		}
		int framesize=tchannels*maxslices;
		int nframes=(int)(totsize/framesize);
		ImagePlus imp3=new ImagePlus("Merged Stacks",newstack);
		imp3.setOpenAsHyperStack(true);
		imp3.setDimensions(tchannels,maxslices,nframes);
		jutils.set_psize(imp3,(double)psize);
		if(tchannels>1) new CompositeImage(imp3,CompositeImage.COLOR).show();
		else imp3.show();
	}

}
