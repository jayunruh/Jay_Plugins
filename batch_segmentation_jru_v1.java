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
import jguis.*;
import jalgs.jseg.*;
import jalgs.*;
import java.io.*;
import ij.io.*;

public class batch_segmentation_jru_v1 implements PlugIn {
	// this is a simple batch DAPI segmentation that uses the ImageJ default auto-threshold
	boolean overwrite;

	public void run(String arg) {
		DirectoryChooser od=new DirectoryChooser("Choose Directory");
		String dir=od.getDirectory();
		if(dir==null || dir==""){return;}
		segment_directory(dir);
		overwrite=false;
	}

	public boolean fileexists(String path){
		return (new File(path)).exists();
	}

	public void segment_directory(String dir){
		String[] list=new File(dir).list();
		int nfiles=list.length;
		for(int i=0;i<nfiles;i++){
			//if(list[i].endsWith(".tiff") && list[i].contains("rc1")){
			if(list[i].endsWith(".tif") && list[i].contains("001.tif")){
				//ImagePlus imp=(new LOCI_file_reader()).get_loci_imp(dir,list[i]);
				int dotpos=list[i].lastIndexOf(".");
				String newpath=dir+list[i].substring(0,dotpos)+"_mask.tif";
				if(!overwrite && fileexists(newpath)) continue;
				ImagePlus imp=(new Opener()).openTiff(dir,list[i]);
				if(imp!=null){
					int width=imp.getWidth(); int height=imp.getHeight();
					float[] image=(float[])imp.getStack().getProcessor(1).convertToFloat().getPixels();
					//image=jsmooth.bin2D(image,width,height,2,true);
					//width/=2;
					//height/=2;
					//byte[] mask=segment_yeast_trans.segment_image2(image,width,height,99.4f,50,1000000);
					//byte[] mask=segment_yeast_trans.segment_image(image,width,height,450.0f,25,1000000);
					byte[] mask=threshimage(image);
					//findblobs3 fb=new findblobs3(width,height);
					//float[] objects=fb.dofindblobs(mask);
					//fb.filter_clusters(objects,2,2,2);
					//mask=fb.tobinary(objects,false);
					ImagePlus masked=new ImagePlus("Masked",new ByteProcessor(width,height,mask,null));
					FileSaver fs=new FileSaver(masked);
					IJ.log(newpath);
					fs.saveAsTiff(newpath);
				}
			} else {
				File temp=new File(dir+File.separator+list[i]);
				if(temp.isDirectory()){
					segment_directory(temp.getAbsolutePath()+File.separator);
				}
			}
			IJ.showProgress(i,nfiles);
		}
	}

	public byte[] threshimage(Object pix){
		float[] minmax={jstatistics.getstatistic("Min",pix,null),jstatistics.getstatistic("Max",pix,null)};
		int[] hist=jstatistics.histogram(pix,new float[]{256.0f,minmax[0],minmax[1]});
		int index=(new AutoThresholder()).getThreshold("Default",hist);
		float thresh=(((float)index)/255.0f)*(minmax[1]-minmax[0])+minmax[0];
		return findblobs3.threshimage(pix,thresh);
	}

}
