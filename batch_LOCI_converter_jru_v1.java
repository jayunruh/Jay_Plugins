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
import ij.io.*;
import jguis.*;
import jalgs.*;

public class batch_LOCI_converter_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		int nseries=(new LOCI_file_reader()).getNSeries(directory,fname);
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Avoid_Metadata",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		boolean nometa=gd2.getNextBoolean();
		LOCI_file_reader lfr=new LOCI_file_reader();
		lfr.nometa=nometa;
		String[] names=lfr.getSeriesNames(directory,fname);
		GenericDialog gd=new GenericDialog("Options");
		int minsize=5;
		int maxsize=1000000;
		int binz=1;
		gd.addNumericField("Min_stack_size (all frames)",minsize,0);
		gd.addNumericField("Max_stack_size (all frames)",maxsize,0);
		gd.addNumericField("Z_bin_size",binz,0);
		gd.addCheckbox("Max_Z_Project?",true);
		gd.addCheckbox("Different_Output_Directory?",false);
		gd.addCheckbox("Convert_All_Channels?",true);
		gd.addNumericField("Channel_To_Convert",1,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		minsize=(int)gd.getNextNumber();
		maxsize=(int)gd.getNextNumber();
		binz=(int)gd.getNextNumber();
		boolean zproj=gd.getNextBoolean();
		boolean diffout=gd.getNextBoolean();
		boolean allchan=gd.getNextBoolean();
		int convchan=(int)gd.getNextNumber()-1;
		String outdir=directory;
		if(diffout){
			DirectoryChooser dc=new DirectoryChooser("Save Directory");
			outdir=dc.getDirectory();
			if(outdir==null) return;
		}
		for(int i=0;i<nseries;i++){
			int[] sizes=lfr.get_imp_sizes(directory,fname,i);
			int totsize=sizes[2]*sizes[3]*sizes[4];
			IJ.log("series = "+i+", size = "+totsize);
			if(totsize>minsize && totsize<=maxsize){
				ImagePlus imp=lfr.get_loci_imp(directory,fname,false,i);
				if(imp!=null){
					ImageStack stack=imp.getStack();
					int width=imp.getWidth(); int height=imp.getHeight();
					int frames=imp.getNFrames();
					int slices=imp.getNSlices();
					int channels=imp.getNChannels();
					ImageStack binstack=new ImageStack(width,height);
					int newslices=0;
					for(int f=0;f<frames;f++){
						if(allchan){
							Object[][] binned=new Object[channels][];
							for(int j=0;j<channels;j++){
								Object[] tempstack=jutils.get3DZSeries(stack,j,f,frames,slices,channels);
								if(zproj){
									binned[j]=new Object[]{algutils.get_stack_proj_stat("Max",tempstack,width,height,tempstack.length,null)};
								} else {
									if(binz<=1) binned[j]=tempstack;
									else binned[j]=algutils.bin_stack(tempstack,width,height,binz);
								}
							}
							newslices=binned[0].length;
							for(int j=0;j<newslices;j++){
								for(int k=0;k<channels;k++){
									binstack.addSlice("",binned[k][j]);
								}
							}
						} else {
							Object[][] binned=new Object[1][];
							Object[] tempstack=jutils.get3DZSeries(stack,convchan,f,frames,slices,channels);
							if(zproj){
								binned[0]=new Object[]{algutils.get_stack_proj_stat("Max",tempstack,width,height,tempstack.length,null)};
							} else {
								if(binz<=1) binned[0]=tempstack;
								else binned[0]=algutils.bin_stack(tempstack,width,height,binz);
							}
							newslices=binned[0].length;
							for(int j=0;j<newslices;j++){
								binstack.addSlice("",binned[0][j]);
							}
						}
					}
					int[] colors=new int[1];
					if(allchan) colors=new int[channels];
					colors[0]=2;
					if(colors.length>1) colors[1]=3;
					if(colors.length>2) colors[2]=0;
					for(int c=3;c<colors.length;c++) colors[c]=1;
					LUT[] luts=new LUT[colors.length];
					for(int c=0;c<colors.length;c++) luts[c]=jutils.get_lut_for_color(jutils.colors[colors[c]]);
					if(binstack.getSize()>1){
						CompositeImage imp2=new CompositeImage(new ImagePlus("temp",binstack),CompositeImage.COMPOSITE);
						imp2.setDimensions(colors.length,newslices,frames);
						imp2.setLuts(luts);
						imp2.copyScale(imp);
						FileSaver fs=new FileSaver(imp2);
						String outname=names[i];
						outname=outname.replace("/","_");
						outname=outname.replace("\\","_");
						fs.saveAsTiffStack(outdir+names[i]+".tif");
						imp2.close();
						IJ.showStatus(names[i]+" exported");
					}
				}
				imp.close();
			}
			if(IJ.escapePressed()){break;}
		}
	}

}
