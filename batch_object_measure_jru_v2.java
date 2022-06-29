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
import ij.plugin.filter.BackgroundSubtracter;

public class batch_object_measure_jru_v2 implements PlugIn {
	//this version is for batch nuclear vs. cytoplasmic intensity screening
	//I used it for Parama's HDAC5 screen
	String header;
	//float back;
	int backsubindex;
	float ballrad;
	boolean overwrite,output;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] backsubs={"None","Min","Rolling_Ball"};
		gd.addChoice("Back Sub Options",backsubs,backsubs[1]);
		gd.addNumericField("Rolling_Ball_Radius(pixels)",50.0f,5,15,null);
		//gd.addCheckbox("Overwrite Files",false);
		gd.addCheckbox("Output_Overlay",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		backsubindex=gd.getNextChoiceIndex();
		ballrad=(float)gd.getNextNumber();
		//overwrite=gd.getNextBoolean();
		output=gd.getNextBoolean();
		header="file\tid\tnucavg\tcytoavg\tratio\tback\tratio2\trow\tcol\twellid";
		//back=0.0f;
		DirectoryChooser od=new DirectoryChooser("Choose Directory");
		String dir=od.getDirectory();
		if(dir==null) return;
		measure_directory(dir);
		IJ.log("analysis finished");
	}

	public void measure_directory(String dir){
		String[] list=new File(dir).list();
		int nfiles=list.length;
		StringBuffer sb=new StringBuffer();
		sb.append(header+"\n");
		int entries=0;
		for(int i=0;i<nfiles;i++){
			if(list[i].endsWith("001.tif")){
				//at some point we should try to read in the Summary.xls if it exists and only analyze new files
				//for now just reanalyze everything
				String beginning=list[i].substring(0,17);
				//String ending=list[i].substring(15,28);
				String gfpname=beginning+"3.tif";
				if(new File(dir+gfpname).exists()){
					ImagePlus imp=(new Opener()).openTiff(dir,gfpname);
					if(imp==null) continue;
					int width=imp.getWidth(); int height=imp.getHeight();
					float[] image=algutils.convert_arr_float(imp.getProcessor().getPixels());
					//optionally subtract a smoothed background
					float back=0.0f;
					if(backsubindex==1) back=jstatistics.getstatistic("not0min",image,null); //need to avoid dark pixels
					if(backsubindex==2){
						//float[] blurred=image.clone();
						//jsmooth.blur2D(blurred,50.0f,width,height);
						//for(int j=0;j<image.length;j++) image[j]-=blurred[j];
						FloatProcessor fp2=new FloatProcessor(width,height,image,null);
						fp2.snapshot();
						BackgroundSubtracter bs=new BackgroundSubtracter();
						bs.rollingBallBackground(fp2,100.0,false,false,false,true,true);
					}
					ImagePlus maskimp=(new Opener()).openTiff(dir,list[i]);
					if(maskimp==null) continue;
					float[] dapiimage=algutils.convert_arr_float(maskimp.getProcessor().getPixels());
					//byte[] mask=(byte[])maskimp.getProcessor().getPixels();
					byte[] mask=threshimage(dapiimage);
					findblobs3 fb=new findblobs3(width,height);
					float[] objects=fb.dofindblobs(mask);
					fb.clear_edges(objects,true);
					fb.filter_area(objects,new int[]{4,1000000},true);
					int nobjects=fb.nobjects;
					int[][] objlims=fb.getallfilllimits(objects);
					float[] circ=fb.get_circ(objects,4);
					int[][] circlims=fb.getallfilllimits(circ);
					for(int k=0;k<nobjects;k++){
						float nucavg=fb.get_object_stats(objects,k+1,image,objlims[k],"Avg");
						float cytoavg=fb.get_object_stats(circ,k+1,image,circlims[k],"Avg");
						float ratio=nucavg/cytoavg;
						float ratio2=(nucavg-back)/(cytoavg-back);
						String row=beginning.substring(0,3);
						String col=beginning.substring(3,6);
						String wellid=beginning.substring(0,6);
						sb.append(""+beginning+"\t"+(k+1)+"\t"+nucavg+"\t"+cytoavg+"\t"+ratio+"\t"+back+"\t"+ratio2+"\t"+row+"\t"+col+"\t"+wellid+"\n");
					}
					//now optionally output a measurement image that has the circ as a second channel
					if(output){
						ImageStack dispstack=new ImageStack(width,height);
						dispstack.addSlice("",image);
						dispstack.addSlice("",circ);
						ImagePlus dispimp=new ImagePlus("Analysis",dispstack);
						dispimp.setOpenAsHyperStack(true);
						dispimp.setDimensions(2,1,1);
						//dispimp.show();
						//FileSaver fs=new FileSaver(new CompositeImage(dispimp,CompositeImage.COMPOSITE));
						FileSaver fs=new FileSaver(dispimp);
						String newpath=dir+beginning+"3_overlay.tif";
						IJ.log(newpath);
						fs.saveAsTiffStack(newpath);
					} else {
						IJ.log(list[i]);
					}
					entries++;
					//return;
				}
			} else {
				File temp=new File(dir+list[i]);
				if(temp.isDirectory()){
					measure_directory(temp.getAbsolutePath()+File.separator);
				}
			}

			IJ.showProgress(i,nfiles);
		}
		if(entries>0){
			//here we need to write the string to file
			try{
				File summary=new File(dir+"Summary.xls");
				BufferedWriter d=new BufferedWriter(new FileWriter(summary));
				d.write(sb.toString());
				d.close();
			} catch(IOException e){
				IJ.showMessage("error writing file");
			}
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
