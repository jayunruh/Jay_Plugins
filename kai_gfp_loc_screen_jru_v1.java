/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import jalgs.*;
import jguis.*;
import jalgs.jseg.*;
import ij.text.*;
import ij.io.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.*;

public class kai_gfp_loc_screen_jru_v1 implements PlugIn {
	List<List<String>> listtable;

	public void run(String arg) {
		/*ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Puncta","Membranes","Hoechst"});
		if(imps==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Output_Images?",false);
		gd.addCheckbox("Close_Images",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean output=gd.getNextBoolean();
		boolean close=gd.getNextBoolean();
		//exec(imps[0],imps[1],output);
		measure_gfp mp=new measure_gfp(imps[0],imps[1],output,close);
		Thread thread=new Thread(mp);
		thread.start();*/
		DirectoryChooser dc=new DirectoryChooser("Choose Directory");
		String dir=dc.getDirectory();
		if(dir==null) return;
		String foldername=new File(dir).getName();
		IJ.log("folder name = "+foldername);
		String[] list=new File(dir).list();
		jsort.javasort_order(list);
		listtable=new ArrayList<List<String>>();
		ExecutorService executor=Executors.newFixedThreadPool(5);
		//channels are 1:hoechst, 2:gfp, 4:gfplong, 3:fm4-64
		for(int i=0;i<list.length;i++){
			if(list[i].endsWith("3.tif")){
				int rep=Integer.parseInt(list[i].substring(7,8));
				//if(rep>3) continue;
				int pos=list[i].length()-5;
				String gfpname=list[i].substring(0,pos);
				gfpname=gfpname+"4.tif";
				String hoechstname=list[i].substring(0,pos);
				hoechstname=hoechstname+"1.tif";
				//IJ.log(gfpname);
				//measure_puncta mp=new measure_puncta(dir+gfpname,dir+list[i],false,true);
				Runnable worker=new measure_gfp(dir+gfpname,dir+list[i],dir+hoechstname,false,true,listtable);
				//Thread thread=new Thread(mp);
				//thread.start();
				executor.execute(worker);
			}
			if(IJ.escapePressed()){IJ.error("analysis terminated"); executor.shutdownNow(); return;}
			//IJ.showProgress(i,list.length);
			if(i>0 && (i%200)==0){
				executor.shutdown(); //not sure if this is a good idea
				while(!executor.isTerminated()){;}
				IJ.showStatus("Image "+i+" of "+list.length+" complete");
				System.gc();
				executor=Executors.newFixedThreadPool(5);
				//(new WaitForUserDialog("click ok to move on")).show();
			}
		}
		IJ.log("process completed");
		(new WaitForUserDialog("Wait for threads")).show();
		table_tools.create_table(foldername+"_results",listtable,new String[]{"image","object","area","avg","stdev","puntarea","punctavg","membavg","hoechstavg","nucarea","nucavg","nucpuntarea","nucpuntavg","cellmembarea","cellmembavg","nucmembarea","nucmembavg"});
	}

}

class measure_gfp implements Runnable{
	private ImagePlus gfpimp;
	private ImagePlus membimp;
	private ImagePlus hoechstimp;
	private boolean output,close;
	private List<List<String>> listtable;

	public measure_gfp(String gfpname,String membname,String hoechstname,boolean output,boolean close){
		gfpimp=IJ.openImage(gfpname);
		membimp=IJ.openImage(membname);
		hoechstimp=IJ.openImage(hoechstname);
		this.output=output;
		this.close=close;
	}

	public measure_gfp(String gfpname,String membname,String hoechstname,boolean output,boolean close,List<List<String>> listtable){
		gfpimp=IJ.openImage(gfpname);
		membimp=IJ.openImage(membname);
		hoechstimp=IJ.openImage(hoechstname);
		this.output=output;
		this.close=close;
		this.listtable=listtable;
	}

	public measure_gfp(ImagePlus gfpimp,ImagePlus membimp,ImagePlus hoechstimp,boolean output,boolean close){
		this.gfpimp=gfpimp;
		this.membimp=membimp;
		this.hoechstimp=hoechstimp;
		this.output=output;
		this.close=close;
	}

	public measure_gfp(ImagePlus gfpimp,ImagePlus membimp,ImagePlus hoechstimp,boolean output,boolean close,List<List<String>> listtable){
		this.gfpimp=gfpimp;
		this.membimp=membimp;
		this.hoechstimp=hoechstimp;
		this.output=output;
		this.close=close;
		this.listtable=listtable;
	}

	public void run(){
		ImagePlus[] imps={gfpimp,membimp,hoechstimp};
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		//start by finding the cell boundaries from the membrane dye
		float[] membimage=(float[])imps[1].getProcessor().convertToFloat().getPixels();
		findblobs3 fb=new findblobs3(width,height);
		float[] cellobj=segment_yeast_trans.segment_image4(membimage,width,height,0.05f,"Max",20,700);
		//float[] cellobj=segment_yeast_trans.segment_image4(membimage,width,height,0.015f,"Max",25,700);
		fb.set_objects(cellobj);
		if(output) (new ImagePlus("cellmask",new ByteProcessor(width,height,fb.tobinary(cellobj,true),null))).show();
		float[] gfpimage=(float[])imps[0].getProcessor().convertToFloat().getPixels();
		float[] hoechstimage=(float[])imps[2].getProcessor().convertToFloat().getPixels();
		//now get the cell stats (avg,var,area)
		int ncells=fb.nobjects;
		int[] areas=fb.get_areas(cellobj);
		int[][] lims=fb.getallfilllimits(cellobj);
		float[] avgs=new float[ncells];
		float[] calcavgs=new float[ncells];
		float[] hoechstavgs=new float[ncells];
		float[] hoechstmins=new float[ncells];
		float[] hoechstmaxs=new float[ncells];
		float[] vars=new float[ncells];
		for(int i=0;i<ncells;i++){
			avgs[i]=fb.get_object_stats(cellobj,i+1,gfpimage,lims[i],"Avg");
			calcavgs[i]=fb.get_object_stats(cellobj,i+1,membimage,lims[i],"Avg");
			hoechstavgs[i]=fb.get_object_stats(cellobj,i+1,hoechstimage,lims[i],"Avg");
			hoechstmins[i]=fb.get_object_stats(cellobj,i+1,hoechstimage,lims[i],"Min");
			hoechstmaxs[i]=fb.get_object_stats(cellobj,i+1,hoechstimage,lims[i],"Max");
			vars[i]=fb.get_object_stats(cellobj,i+1,gfpimage,lims[i],"StDev");
		}

		//now segment the potential GFP puncta
		//start by blurring with a gaus stdev of 10
		float[] blurred=gfpimage.clone();
		jsmooth.blur2D(blurred,10.0f,width,height);
		//now divide by the blurred image
		for(int i=0;i<blurred.length;i++) blurred[i]=gfpimage[i]/blurred[i];
		//subtract a rolling ball of 4
		blurred=jutils.sub_roll_ball_back(blurred,4.0f,width,height);
		if(output){
			(new ImagePlus("gfpmask",new ByteProcessor(width,height,findblobs3.threshimage(blurred,0.5f),null))).show();
		}
		//now count the number of pixels for each object with value over 0.5
		float[] gfpareas=new float[ncells];
		float[] gfpavgs=new float[ncells];
		for(int i=0;i<width*height;i++){
			if(cellobj[i]>0.0f){
				if(blurred[i]>=0.5f){
					int id=(int)cellobj[i]-1;
					gfpareas[id]+=1.0f;
					gfpavgs[id]+=gfpimage[i];
				}
			}
		}
		for(int i=0;i<ncells;i++){
			if(gfpareas[i]>0.0f) gfpavgs[i]/=gfpareas[i];
		}

		//now segment the hoechst within each cell
		//first smooth the hoechst signal
		float[] hblurred=hoechstimage.clone();
		jsmooth.blur2D(hblurred,1.0f,width,height);
		//now find the mins and maxs for each cell
		for(int i=0;i<ncells;i++){
			hoechstmins[i]=fb.get_object_stats(cellobj,i+1,hblurred,lims[i],"Min");
			hoechstmaxs[i]=fb.get_object_stats(cellobj,i+1,hblurred,lims[i],"Max");
		}
		//next threshold halfway between the min and max
		float[] nucgfpavgs=new float[ncells];
		float[] nucgfpareas=new float[ncells];
		float[] nucpgfpavgs=new float[ncells];
		float[] nucpgfpareas=new float[ncells];
		float[] nucobj=new float[width*height];
		for(int i=0;i<width*height;i++){
			if(cellobj[i]>0.0f){
				int id=(int)cellobj[i]-1;
				float thresh=hoechstmins[id]+0.5f*(hoechstmaxs[id]-hoechstmins[id]);
				if(hblurred[i]>thresh){
					nucobj[i]=(float)(id+1);
					nucgfpareas[id]+=1.0f;
					nucgfpavgs[id]+=gfpimage[i];
					if(blurred[i]>=0.5f){
						nucpgfpareas[id]+=1.0f;
						nucpgfpavgs[id]+=gfpimage[i];
					}
				}
			}
		}
		for(int i=0;i<ncells;i++){
			if(nucgfpareas[i]>0.0f) nucgfpavgs[i]/=nucgfpareas[i];
			if(nucpgfpareas[i]>0.0f) nucpgfpavgs[i]/=nucpgfpareas[i];
		}

		//need to generate membrane and nuc membrane compartments
		fb.erodeobjects3(cellobj); //this version doesn't renumber
		//fb.erodeobjects(cellobj);
		float[] cellcirc=fb.get_circ(cellobj,2);
		int[][] cellcirclims=fb.getallfilllimits(cellcirc);
		fb.erodeobjects3(nucobj); //this version doesn't renumber
		float[] nuccirc=fb.get_circ(nucobj,2);
		int[][] nuccirclims=fb.getallfilllimits(nuccirc);
		//IJ.log(""+cellcirclims.length+" , "+nuccirclims.length+" , "+fb.nobjects);
		float[] cellcircavgs=new float[ncells];
		float[] cellcircareas=new float[ncells];
		float[] nuccircavgs=new float[ncells];
		float[] nuccircareas=new float[ncells];
		for(int i=0;i<ncells;i++){
			cellcircavgs[i]=fb.get_object_stats(cellcirc,i+1,gfpimage,cellcirclims[i],"Avg");
			cellcircareas[i]=fb.get_object_stats(cellcirc,i+1,gfpimage,cellcirclims[i],"Count");
		 	nuccircavgs[i]=fb.get_object_stats(nuccirc,i+1,gfpimage,nuccirclims[i],"Avg");
			nuccircareas[i]=fb.get_object_stats(nuccirc,i+1,gfpimage,nuccirclims[i],"Count");
		}
		
		//now output the results to the table
		for(int i=0;i<ncells;i++){
			List<String> row=new ArrayList<String>();
			row.add(imps[0].getTitle()); row.add(""+i); row.add(""+areas[i]); row.add(""+avgs[i]); row.add(""+vars[i]); row.add(""+gfpareas[i]); row.add(""+gfpavgs[i]); row.add(""+calcavgs[i]);
			row.add(""+hoechstavgs[i]); row.add(""+nucgfpareas[i]); row.add(""+nucgfpavgs[i]); row.add(""+nucpgfpareas[i]); row.add(""+nucpgfpavgs[i]);
			row.add(""+cellcircareas[i]); row.add(""+cellcircavgs[i]); row.add(""+nuccircareas[i]); row.add(""+nuccircavgs[i]);
			listtable.add(row);
		}
		//IJ.log(membimp.getTitle()+" finished");
		System.out.println(membimp.getTitle()+" finished");
		if(close){
			gfpimp.close();
			membimp.close();
			hoechstimp.close();
		}
	}

}
