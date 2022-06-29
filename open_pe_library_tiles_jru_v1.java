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
import java.awt.Frame;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.io.*;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class open_pe_library_tiles_jru_v1 implements PlugIn,gui_interface {
	//here we open a subset of a perkin elmer mvd2 library for stitching as a time series

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		LOCI_file_reader r=new LOCI_file_reader();
		int nseries=r.getNSeries(directory,fname);
		String[] names=r.getSeriesNames(directory,fname);
		String[] keys={"X Location","Y Location"};
		String[][] vals=r.batch_get_series_metadata_value(directory,fname,keys);
		//now we group names into tiled sets (images end with '(raw tile x--)')
		String[] parentnames=new String[names.length];
		int[] tilenum=new int[names.length];
		float[][] icoords=new float[2][names.length];
		IJ.log("series name dump");
		for(int i=0;i<names.length;i++){
			int pos=names[i].indexOf("(raw tile ");
			if(pos<0){
				parentnames[i]="";
			}else{
				if(pos==0){
					parentnames[i]="unnamed";
				} else {
					parentnames[i]=names[i].substring(0,pos-1);
				}
				int pos2=names[i].indexOf(")",pos+8);
				String temp=names[i].substring(pos+10,pos2);
				tilenum[i]=Integer.parseInt(temp);
				icoords[0][i]=Float.parseFloat(vals[0][i]);
				icoords[1][i]=Float.parseFloat(vals[1][i]);
			}
			IJ.log(parentnames[i]+" , "+tilenum[i]+" , "+icoords[0][i]+" , "+icoords[1][i]);
		}
		//now get the unique names
		Object[] temptilenames=getUniqueNames(parentnames,tilenum);
		String[] tilenames=(String[])temptilenames[0];
		String[] disptilenames=(String[])temptilenames[2];
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Select Image",disptilenames,disptilenames[0]);
		gd.addNumericField("Start_C",0,0);
		gd.addNumericField("End_C",-1,0);
		gd.addNumericField("Start_Z",0,0);
		gd.addNumericField("End_Z",-1,0);
		gd.addNumericField("Start_T",0,0);
		gd.addNumericField("End_T",0,0);
		gd.addChoice("Z_Proj_Stat",jstatistics.stats,jstatistics.stats[2]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int selindex=gd.getNextChoiceIndex();
		String impname=disptilenames[selindex];
		int cstart=(int)gd.getNextNumber();
		int cend=(int)gd.getNextNumber();
		int zstart=(int)gd.getNextNumber();
		int zend=(int)gd.getNextNumber();
		int tstart=(int)gd.getNextNumber();
		int tend=(int)gd.getNextNumber();
		String projstat=jstatistics.stats[gd.getNextChoiceIndex()];
		
		//now that we have the parent name, make our series list
		int ntiles=((int[])temptilenames[1])[selindex];
		int[] selseries=new int[ntiles];
		int[] seltilenum=new int[ntiles];
		int start=((int[])temptilenames[3])[selindex];
		for(int i=0;i<ntiles;i++){
			selseries[i]=start+i;
			seltilenum[i]=tilenum[i];
		}
		int[] order=jsort.get_javasort_order(seltilenum);
		
		float[][] coords=new float[2][ntiles];
		for(int i=0;i<ntiles;i++){
			coords[0][i]=icoords[0][selseries[order[i]]]; coords[1][i]=icoords[1][selseries[order[i]]];
		}
		new PlotWindow4("Stitching_Coords","x (um)","y (um)",coords[0],coords[1]).draw();

		//now read in the huge image stack
		Object[] hugestack=null;
		int width=0; int height=0; int nchan=0; float psize=0.0f; int stacksize=0; int nslices=0; int nframes=0;
		int[] limits={cstart,cend,zstart,zend,tstart,tend};
		IJ.showStatus("reading and projecting tiles");
		for(int i=0;i<ntiles;i++){
			ImagePlus imp=(new LOCI_file_reader()).get_loci_subimp(directory,fname,false,selseries[order[i]],true,projstat,-1,limits);
			ImageStack tstack=imp.getStack();
			if(hugestack==null){
				width=imp.getWidth(); height=imp.getHeight(); nchan=imp.getNChannels();
				psize=(float)jutils.get_psize(imp); stacksize=tstack.getSize();
				nslices=imp.getNSlices(); nframes=imp.getNFrames();
				//if(nslices==1){nslices=nframes; nframes=1;}
				//IJ.log(""+nchan+" , "+nslices+" , "+nframes);
				hugestack=new Object[stacksize*ntiles];
			}
			for(int j=0;j<nframes;j++){
				for(int k=0;k<nslices*nchan;k++){
					hugestack[j*nchan*nslices*ntiles+i*nchan*nslices+k]=tstack.getPixels(j*nslices*nchan+k+1);
				}
			}
			System.gc();
		}
		ImagePlus hugeimp=jutils.create_hyperstack(impname,jutils.array2stack(hugestack,width,height),ntiles*nframes,nslices,nchan,true,null);
		jutils.set_psize(hugeimp,psize);
		FileInfo fi=hugeimp.getOriginalFileInfo();
		if(fi==null) {
			fi=new FileInfo(); fi.width=width; fi.height=height;
		}
		fi.directory=directory;
		fi.fileName=fname;
		hugeimp.setFileInfo(fi);
		hugeimp.show();
	}

	public static String getStitchedNames(String directory,String fname){
		int nseries=(new LOCI_file_reader()).getNSeries(directory,fname);
		String[] names=(new LOCI_file_reader()).getSeriesNames(directory,fname);
		//now we group names into tiled sets (images end with '(raw tile x--)')
		String[] parentnames=new String[names.length];
		int[] tilenum=new int[names.length];
		for(int i=0;i<names.length;i++){
			int pos=names[i].indexOf("(raw tile ");
			if(pos<0){
				parentnames[i]="";
			}else{
				parentnames[i]=names[i].substring(0,pos-1);
				int pos2=names[i].indexOf(")",pos+8);
				String temp=names[i].substring(pos+10,pos2);
				tilenum[i]=Integer.parseInt(temp);
			}
			//IJ.log(parentnames[i]+" , "+tilenum[i]);
		}
		//now get the unique names
		Object[] tilenames=getUniqueNames(parentnames);
		return table_tools.print_string_array((String[])tilenames[2]);
		//return tilenames;
	}

	public int[] trunc_array(int[] arr,int len){
		int[] temp=new int[len];
		System.arraycopy(arr,0,temp,0,len);
		return temp;
	}

	public static Object[] getUniqueNames(String[] list1){
		String[] list=list1.clone();
		//sort the list
		jsort.javasort_order(list);
		List<String> celllist=new ArrayList<String>();
		//first skip all of the file names that aren't tiles (set to null)
		int pos=0;
		while(list[pos].length()==0){
			pos++;
		}
		String currcell=list[pos];
		celllist.add(currcell);
		int[] counts=new int[list1.length];
		counts[0]=1;
		for(int i=pos+1;i<list.length;i++){
			if(!list[i].equals(currcell)){
				celllist.add(list[i]);
				currcell=list[i];
			}
			counts[celllist.size()-1]++;
		}
		String[] temp=new String[celllist.size()];
		int[] tempcounts=new int[celllist.size()];
		String[] temp2=new String[celllist.size()];
		for(int i=0;i<celllist.size();i++){
			temp[i]=celllist.get(i);
			tempcounts[i]=counts[i];
			temp2[i]=temp[i]+" ("+tempcounts[i]+" tiles)";
		}
		return new Object[]{temp,tempcounts,temp2};
	}

	public static Object[] getUniqueNames(String[] list1,int[] id){
		//find the unique names: assume an ordered list
		//if multiple names with the id, 1, repeat them
		//non-tile names are null
		//start by finding the first non-null name
		int pos=0;
		while(list1[pos].length()==0 && pos<list1.length){
			pos++;
		}
		List<String> namelist=new ArrayList<String>();
		int[] counts=new int[list1.length];
		int[] starts=new int[list1.length];
		String currname=list1[pos];
		namelist.add(currname);
		int nnames=0;
		counts[nnames]=1;
		starts[nnames]=pos;
		for(int i=pos+1;i<list1.length;i++){
			if(list1[i].equals(currname) && id[i]!=1){
				counts[nnames]++;
			} else {
				nnames++;
				currname=list1[i];
				namelist.add(currname);
				counts[nnames]=1;
				starts[nnames]=i;
			}
		}
		nnames++;
		String[] temp=new String[nnames];
		int[] tempcounts=new int[nnames];
		int[] tempstarts=new int[nnames];
		String[] temp2=new String[nnames];
		for(int i=0;i<nnames;i++){
			temp[i]=namelist.get(i);
			tempcounts[i]=counts[i];
			temp2[i]=temp[i]+" ("+tempcounts[i]+" tiles)";
			tempstarts[i]=starts[i];
		}
		return new Object[]{temp,tempcounts,temp2,tempstarts};
	}

	public void showMessage(String message){ IJ.showMessage(message);}

	public void showProgress(int currpos,int finalpos){ IJ.showProgress(currpos,finalpos);}

}
