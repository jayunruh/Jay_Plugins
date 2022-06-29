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

public class stitch_pe_library_jru_v2 implements PlugIn,gui_interface,FrameInterface {

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
		//gd.addCheckbox("Subtract Avg Back",true);
		//gd.addCheckbox("Smooth Back",false);
		//gd.addNumericField("Smooth_Stdev",100,5,15,null);
		gd.addNumericField("Start_C",0,0);
		gd.addNumericField("End_C",-1,0);
		gd.addNumericField("Start_Z",0,0);
		gd.addNumericField("End_Z",-1,0);
		gd.addNumericField("Start_T",0,0);
		gd.addNumericField("End_T",-1,0);
		gd.addCheckbox("Z_Project",false);
		gd.addChoice("Z_Proj_Stat",jstatistics.stats,jstatistics.stats[2]);
		gd.addCheckbox("Plot_Has_Units",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int selindex=gd.getNextChoiceIndex();
		//boolean sub=gd.getNextBoolean();
		//boolean smoothback=gd.getNextBoolean();
		//float sbstdev=(float)gd.getNextNumber();
		int cstart=(int)gd.getNextNumber();
		int cend=(int)gd.getNextNumber();
		int zstart=(int)gd.getNextNumber();
		int zend=(int)gd.getNextNumber();
		int tstart=(int)gd.getNextNumber();
		int tend=(int)gd.getNextNumber();
		boolean proj=gd.getNextBoolean();
		String projstat=jstatistics.stats[gd.getNextChoiceIndex()];
		boolean plotunits=gd.getNextBoolean();
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

		SaveDialog sd=new SaveDialog("Save As","",".tif");
		String outdir=sd.getDirectory();
		String outfname=sd.getFileName();
		if(outfname==null || outfname.length()==0) return;

		//now read in the huge image stack, stitching and saving as we go
		int width=0; int height=0; int nchan=0; float psize=0.0f; int stacksize=0; int nslices=0; int nframes=0;
		int[] limits={cstart,cend,zstart,zend,tstart,tend};
		int[] sizes=(new LOCI_file_reader()).get_imp_sizes(directory,fname,selseries[order[0]]);
		width=sizes[0]; height=sizes[1]; nchan=sizes[2]; nslices=sizes[3]; nframes=sizes[4];
		if(limits[0]<0) limits[0]=0;
		if(limits[1]<0 || limits[1]>=nchan) limits[1]=nchan-1;
		if(limits[2]<0) limits[2]=0;
		if(limits[3]<0 || limits[3]>=nslices) limits[3]=nslices-1;
		if(limits[4]<0) limits[4]=0;
		if(limits[5]<0 || limits[5]>=nframes) limits[5]=nframes-1;
		nslices=limits[3]-limits[2]+1;
		if(proj) nslices=1;
		nframes=limits[5]-limits[4]+1;
		nchan=limits[1]-limits[0]+1;
		LOCI_series_reader lsr=new LOCI_series_reader(directory,fname,false);
		//ImagePlus testimp=(new LOCI_file_reader()).get_loci_subimp(directory,fname,false,selseries[order[0]],false,"Max",-1,new int[]{0,-1,0,0,0,0});
		ImagePlus testimp=lsr.getSubImp(selseries[order[0]],false,"Max",new int[]{0,-1,0,0,0,0});
		psize=(float)jutils.get_psize(testimp);

		ImageWindow[] iw=jutils.selectPlots(false,1,new String[]{"coord_plot (pixels)"}); //automatically assume plot has pixel units--maybe shouldn't
		if(iw==null) return;
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
		int sel=(Integer)jutils.runPW4VoidMethod(iw[0],"getSelected"); if(sel<0) sel=0;
		float[] xcoords=xvals[sel].clone();
		float[] ycoords=yvals[sel].clone();
		if(plotunits){
			float minx=jstatistics.getstatistic("Min",xcoords,null);
			float miny=jstatistics.getstatistic("Min",ycoords,null);
			for(int i=0;i<xcoords.length;i++){
				xcoords[i]=(xcoords[i]-minx)/psize;
				ycoords[i]=(ycoords[i]-miny)/psize;
			}
			//new PlotWindow4("test_coords","x","y",xcoords,ycoords).draw();
		}
		float[] overlap=stitching.getAvgOverlap(new float[][]{xcoords,ycoords},width,height);

		//now get the flatness correction if desired--should have nchan channels but may only have one
		ImagePlus[] flatback=jutils.selectImages(true,1,new String[]{"background image"});
		float[][] backpix=null;
		if(flatback!=null && flatback[0]!=null){
			Object[] temp=jutils.stack2array(flatback[0].getStack());
			backpix=new float[nchan][];
			backpix[0]=algutils.convert_arr_float(temp[0]);
			float[] backavg=new float[backpix.length];
			backavg[0]=jstatistics.getstatistic("Avg",backpix[0],null);
			for(int i=1;i<nchan;i++){
				if(temp.length>1) backpix[i]=algutils.convert_arr_float(temp[i]);
				else backpix[i]=backpix[0].clone();
				backavg[i]=jstatistics.getstatistic("Avg",backpix[i],null);
			}
			//normalize the background so that its average intensity is 1
			for(int i=0;i<nchan;i++){
				for(int j=0;j<backpix[i].length;j++) backpix[i][j]/=backavg[i];
			}
		}

		int typeindex=algutils.get_array_type(testimp.getStack().getPixels(1));
		if(proj) typeindex=2; //float type
		stitching sclass=new stitching(width,height,xcoords,ycoords,1.0f,typeindex,this);
		//create a tiny version of our image to initialize the tiff writer
		Object[] tinystack=new Object[nchan*nframes*nslices];
		for(int i=0;i<nchan*nframes*nslices;i++){tinystack[i]=algutils.create_array(5*5,typeindex);}
		ImagePlus tinyimp=jutils.create_hyperstack("tiny image",jutils.array2stack(tinystack,5,5),testimp,nframes,nslices,nchan);
		Tiff_Writer tw=new Tiff_Writer(tinyimp,sclass.newwidth,sclass.newheight,nframes*nslices*nchan,this);
		tw.saveAsTiffStack2(outdir+outfname);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nslices;j++){
				Object[][] chanstack=new Object[nchan][ntiles];
				for(int k=0;k<ntiles;k++){
					int[] limits2={limits[0],limits[1],limits[2],limits[3],limits[4]+i,limits[4]+i};
					if(!proj){limits2[2]=limits[2]+j; limits2[3]=limits[2]+j;} //if not a projection, read a single slice
					//ImagePlus imp=(new LOCI_file_reader()).get_loci_subimp(directory,fname,false,selseries[order[k]],proj,projstat,-1,limits2);
					ImagePlus imp=lsr.getSubImp(selseries[order[k]],proj,projstat,limits2);
					Object[] tempstack=jutils.stack2array(imp.getStack());
					for(int l=0;l<nchan;l++){
						chanstack[l][k]=tempstack[l];
						if(backpix!=null){
							chanstack[l][k]=divBack(chanstack[l][k],backpix[l],0.1f);
						}
					}
				}
				for(int k=0;k<nchan;k++){
					Object stitched=sclass.stitch_frame(chanstack[k],true,overlap[0],overlap[1]);
					tw.saveFrame(stitched);
				}
				IJ.showStatus("frame "+(i+1)+ ", slice "+(j+1)+" stitched");
			}
		}
		lsr.dispose();
		tw.endWrite();
		IJ.showStatus("done stitching");
	}

	public Object divBack(Object image,float[] back,float mindiv){
		int type=algutils.get_array_type(image);
		float[] fimg=algutils.convert_arr_float2(image);
		for(int i=0;i<fimg.length;i++){if(back[i]>mindiv) fimg[i]/=back[i];}
		return algutils.convert_array2(fimg,type);
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

	public Object getNextFrame(){
		return null;
	}

	public void showMessage(String message){ IJ.showMessage(message);}

	public void showProgress(int currpos,int finalpos){ IJ.showProgress(currpos,finalpos);}

}
