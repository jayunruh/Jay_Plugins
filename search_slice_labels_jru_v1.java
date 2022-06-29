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
import ij.plugin.*;
import ij.plugin.frame.*;
import jguis.*;
import jalgs.*;
import ij.text.*;
import java.util.*;

public class search_slice_labels_jru_v1 implements PlugIn {
	//this version searches for substrings, not exact matches

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		GenericDialog gd=new GenericDialog("Name List");
		gd.addCheckbox("Use_Table_Column",false);
		gd.addTextAreas("",null,10,20);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean usetable=gd.getNextBoolean();
		String input=gd.getNextText();
		if(usetable){
			TextWindow[] tw=jutils.selectTables(false,1);
			if(tw!=null && tw[0]!=null){
				TextPanel tp=tw[0].getTextPanel();
				String[] col_labels=table_tools.getcollabels(tp);
				GenericDialog gd2=new GenericDialog("options");
				gd2.addChoice("Search_Column",col_labels,col_labels[0]);
				gd2.showDialog(); if(gd2.wasCanceled()) return;
				int colindex=gd2.getNextChoiceIndex();
				List<List<String>> listtable=table_tools.table2listtable(tp);
				String[] colvals=table_tools.get_listtable_column(listtable,colindex);
				input=table_tools.print_string_array(colvals,3);
			}
		}
		//IJ.log("input = \n"+input);
		int nchan=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		if(nframes==1){
			nframes=nslices;
			nslices=1;
		}
		ImageStack stack=imp.getStack();
		String[] labels=new String[nframes];
		for(int i=0;i<nframes;i++){
			int firstframe=i*nchan*nslices;
			labels[i]=stack.getSliceLabel(firstframe+1);
			if(labels[i]==null) labels[i]=""+(i+1);
		}
		int[] foundindices=new int[stack.getSize()];
		int nfound=0;
		String[] list=(new delimit_string(' ')).getrows(input);
		for(int i=0;i<list.length;i++){
			list[i]=list[i].trim();
			//IJ.log(""+list[i]);
		}
		//for(int i=0;i<list.length;i++) list[i]=list[i].toLowerCase();
		for(int i=0;i<list.length;i++){
			boolean found=false;
			for(int j=0;j<labels.length;j++){
				if(labels[j].indexOf(list[i])>=0){
					foundindices[nfound]=j;
					nfound++;
					found=true;
					break;
				}
			}
			if(!found){
				IJ.log(list[i]+" not found");
			}
		}
		IJ.log("found "+nfound+" labels, assembling stack");
		ImageStack stack2=new ImageStack(imp.getWidth(),imp.getHeight());
		for(int i=0;i<nfound;i++){
			if(nchan>1){
				int foundslice=(int)(foundindices[i]/(nchan*nslices));
				for(int j=0;j<nchan*nslices;j++){
					int temp=j+foundslice*nchan*nslices;
					stack2.addSlice(labels[temp],stack.getPixels(temp+1));
				}
			} else {
				int foundslice=(int)(foundindices[i]/nslices);
				for(int j=0;j<nslices;j++){
					int temp=j+foundslice*nslices;
					stack2.addSlice(labels[temp],stack.getPixels(temp+1));
				}
			}
		}
		if(nchan==1) new ImagePlus("Search Results",stack2).show();
		else{
			ImagePlus outimp=new ImagePlus("Search Results",stack2);
			outimp.setOpenAsHyperStack(true);
			outimp.setDimensions(nchan,nslices,nfound);
			new CompositeImage(outimp,CompositeImage.COLOR).show();
		}
	}

}
