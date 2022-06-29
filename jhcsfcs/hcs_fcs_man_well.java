/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

import jalgs.jsort;

import java.io.File;

public class hcs_fcs_man_well extends hcs_fcs_well{

	public hcs_fcs_man_well(String well,hcs_fcs_plate parent1){
		super(null,null,well,0,parent1);
		// get the list of files for this well
		String[] names=(new File(parent.dir+File.separator+foldername)).list();
		int numfiles=0;
		String[] rawnames=null;
		if(names!=null){
			// figure out what file names we have
			rawnames=new String[names.length];
			for(int i=0;i<names.length;i++){
				if(names[i].indexOf(".raw")>0&&names[i].indexOf("Marker")<0){
					numfiles++;
					rawnames[i]=(new parse_filenames()).get_raw_fcs_name(names[i]);
				}else{
					rawnames[i]=null;
				}
			}
		}
		if(numfiles>0){
			// weed out the names we don't need
			String[] names2=new String[numfiles];
			String[] cellnames=new String[numfiles];
			int[] cellfilenum=new int[numfiles];
			float[] cellid=new float[numfiles];
			int counter=0;
			int ncells=0;
			for(int i=0;i<names.length;i++){
				if(rawnames[i]!=null){
					names2[counter]=names[i].substring(0);
					cellid[counter]=-1.0f;
					String tempcell=rawnames[i].substring(0);
					for(int j=0;j<ncells;j++){
						if(cellnames[j].equals(tempcell)){
							cellid[counter]=j;
							cellfilenum[j]++;
							break;
						}
					}
					if(cellid[counter]<0.0f){
						cellnames[ncells]=tempcell.substring(0);
						cellfilenum[ncells]=1;
						cellid[counter]=ncells;
						ncells++;
					}
					counter++;
				}
			}
			// now we need to sort the file names by cell number
			(new jsort()).sort_string_by(names2,cellid);
			// find where each cell starts
			int[] cellstartindex=new int[ncells];
			cellstartindex[0]=0;
			counter=0;
			for(int i=1;i<numfiles;i++){
				if(cellid[i]>cellid[cellstartindex[counter]]){
					counter++;
					cellstartindex[counter]=i;
				}
			}
			// finally, initialize a cell object for each set of cell file names
			cells=new hcs_fcs_man_cell[ncells];
			for(int i=0;i<(ncells-1);i++){
				String[] cellfiles=new String[cellfilenum[i]];
				for(int j=cellstartindex[i];j<cellstartindex[i+1];j++){
					cellfiles[j-cellstartindex[i]]=names2[j].substring(0);
				}
				cells[i]=new hcs_fcs_man_cell(cellfiles,this,0,0);
			}
			String[] cellfiles=new String[cellfilenum[ncells-1]];
			for(int j=cellstartindex[ncells-1];j<numfiles;j++){
				cellfiles[j-cellstartindex[ncells-1]]=names2[j].substring(0);
			}
			cells[ncells-1]=new hcs_fcs_man_cell(cellfiles,this,0,0);
		}else{
			cells=null;
		}
	}

	public hcs_fcs_man_well(String[] filelist,String well,hcs_fcs_plate parent1){
		super(null,null,well,0,parent1);
		// get the list of files for this well
		int numfiles=filelist.length;
		String[] names=new String[numfiles];
		for(int i=0;i<numfiles;i++){
			File temp=new File(filelist[i]);
			names[i]=temp.getName();
		}
		String[] rawnames=new String[numfiles];
		for(int i=0;i<numfiles;i++){
			rawnames[i]=(new parse_filenames()).get_raw_fcs_name(names[i]);
		}
		// weed out the names we don't need
		String[] names2=new String[numfiles];
		String[] cellnames=new String[numfiles];
		int[] cellfilenum=new int[numfiles];
		float[] cellid=new float[numfiles];
		int counter=0;
		int ncells=0;
		for(int i=0;i<names.length;i++){
			if(rawnames[i]!=null){
				names2[counter]=names[i].substring(0);
				cellid[counter]=-1.0f;
				String tempcell=rawnames[i].substring(0);
				for(int j=0;j<ncells;j++){
					if(cellnames[j].equals(tempcell)){
						cellid[counter]=j;
						cellfilenum[j]++;
						break;
					}
				}
				if(cellid[counter]<0.0f){
					cellnames[ncells]=tempcell.substring(0);
					cellfilenum[ncells]=1;
					cellid[counter]=ncells;
					ncells++;
				}
				counter++;
			}
		}
		// now we need to sort the file names by cell number
		(new jsort()).sort_string_by(names2,cellid);
		// find where each cell starts
		int[] cellstartindex=new int[ncells];
		cellstartindex[0]=0;
		counter=0;
		for(int i=1;i<numfiles;i++){
			if(cellid[i]>cellid[cellstartindex[counter]]){
				counter++;
				cellstartindex[counter]=i;
			}
		}
		// finally, initialize a cell object for each set of cell file names
		cells=new hcs_fcs_man_cell[ncells];
		for(int i=0;i<(ncells-1);i++){
			String[] cellfiles=new String[cellfilenum[i]];
			for(int j=cellstartindex[i];j<cellstartindex[i+1];j++){
				cellfiles[j-cellstartindex[i]]=names2[j].substring(0);
			}
			cells[i]=new hcs_fcs_man_cell(cellfiles,this,0,0);
		}
		String[] cellfiles=new String[cellfilenum[ncells-1]];
		for(int j=cellstartindex[ncells-1];j<numfiles;j++){
			cellfiles[j-cellstartindex[ncells-1]]=names2[j].substring(0);
		}
		cells[ncells-1]=new hcs_fcs_man_cell(cellfiles,this,0,0);
	}

}
