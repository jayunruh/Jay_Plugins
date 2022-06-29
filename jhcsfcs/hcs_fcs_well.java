/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

import jalgs.jsort;
import jalgs.jstatistics;

import java.io.File;

public class hcs_fcs_well{
	// this is a container class for each hcs_fcs well
	public hcs_fcs_cell[] cells;
	public String orfname,foldername;
	public int wellnumber;
	public hcs_fcs_plate parent;
	public double[] params;
	public boolean valid;

	public hcs_fcs_well(hcs_fcs_cell[] cells,String orfname,String foldername,int wellnumber,hcs_fcs_plate parent){
		this.cells=cells;
		this.orfname=orfname;
		this.foldername=foldername;
		this.wellnumber=wellnumber;
		this.parent=parent;
	}

	public hcs_fcs_well(String well,hcs_fcs_plate parent1){
		foldername=well;
		parent=parent1;
		wellnumber=new parse_filenames().get_folder_well_number(well);
		// get the list of files for this well
		String[] names=(new File(parent.dir+File.separator+well)).list();
		// figure out what file names we have
		int numfiles=0;
		int[][] rawfileinfo=null;
		if(names!=null){
			rawfileinfo=new int[names.length][];
			for(int i=0;i<names.length;i++){
				if(names[i].indexOf(".raw")>0&&names[i].indexOf("Marker")<0){
					numfiles++;
					rawfileinfo[i]=(new parse_filenames()).get_raw_fcs_info(names[i]);
				}else{
					rawfileinfo[i]=null;
				}
			}
		}
		if(numfiles>0){
			// weed out the names we don't need
			float[] cellnumber=new float[numfiles];
			String[] names2=new String[numfiles];
			int counter=0;
			for(int i=0;i<names.length;i++){
				if(rawfileinfo[i]!=null){
					names2[counter]=names[i].substring(0);
					// the cell number is the overview number*100 plus the image
					// number
					cellnumber[counter]=rawfileinfo[i][1]*100+rawfileinfo[i][2];
					counter++;
				}
			}
			orfname=(new parse_filenames()).get_raw_fcs_orf(names2[0]);
			// now we need to sort the file names by cell number
			(new jsort()).sort_string_by(names2,cellnumber);
			// now sort the cellnumber array itself
			jsort.javasort_order(cellnumber);
			// find out how many cells we have
			int currcell=(int)cellnumber[0];
			int numcells=1;
			for(int i=1;i<numfiles;i++){
				if((int)cellnumber[i]>currcell){
					currcell=(int)cellnumber[i];
					numcells++;
				}
			}
			// find out where each cell starts and how many files it has
			int[] cellfilenum=new int[numcells];
			cellfilenum[0]=1;
			int[] cellstartindex=new int[numcells];
			cellstartindex[0]=0;
			currcell=(int)cellnumber[0];
			int cellcounter=0;
			for(int i=1;i<numfiles;i++){
				if((int)cellnumber[i]>currcell){
					currcell=(int)cellnumber[i];
					cellcounter++;
					cellfilenum[cellcounter]=1;
					cellstartindex[cellcounter]=i;
				}else{
					cellfilenum[cellcounter]++;
				}
			}
			// finally, initialize a cell object for each set of cell file names
			cells=new hcs_fcs_cell[numcells];
			for(int i=0;i<(numcells-1);i++){
				String[] cellfiles=new String[cellfilenum[i]];
				for(int j=cellstartindex[i];j<cellstartindex[i+1];j++){
					cellfiles[j-cellstartindex[i]]=names2[j].substring(0);
				}
				cells[i]=new hcs_fcs_cell(cellfiles,this,0,0);
			}
			String[] cellfiles=new String[cellfilenum[numcells-1]];
			for(int j=cellstartindex[numcells-1];j<numfiles;j++){
				cellfiles[j-cellstartindex[numcells-1]]=names2[j].substring(0);
			}
			cells[numcells-1]=new hcs_fcs_cell(cellfiles,this,0,0);
		}else{
			cells=null;
		}
	}

	public void update_params(int stat){
		int nparams=cells[0].params.length;
		params=new double[nparams];
		int nvalid=0;
		for(int i=0;i<cells.length;i++){
			if(cells[i].valid){
				nvalid++;
			}
		}
		if(nvalid==0){
			valid=false;
			return;
		}
		for(int i=0;i<nparams;i++){
			float[] temp=new float[nvalid];
			int counter=0;
			for(int j=0;j<cells.length;j++){
				if(cells[j].valid){
					temp[counter]=(float)cells[j].params[i];
					counter++;
				}
			}
			params[i]=jstatistics.getstatistic(jstatistics.stats[stat],temp,null);
		}
	}

}
