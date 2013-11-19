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

public class hcs_fcs_plate{
	public hcs_fcs_well[] wells;
	public String dir;

	public hcs_fcs_plate(hcs_fcs_well[] wells,String directory){
		this.wells=wells;
		this.dir=directory;
	}

	public hcs_fcs_plate(String directory){
		// here we generate a nested table of hcsfcs filenames
		// first get the list of well folders
		dir=directory.substring(0);
		String[] wells1=(new File(directory)).list();
		int[] wellnumber=new int[wells1.length];
		int numwells=0;
		for(int i=0;i<wells1.length;i++){
			if(wells1[i].indexOf("_Well_")>0&&countraw(wells1[i])>0){
				numwells++;
				wellnumber[i]=(new parse_filenames()).get_folder_well_number(wells1[i]);
			}else{
				wellnumber[i]=-1;
			}
		}
		// now copy the valid well names to a new array
		int counter=0;
		String[] wells2=new String[numwells];
		float[] wellnumber2=new float[numwells];
		for(int i=0;i<wells1.length;i++){
			if(wellnumber[i]>=0){
				wells2[counter]=wells1[i].substring(0);
				wellnumber2[counter]=(float)wellnumber[i];
				counter++;
			}
		}
		// now sort the well names using the well number array
		(new jsort()).sort_string_by(wells2,wellnumber2);
		// now create well objects for each well
		wells=new hcs_fcs_well[numwells];
		for(int i=0;i<numwells;i++){
			wells[i]=new hcs_fcs_well(wells2[i],this);
		}
	}

	public int getnfiles(){
		int nfiles=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					for(int l=0;l<runs[k].filenames.length;l++){
						nfiles++;
					}
				}
			}
		}
		return nfiles;
	}

	public File[] getfilelist(){
		File[] list=new File[getnfiles()];
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					for(int l=0;l<runs[k].filenames.length;l++){
						list[counter]=runs[k].getfile(l+1);
						counter++;
					}
				}
			}
		}
		return list;
	}

	public String[] getpathlist(){
		String[] list=new String[getnfiles()];
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					for(int l=0;l<runs[k].filenames.length;l++){
						list[counter]=runs[k].getpath(l+1);
						counter++;
					}
				}
			}
		}
		return list;
	}

	public String[][] getplatetable(){
		// the table columns are well name, cell name, run#, and path
		String[][] list=new String[getnfiles()][4];
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					for(int l=0;l<runs[k].filenames.length;l++){
						list[counter][0]=wells[i].foldername.substring(0);
						list[counter][1]=cells[j].cellname.substring(0);
						list[counter][2]=""+runs[k].runnumber;
						list[counter][3]=runs[k].getpath(l+1);
						counter++;
					}
				}
			}
		}
		return list;
	}

	public void generate_ids(){
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					runs[k].id=counter;
					counter++;
				}
			}
		}
	}

	public hcs_fcs_run getrunforid(int id){
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					if(counter==id){
						return runs[k];
					}
					counter++;
				}
			}
		}
		return null;
	}

	public int[][] getplatenumbers(){
		// the table columns are well#, cell#, and run#
		int[][] list=new int[getnfiles()][3];
		int counter=0;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					for(int l=0;l<runs[k].filenames.length;l++){
						list[counter][0]=wells[i].wellnumber;
						list[counter][1]=cells[j].overview*10+cells[j].image;
						list[counter][2]=runs[k].runnumber;
						counter++;
					}
				}
			}
		}
		return list;
	}

	public void filter_data(int paramnumber,int filtertype,double filtervalue){
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					if(runs[k].valid){
						if(filtertype==0){
							runs[k].valid=(runs[k].params[paramnumber]>filtervalue);
						}else{
							runs[k].valid=(runs[k].params[paramnumber]<filtervalue);
						}
					}
				}
			}
		}
	}

	public void unfilter_data(int paramnumber,int filtertype,double filtervalue){
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					if(!runs[k].valid){
						if(filtertype==0){
							runs[k].valid=!(runs[k].params[paramnumber]>filtervalue);
						}else{
							runs[k].valid=!(runs[k].params[paramnumber]<filtervalue);
						}
					}
				}
			}
		}
	}

	public int countraw(String foldername){
		String[] names=new File(dir+File.separator+foldername).list();
		int nraw=0;
		for(int i=0;i<names.length;i++){
			if(names[i].endsWith(".raw")&&names[i].indexOf("Marker")<0){
				nraw++;
			}
		}
		return nraw;
	}
}
