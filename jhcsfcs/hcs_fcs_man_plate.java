/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

import jalgs.jdataio;

import java.io.File;

public class hcs_fcs_man_plate extends hcs_fcs_plate{

	public hcs_fcs_man_plate(String directory,boolean singlewell){
		super(null,directory);
		// here we generate a nested table of hcsfcs filenames
		if(!singlewell){
			// first get the list of well folders
			dir=directory.substring(0);
			String[] wellnames=(new jdataio()).get_sorted_string_list(dir,null);
			boolean[] valid=new boolean[wellnames.length];
			int nvalid=0;
			for(int i=0;i<wellnames.length;i++){
				if((new File(dir+File.separator+wellnames[i])).isDirectory()&&countraw(wellnames[i])>0){
					valid[i]=true;
					nvalid++;
				}else{
					valid[i]=false;
				}
			}
			// eliminate all non-directory folders
			int counter=0;
			String[] wellnames2=new String[nvalid];
			for(int i=0;i<wellnames.length;i++){
				if(valid[i]){
					wellnames2[counter]=wellnames[i].substring(0);
					counter++;
				}
			}
			// now create well objects for each well
			wells=new hcs_fcs_man_well[wellnames2.length];
			for(int i=0;i<wellnames2.length;i++){
				wells[i]=new hcs_fcs_man_well(wellnames2[i],this);
			}
		}else{
			File dirfile=new File(directory);
			dir=dirfile.getParent();
			String name=dirfile.getName();
			wells=new hcs_fcs_man_well[1];
			wells[0]=new hcs_fcs_man_well(name,this);
		}
	}

	public hcs_fcs_man_plate(String[] filelist){
		super(null,null);
		File temp=new File(filelist[0]);
		File dirfile=temp.getParentFile();
		dir=dirfile.getParent();
		String name=dirfile.getName();
		wells=new hcs_fcs_man_well[1];
		wells[0]=new hcs_fcs_man_well(filelist,name,this);
	}

}
