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

public class hcs_fcs_cell{
	// this class represents an hcs_fcs_cell
	// the cell has a single zoom scan, an fcs position, and eight raw fcs data
	// files
	public hcs_fcs_run[] runs;
	public String imagename,cellname;
	public int posx,posy,overview,image;
	public hcs_fcs_well parent;
	public double[] params;
	public boolean valid;

	public hcs_fcs_cell(hcs_fcs_run[] runs,String imagename,int posx,int posy,int overview,int image,String cellname,hcs_fcs_well parent){
		this.runs=runs;
		this.imagename=imagename;
		this.posx=posx;
		this.posy=posy;
		this.overview=overview;
		this.image=image;
		this.cellname=cellname;
		this.parent=parent;
	}

	public hcs_fcs_cell(String[] curvenames1,hcs_fcs_well parent1,int posx1,int posy1){
		parent=parent1;
		int ncurves=curvenames1.length;
		float[] runnumber=new float[ncurves];
		// first sort by run number
		for(int i=0;i<ncurves;i++){
			runnumber[i]=((new parse_filenames()).get_raw_fcs_info(curvenames1[i])[3]);
		}
		(new jsort()).sort_string_by(curvenames1,runnumber);
		// find out how many runs there are
		int currrun=(int)runnumber[0];
		int numruns=1;
		for(int i=1;i<ncurves;i++){
			if((int)runnumber[i]>currrun){
				numruns++;
				currrun=(int)runnumber[i];
			}
		}
		// find out where each run starts and how many files it has
		int[] runfilenum=new int[numruns];
		runfilenum[0]=1;
		int[] runstartindex=new int[numruns];
		runstartindex[0]=0;
		currrun=(int)runnumber[0];
		int runcounter=0;
		for(int i=1;i<ncurves;i++){
			if((int)runnumber[i]>currrun){
				currrun=(int)runnumber[i];
				runcounter++;
				runfilenum[runcounter]=1;
				runstartindex[runcounter]=i;
			}else{
				runfilenum[runcounter]++;
			}
		}
		// finally, initialize a run object for each set of run file names
		runs=new hcs_fcs_run[numruns];
		for(int i=0;i<(numruns-1);i++){
			String[] runfiles=new String[runfilenum[i]];
			for(int j=runstartindex[i];j<runstartindex[i+1];j++){
				runfiles[j-runstartindex[i]]=curvenames1[j].substring(0);
			}
			runs[i]=new hcs_fcs_run(runfiles,this);
		}
		String[] runfiles=new String[runfilenum[numruns-1]];
		for(int j=runstartindex[numruns-1];j<ncurves;j++){
			runfiles[j-runstartindex[numruns-1]]=curvenames1[j].substring(0);
		}
		runs[numruns-1]=new hcs_fcs_run(runfiles,this);
		posx=posx1;
		posy=posy1;
		overview=new parse_filenames().get_raw_fcs_info(curvenames1[0])[1];
		image=new parse_filenames().get_raw_fcs_info(curvenames1[0])[2];
		imagename=new parse_filenames().generate_zname_from_rawname(curvenames1[0]);
		cellname="Cell"+(overview*10+image);
	}

	public void update_params(int stat){
		int nparams=runs[0].params.length;
		params=new double[nparams];
		int nvalid=0;
		for(int i=0;i<runs.length;i++){
			if(runs[i].valid){
				nvalid++;
			}
		}
		if(nvalid==0){
			valid=false;
			return;
		}
		valid=true;
		for(int i=0;i<nparams;i++){
			float[] temp=new float[nvalid];
			int counter=0;
			for(int j=0;j<runs.length;j++){
				if(runs[j].valid){
					temp[counter]=(float)runs[j].params[i];
					counter++;
				}
			}
			params[i]=jstatistics.getstatistic(jstatistics.stats[stat],temp,null);
			// params[i]=1.0;
		}
	}

}
