/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

import jalgs.jsort;

public class hcs_fcs_man_cell extends hcs_fcs_cell{

	public hcs_fcs_man_cell(String[] curvenames1,hcs_fcs_well parent1,int posx1,int posy1){
		super(null,null,posx1,posy1,0,0,null,parent1);
		int ncurves=curvenames1.length;
		float[] runnumber=new float[ncurves];
		// first sort by run number
		for(int i=0;i<ncurves;i++){
			runnumber[i]=((new parse_filenames()).get_raw_fcs_run(curvenames1[i]));
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
		cellname=new parse_filenames().get_raw_fcs_name(curvenames1[0]);
	}

}
