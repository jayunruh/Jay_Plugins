/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

import java.io.File;

public class hcs_fcs_run{
	public String[] filenames;
	public int runnumber;
	public double[] params;
	public int id;
	public hcs_fcs_cell parent;
	public boolean valid;

	public hcs_fcs_run(String[] filenames,int runnumber,hcs_fcs_cell parent){
		this.filenames=filenames;
		this.runnumber=runnumber;
		this.parent=parent;
		valid=true;
	}

	public hcs_fcs_run(String[] filenames1,hcs_fcs_cell parent1){
		parent=parent1;
		runnumber=new parse_filenames().get_raw_fcs_run(filenames1[0]);
		if(filenames1.length<2){
			if(filenames1==null){
				filenames=null;
			}else{
				filenames=new String[1];
				filenames[0]=filenames1[0].substring(0);
			}
		}else{
			filenames=new String[2];
			int ch1=(new parse_filenames()).get_raw_fcs_channel(filenames1[0]);
			if(ch1==1){
				filenames[0]=filenames1[0].substring(0);
				filenames[1]=filenames1[1].substring(0);
			}else{
				filenames[1]=filenames1[0].substring(0);
				filenames[0]=filenames1[1].substring(0);
			}
		}
		valid=true;
	}

	public File getfile(int channel){
		// note that the channel here is one based
		String sep=File.separator;
		if(channel==1){
			return new File(parent.parent.parent.dir+sep+parent.parent.foldername+sep+filenames[channel-1]);
		}else{
			if(filenames.length>1){
				return new File(parent.parent.parent.dir+sep+parent.parent.foldername+sep+filenames[channel-1]);
			}else{
				return null;
			}
		}
	}

	public String getpath(int channel){
		// note that the channel here is one based
		String sep=File.separator;
		if(channel==1){
			return parent.parent.parent.dir+sep+parent.parent.foldername+sep+filenames[channel-1];
		}else{
			if(filenames.length>1){
				return parent.parent.parent.dir+sep+parent.parent.foldername+sep+filenames[channel-1];
			}else{
				return null;
			}
		}
	}

}
