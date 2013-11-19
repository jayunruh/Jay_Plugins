/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jhcsfcs;

public class parse_filenames{

	public int get_folder_well_number(String folder){
		int location=folder.indexOf("_Well_");
		location+=6;
		int location2=folder.indexOf(".",location+1);
		String wellstring=folder.substring(location,location2);
		return Integer.parseInt(wellstring);
	}

	public int[] get_raw_fcs_info(String rawfilename){
		// here we output an integer array with well#, overview#, image#, run#,
		// and channel#
		// the date MUST be in format mon_da_year or at least have 11 characters
		// or there must be no date
		// the general format is
		// Name_fcs_Mon_da_year_well#_overview#_image#_zeiss's
		// junk_R#_P#_K#_Ch#.raw
		int[] outarray=new int[5];
		int location=rawfilename.indexOf("_fcs_");
		if(Character.isLetter(rawfilename.charAt(location+5))){
			location+=17;
		}else{
			location+=5;
		}
		int location2=rawfilename.indexOf("_",location+1);
		outarray[0]=Integer.parseInt(rawfilename.substring(location,location2));
		location2++;
		int location3=rawfilename.indexOf("_",location2+1);
		outarray[1]=Integer.parseInt(rawfilename.substring(location2,location3));
		location3++;
		int location4=rawfilename.indexOf("_",location3+1);
		outarray[2]=Integer.parseInt(rawfilename.substring(location3,location4));
		int location5=rawfilename.indexOf("_",location4+3);
		location5+=2;
		int location6=rawfilename.indexOf("_",location5+1);
		outarray[3]=Integer.parseInt(rawfilename.substring(location5,location6));
		int location7=rawfilename.length()-5;
		int location8=rawfilename.indexOf(".",location7+1);
		outarray[4]=Integer.parseInt(rawfilename.substring(location7,location8));
		return outarray;
	}

	public String get_raw_fcs_date(String rawfilename){
		int location=rawfilename.indexOf("_fcs_");
		location+=5;
		if(Character.isLetter(rawfilename.charAt(location))){
			return rawfilename.substring(location,location+11);
		}else{
			return null;
		}
	}

	public String get_raw_fcs_orf(String rawfilename){
		int location=rawfilename.indexOf("_fcs_");
		return rawfilename.substring(0,location);
	}

	public String generate_zname_from_rawname(String rawfilename){
		String orfname=get_raw_fcs_orf(rawfilename);
		String date=get_raw_fcs_date(rawfilename);
		int[] indices=get_raw_fcs_info(rawfilename);
		if(date!=null){
			return "cr_"+orfname+"_zoomf_"+date+"_"+indices[0]+"_"+indices[1]+"_"+indices[2]+".lsm";
		}else{
			return "cr_"+orfname+"_zoomf_"+indices[0]+"_"+indices[1]+"_"+indices[2]+".lsm";
		}
	}

	public int get_raw_fcs_run(String rawfilename){
		// work backwards from the end of the name to get the run number
		String[] delim=rawfilename.split("_");
		int length=delim.length;
		return Integer.parseInt(delim[length-4].substring(1));
	}

	public int get_raw_fcs_channel(String rawfilename){
		String[] delim=rawfilename.split("_");
		String[] delim2=delim[delim.length-1].split("\\.");
		return Integer.parseInt(delim2[0].substring(2));
	}

	public String get_raw_fcs_name(String rawfilename){
		// here we work backwards from the end of the file name to get
		// everything before the zeiss id string
		String[] delim=rawfilename.split("_");
		int length=delim.length;
		StringBuffer sb=new StringBuffer();
		sb.append(delim[0]);
		for(int i=1;i<(length-5);i++){
			sb.append("_"+delim[i]);
		}
		return sb.toString();
	}
}
