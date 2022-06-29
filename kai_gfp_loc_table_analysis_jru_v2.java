/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.text.*;
import java.util.*;
import jguis.*;
import jalgs.*;
import jalgs.jsim.*;

public class kai_gfp_loc_table_analysis_jru_v2 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		//columns are 0image,1object,2area,3avg,4stdev,5puntarea,6punctavg,7membavg,8hoechstavg,9nucarea,10nucavg,11nucpuntarea,12nucpuntavg,13cellmembarea,14cellmembavg,15nucmembarea,16nucmembavg
		//need to calculate 17area ratio, 18wellid, 19row, 20col,21temp,22strain,23autoavg,24cytoavg,25nuc_ratio,26memb_ratio,27nucmemb_ratio,28cellcount
		String[] newlabels={"area_ratio","wellid","row","col","temp","strain","autoavg","cytoavg","nuc_ratio","memb_ratio","nucmemb_ratio","cellct"};
		String[] newlabels2={"area_ratio","wellid","row","col","temp","strain","autoavg","cytoavg"};
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Auto_23_deg",0.0,5,15,null);
		gd.addNumericField("Auto_35_deg",0.0,5,15,null);
		gd.addNumericField("Min_Area",150,5,15,null);
		gd.addNumericField("Max_Area",600,5,15,null);
		gd.addNumericField("Min_Area_Ratio",0.05,5,15,null);
		gd.addNumericField("Max_Area_Ratio",0.4,5,15,null);
		gd.addNumericField("Ref_23_row",3,0);
		gd.addNumericField("Ref_35_row",4,0);
		gd.addNumericField("Ref_23_col",1,0);
		gd.addNumericField("Ref_35_col",1,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float auto23=(float)gd.getNextNumber();
		float auto35=(float)gd.getNextNumber();
		float minarea=(float)gd.getNextNumber();
		float maxarea=(float)gd.getNextNumber();
		float minarearatio=(float)gd.getNextNumber();
		float maxarearatio=(float)gd.getNextNumber();
		int refrow23=(int)gd.getNextNumber()-1;
		int refrow35=(int)gd.getNextNumber()-1;
		int refcol23=(int)gd.getNextNumber()-1;
		int refcol35=(int)gd.getNextNumber()-1;
		//start by making the nuc area ratio column
		float[] area=table_tools.get_column_array(listtable,2);
		float[] nucarea=table_tools.get_column_array(listtable,9);
		float[] area_ratio=new float[area.length];
		List<List<String>> filtered=new ArrayList<List<String>>();
		for(int i=0;i<area.length;i++){
			area_ratio[i]=nucarea[i]/area[i];
			if(area[i]>minarea && area[i]<maxarea && area_ratio[i]>minarearatio && area_ratio[i]<maxarearatio){
				List<String> row=listtable.get(i);
				row.add(""+area_ratio[i]);
				String wellid=row.get(0).substring(0,6);
				int nrow=Integer.parseInt(wellid.substring(0,3));
				int col=Integer.parseInt(wellid.substring(3,6));
				int temp=23+((nrow-1)%2)*12;
				int strain=22*(int)Math.floor((nrow-1)/2)+col-1;
				row.add(wellid);
				row.add(""+nrow);
				row.add(""+col);
				row.add(""+temp);
				row.add(""+strain);
				float autoavg=auto35;
				if(temp==23) autoavg=auto23;
				row.add(""+autoavg);
				float nucavg=table_tools.get_number(listtable,i,10)-autoavg; row.set(10,""+nucavg);
				float avg=table_tools.get_number(listtable,i,3)-autoavg; row.set(3,""+avg);
				float cellmembavg=table_tools.get_number(listtable,i,14)-autoavg; row.set(14,""+cellmembavg);
				float nucmembavg=table_tools.get_number(listtable,i,16)-autoavg; row.set(16,""+nucmembavg);
				float cytoavg=(avg*area[i]-nucavg*nucarea[i])/(area[i]-nucarea[i]);
				row.add(""+cytoavg);
				//now calculate the cell by cell ratios (maybe better to do after averaging?)
				/*float nuc_ratio=nucavg/avg;
				float memb_ratio=cellmembavg*(area[i]-nucarea[i])/(avg*area[i]-nucavg*nucarea[i]); //ratio with respect to cytoplasm
				float nucmemb_ratio=nucmembavg/nucavg; //ratio with respect to nucleus
				row.add(""+nuc_ratio);
				row.add(""+memb_ratio);
				row.add(""+nucmemb_ratio);*/
				filtered.add(row);
			}
		}
		listtable=null;
		table_tools.sort_listtable(filtered,18);
		String[] col_labels3=new String[col_labels.length+newlabels2.length];
		System.arraycopy(col_labels,0,col_labels3,0,col_labels.length);
		System.arraycopy(newlabels2,0,col_labels3,col_labels.length,newlabels2.length);
		String prefix=tw[0].getTitle();
		if(prefix.indexOf(".")>0) prefix=prefix.substring(0,prefix.indexOf("."));
		table_tools.create_table(prefix+"_filtered.xls",filtered,col_labels3);
		List<List<String>> avglist=table_tools.get_cell_stat_list(filtered,18,"Avg",true);
		List<List<String>> semlist=table_tools.get_cell_stat_list(filtered,18,"StErr",true);
		//the sem values aren't right for temp and strain
		for(int i=0;i<avglist.size();i++){
			semlist.get(i).set(21,avglist.get(i).get(21));
			semlist.get(i).set(22,avglist.get(i).get(22));
		}
		//sort by temp and then strain
		table_tools.sort_listtable(avglist,21,true);
		table_tools.sort_listtable(avglist,22,true);
		table_tools.sort_listtable(semlist,21,true);
		table_tools.sort_listtable(semlist,22,true);
		//table_tools.create_table(prefix+"_tavg.xls",avglist,null);
		//table_tools.create_table(prefix+"_tsem.xls",semlist,null);

		//optionally add the nuclear and membrane ratios after the fact
		for(int i=0;i<avglist.size();i++){
			float nucavg=table_tools.get_number(avglist,i,10);
			float avg=table_tools.get_number(avglist,i,3);
			float cellmembavg=table_tools.get_number(avglist,i,14);
			float nucmembavg=table_tools.get_number(avglist,i,16);
			float cytoavg=table_tools.get_number(avglist,i,24);

			float nucavgsem=table_tools.get_number(semlist,i,10);
			float avgsem=table_tools.get_number(semlist,i,3);
			float cellmembavgsem=table_tools.get_number(semlist,i,14);
			float nucmembavgsem=table_tools.get_number(semlist,i,16);
			float cytoavgsem=table_tools.get_number(semlist,i,24);

			float nucarea1=table_tools.get_number(avglist,i,9);
			float area1=table_tools.get_number(avglist,i,2);
			float cellct=table_tools.get_number(avglist,i,25);
			avglist.get(i).remove(25);
			semlist.get(i).remove(25);

			float nuc_ratio=nucavg/avg; 
			float nrsem=propRatioErrs(nucavg,avg,nucavgsem,avgsem);
			avglist.get(i).add(""+nuc_ratio); semlist.get(i).add(""+nrsem);
			float memb_ratio=cellmembavg/cytoavg;
			float mrsem=propRatioErrs(cellmembavg,cytoavg,cellmembavgsem,cytoavgsem);
			avglist.get(i).add(""+memb_ratio); semlist.get(i).add(""+mrsem);
			float nm_ratio=nucmembavg/nucavg;
			float nmrsem=propRatioErrs(nucmembavg,nucavg,nucmembavgsem,nucavgsem);
			avglist.get(i).add(""+nm_ratio); semlist.get(i).add(""+nmrsem);
			avglist.get(i).add(""+cellct);
			semlist.get(i).add(""+cellct);
		}

		String[] col_labels2=new String[col_labels.length+newlabels.length];
		System.arraycopy(col_labels,0,col_labels2,0,col_labels.length);
		System.arraycopy(newlabels,0,col_labels2,col_labels.length,newlabels.length);
		table_tools.create_table(prefix+"_avg.xls",avglist,col_labels2);
		table_tools.create_table(prefix+"_sem.xls",semlist,col_labels2);
		filtered=null;
		//columns are 0image,1object,2area,3avg,4stdev,5puntarea,6punctavg,7membavg,8hoechstavg,9nucarea,10nucavg,11nucpuntarea,12nucpuntavg,13cellmembarea,14cellmembavg,15nucmembarea,16nucmembavg
		//17area ratio, 18wellid, 19row, 20col,21temp,22strain,23autoavg,24cytoavg,25nuc_ratio,26memb_ratio,27nucmemb_ratio,28cellcount
		//now calculate the temperature comparisons
		//this table needs name, wellid, strain, expression diff, nuclear enrichment diff, puncta area diff, membrane enrichment diff, and nuclear membrane enrichment diff
		//as well as propagated errors for these
		String[] ratiolabels={"name23","wellid23","name35","wellid35","strain","avg25","semavg","expdiff","sem1","nucendiff","sem2","punctareadiff","sem3","membendiff","sem4","nmendiff","sem5"};
		List<List<String>> ratiolist=new ArrayList<List<String>>();
		for(int i=0;i<avglist.size();i++){
			if(table_tools.get_number(avglist,i,21)==23.0f){
				if(table_tools.get_number(avglist,i+1,21)==35.0f){
					String name23=avglist.get(i).get(0);
					String wellid23=avglist.get(i).get(18);
					String name35=avglist.get(i+1).get(0);
					String wellid35=avglist.get(i+1).get(18);
					String strain=avglist.get(i).get(22);
					List<String> row=new ArrayList<String>();
					row.add(name23); row.add(wellid23); row.add(name35); row.add(wellid35); row.add(strain);
					float[] arr23=table_tools.get_row_array(avglist,i);
					float[] arr35=table_tools.get_row_array(avglist,i+1);
					float[] arr23sem=table_tools.get_row_array(semlist,i);
					float[] arr35sem=table_tools.get_row_array(semlist,i+1);
					int c=3; //expression
					float val=arr35[c]-arr23[c];
					float sem=propSumErrs(arr35[c],arr23[c],arr35sem[c],arr23sem[c]);
					row.add(""+arr23[c]); row.add(""+arr23sem[c]);
					row.add(""+val); row.add(""+sem);
					c=25; //nuclear enrichment
					val=arr35[c]-arr23[c];
					sem=propSumErrs(arr35[c],arr23[c],arr35sem[c],arr23sem[c]);
					row.add(""+val); row.add(""+sem);
					c=5; //puncta area
					val=arr35[c]-arr23[c];
					sem=propSumErrs(arr35[c],arr23[c],arr35sem[c],arr23sem[c]);
					row.add(""+val); row.add(""+sem);
					c=26; //membrane enrichment
					val=arr35[c]-arr23[c];
					sem=propSumErrs(arr35[c],arr23[c],arr35sem[c],arr23sem[c]);
					row.add(""+val); row.add(""+sem);
					c=27; //nm enrichment
					val=arr35[c]-arr23[c];
					sem=propSumErrs(arr35[c],arr23[c],arr35sem[c],arr23sem[c]);
					row.add(""+val); row.add(""+sem);
					ratiolist.add(row);
				}
			}
		}
		//need to fill in the missing strains now
		//strain list should start at 0 and go at least to 175
		//insert a blank strain in the missing slots
		//columns go up two at a time
		for(int i=0;i<ratiolist.size();i++){
			int strain=(int)table_tools.get_number(ratiolist,i,4);
			while(i<strain){
				ratiolist.add(i,make_dummy_strain(i));
				i++;
			}
		}
		table_tools.create_table(prefix+"_derived.xls",ratiolist,ratiolabels);
		IJ.run("heatmap jru v1", "windows="+prefix+"_avg.xls x=col y=row intensity=avg");
		ImagePlus imp=WindowManager.getCurrentImage();
		int temp=imp.getWidth();
		float[] pix=(float[])imp.getProcessor().getPixels();
		float refval23=pix[refrow23*4*temp+refcol23*4];
		float refval35=pix[refrow35*4*temp+refcol35*4];
		IJ.log("ref 23 = "+refval23);
		IJ.log("ref 35 = "+refval35);
	}

	public List<String> make_dummy_strain(int strainnum){
		int col=strainnum%22+1;
		int row=2*(int)Math.floor(strainnum/22)+1;
		ArrayList<String> temp=new ArrayList<String>();
		String srow="0"+row;
		if(row<10) srow="00"+row;
		String scol="0"+col;
		if(col<10) scol="00"+col;
		int row2=row+1;
		String srow2="0"+row2;
		if(row2<10) srow2="00"+row2;
		String wellid1=srow+scol;
		String wellid2=srow2+scol;
		temp.add(wellid1+"-1-001001004.tif");
		temp.add(wellid1);
		temp.add(wellid2+"-1-001001004.tif");
		temp.add(wellid2);
		temp.add(""+strainnum);
		for(int i=0;i<12;i++) temp.add("0");
		return temp;
	}

	public float propRatioErrs(float num,float den,float numsem,float densem){
		rngs random=new rngs();
		float[] temp=new float[100];
		for(int i=0;i<100;i++) temp[i]=(float)random.gasdev(num,numsem)/(float)random.gasdev(den,densem);
		return jstatistics.getstatistic("StDev",temp,null);
	}

	/*public float propRatioErrs(float num,float den,float numsem,float densem){
		return (num/den)*(float)Math.sqrt(numsem*numsem/(num*num)+(densem*densem)/(den*den));
	}*/

	public float propSumErrs(float val1,float val2,float val1sem,float val2sem){
		return (float)Math.sqrt(val1sem*val1sem+val2sem*val2sem);
	}

}
