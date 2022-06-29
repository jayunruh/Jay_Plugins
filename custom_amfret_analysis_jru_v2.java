/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;
import jguis.*;
import jalgs.jfit.*;
import ij.io.*;
import ij.plugin.frame.RoiManager;
import ij.text.*;
import java.io.*;
import java.util.*;

public class custom_amfret_analysis_jru_v2 implements PlugIn {
	//this version loads an fcs file with 

	public void run(String arg) {
		//this plugin loads an fcs file to do custom amfret analysis
		//start by selecting the fcs file
		OpenDialog od = new OpenDialog("Open File",arg);
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		//now select the output directory
		DirectoryChooser dc = new DirectoryChooser("Save Directory");
		dc.setDefaultDirectory(directory);
		String outdir=dc.getDirectory();
		if(outdir==null){return;}
		//now get the analysis paremeters
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_Concentration",1.0,5,15,null);
		gd.addNumericField("Max_Concentration",1000.0,5,15,null);
		gd.addNumericField("Min_AmFRET",-0.2,5,15,null);
		gd.addNumericField("Max_AmFRET",1.0,5,15,null);
		gd.addNumericField("Min_Cells",5,0);
		gd.addNumericField("Start_Crop (conc units)",2.0,5,15,null);
		gd.addNumericField("Min_Fraction_Bimodal",0.02,5,15,null);
		gd.addNumericField("Max_Fraction_Bimodal",0.95,5,15,null);
		gd.addCheckbox("Show_Plots",false);
		gd.addCheckbox("load_gate_roi",false);
		gd.addStringField("gate_roi_path","path_to_roi");
		gd.addStringField("Acceptor_Prefix","Cytosolic");
		gd.addStringField("AmFRET_Prefix","Ratiometric");
		gd.showDialog(); if(gd.wasCanceled()) return;
		float minconc=(float)gd.getNextNumber();
		float maxconc=(float)gd.getNextNumber();
		float minamfret=(float)gd.getNextNumber();
		float maxamfret=(float)gd.getNextNumber();
		int mincells=(int)gd.getNextNumber();
		float startcrop=(float)gd.getNextNumber();
		float minbimodefrac=(float)gd.getNextNumber();
		float maxbimodefrac=(float)gd.getNextNumber();
		boolean showplots=gd.getNextBoolean();
		boolean loadroi=gd.getNextBoolean();
		String roipath=gd.getNextString();
		if(!loadroi) roipath=null;
		String accname=gd.getNextString();
		String fretname=gd.getNextString();
		amfret_utils au=new amfret_utils();
		au.accname=accname;
		au.fretname=fretname;
		List<String> output=au.exec2(directory,name,outdir,roipath,minconc,maxconc,minamfret,maxamfret,mincells,startcrop,minbimodefrac,maxbimodefrac,showplots);
		//output=output.subList(0,20);
		//String[] col_labels={"file","datFile","well","plate","acceptor","c^2","Iter","baseline","amp","EC50","alpha","xshift","EC50_errs","alpha_errs","totcells","fretcells","bimodal_metric","f_gate","delta","delta_errs"};
		String[] col_labels={"file","datFile","well","plate","acceptor","c^2","Iter","baseline","amp","EC50","alpha","xshift","EC50_errs","alpha_errs","totcells","fretcells","bimodal_metric","f_gate","delta","delta_errs","accstdev","nfaccmean","nfaccstdev","faccmean","faccstdev","fretmean","fretstdev"};
		TextWindow tw=jutils.selectTable("Stretched Exp Fits");
		if(tw==null) tw=new TextWindow("Stretched Exp Fits",table_tools.print_string_array(col_labels),"",400,200);
		if(output!=null) tw.append(table_tools.print_string_array(output,0));
	}

}
