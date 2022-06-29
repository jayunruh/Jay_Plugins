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

public class custom_amfret_gate_generation_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin loads an fcs file to do custom amfret analysis
		//start by selecting the fcs file
		OpenDialog od = new OpenDialog("Open File",arg);
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		//now select the output directory for the roi
		SaveDialog sd=new SaveDialog("Save Roi",arg,".roi");
		String outdir=sd.getDirectory();
		if(outdir==null) return;
		String roiname=sd.getFileName();
		//now get the analysis paremeters
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_Concentration (millions)",0.5,5,15,null);
		gd.addNumericField("Max_Concentration (millions)",500.0,5,15,null);
		gd.addNumericField("Min_AmFRET",-0.2,5,15,null);
		gd.addNumericField("Max_AmFRET",1.0,5,15,null);
		gd.addNumericField("Start_Gate (conc units)",1.5,5,15,null);
		gd.addNumericField("Upper_Gate_Conc",99.0,5,15,null);
		gd.addNumericField("Gate_Percentile",99.0,5,15,null);
		gd.addNumericField("Gate_Shift (AmFRET bins)",6,0);
		gd.addStringField("Acceptor_Name","FL10-A");
		gd.addStringField("FRET_Name","FL04-A");
		gd.addCheckbox("Show_Plots",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float minconc=(float)gd.getNextNumber();
		float maxconc=(float)gd.getNextNumber();
		float minamfret=(float)gd.getNextNumber();
		float maxamfret=(float)gd.getNextNumber();
		float startcrop=(float)gd.getNextNumber();
		float uppercentile=(float)gd.getNextNumber();
		float gateper=(float)gd.getNextNumber();
		int gateshift=(int)gd.getNextNumber();
		String accname=gd.getNextString();
		String fretname=gd.getNextString();
		boolean showplots=gd.getNextBoolean();
		amfret_utils au=new amfret_utils();
		au.accname=accname;
		au.fretname=fretname;
		au.getGate(directory,name,outdir,roiname,minconc,maxconc,minamfret,maxamfret,startcrop,uppercentile,gateper,gateshift,showplots);
	}

}
