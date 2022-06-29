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

public class custom_amfret_gate_generation_jru_v2 implements PlugIn {

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
		gd.addNumericField("Min_Concentration",1.0,5,15,null);
		gd.addNumericField("Max_Concentration",1000.0,5,15,null);
		gd.addNumericField("Min_AmFRET",-0.2,5,15,null);
		gd.addNumericField("Max_AmFRET",1.0,5,15,null);
		gd.addNumericField("Start_Gate (conc units)",2.0,5,15,null);
		gd.addNumericField("Upper_Gate_Conc",99.0,5,15,null);
		gd.addNumericField("Gate_Percentile",99.0,5,15,null);
		gd.addNumericField("Gate_Shift (AmFRET bins)",6,0);
		gd.addStringField("Acceptor Prefix","Cytosolic");
		gd.addStringField("AmFRET Prefix","Ratiometric");
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
		boolean showplots=gd.getNextBoolean();
		String accname=gd.getNextString();
		String fretname=gd.getNextString();
		amfret_utils au=new amfret_utils();
		au.accname=accname;
		au.fretname=fretname;
		au.getGate2(directory,name,outdir,roiname,minconc,maxconc,minamfret,maxamfret,startcrop,uppercentile,gateper,gateshift,showplots);
	}

}
