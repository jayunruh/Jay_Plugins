/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import ij.io.*;

public class generate_jru_luts_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] luts={"NICE","NICE_whiteback","log_gray","log_red","log_green","log_blue","greenyellowred","bluegreenred","angle","rand"};
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("LUT",luts,luts[0]);
		gd.addCheckbox("Output All",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int index=gd.getNextChoiceIndex();
		boolean outall=gd.getNextBoolean();
		//first get the directory
		DirectoryChooser.setDefaultDirectory(IJ.getDirectory("luts"));
		DirectoryChooser od=new DirectoryChooser("Choose LUT Directory");
		String dir=od.getDirectory();
		if(dir==null) return;
		if(outall){
			byte[][] lut=lututils.nice_lut(false);
			lututils.write_lut(lut,dir,"NICE.lut");
			lut=lututils.nice_lut(true);
			lututils.write_lut(lut,dir,"NICE_whiteback.lut");
			lut=lututils.log_lut(0);
			lututils.write_lut(lut,dir,"log_gray.lut");
			lut=lututils.log_lut(1);
			lututils.write_lut(lut,dir,"log_red.lut");
			lut=lututils.log_lut(2);
			lututils.write_lut(lut,dir,"log_green.lut");
			lut=lututils.log_lut(3);
			lututils.write_lut(lut,dir,"log_blue.lut");
			lut=lututils.greenyellowred_lut();
			lututils.write_lut(lut,dir,"greenyellowred.lut");
			lut=lututils.bluegreenred_lut();
			lututils.write_lut(lut,dir,"bluegreenred.lut");
			lut=lututils.angle_lut();
			lututils.write_lut(lut,dir,"angle.lut");
			lut=lututils.rand_lut();
			lututils.write_lut(lut,dir,"rand.lut");
		} else {
			byte[][] lut=null;
			switch(index){
				case 0: 	lut=lututils.nice_lut(false);
						lututils.write_lut(lut,dir,"NICE.lut"); break;
				case 1:	lut=lututils.nice_lut(true);
						lututils.write_lut(lut,dir,"NICE_whiteback.lut"); break;
				case 2:	lut=lututils.log_lut(0);
						lututils.write_lut(lut,dir,"log_gray.lut"); break;
				case 3:	lut=lututils.log_lut(1);
						lututils.write_lut(lut,dir,"log_red.lut"); break;
				case 4:	lut=lututils.log_lut(2);
						lututils.write_lut(lut,dir,"log_green.lut"); break;
				case 5:	lut=lututils.log_lut(3);
						lututils.write_lut(lut,dir,"log_blue.lut"); break;
				case 6:	lut=lututils.greenyellowred_lut();
						lututils.write_lut(lut,dir,"greenyellowred.lut"); break;
				case 7:	lut=lututils.bluegreenred_lut();
						lututils.write_lut(lut,dir,"bluegreenred.lut"); break;
				case 8:	lut=lututils.angle_lut();
						lututils.write_lut(lut,dir,"angle.lut"); break;
				case 9:	lut=lututils.rand_lut();
						lututils.write_lut(lut,dir,"rand.lut"); break;
			}
		}
	}

}
