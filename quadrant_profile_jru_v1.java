/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.text.*;
import jguis.*;
import jalgs.jseg.*;

public class quadrant_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Ch1_thresh",1000.0f,5,15,null);
		gd.addNumericField("Ch2_thresh",1000.0f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float ch1thresh=(float)gd.getNextNumber();
		float ch2thresh=(float)gd.getNextNumber();
		ImagePlus[] imp=jutils.selectImages(false,2,new String[]{"mask","composite"});
		ImageStack mstack=imp[0].getStack();
		int width=imp[0].getWidth();
		int height=imp[0].getHeight();
		int frames=mstack.getSize();
		ImageStack cstack=imp[1].getStack();
		int channels=imp[1].getNChannels();
		findblobs3 fb=new findblobs3(width,height);
		int[][] quads=new int[frames][4];
		TextWindow tw=new TextWindow("Quadrant Counts","frame\tnone above\tch1 above\tch2 above\tboth above","",400,200);
		for(int i=0;i<frames;i++){
			Object[] substack=jutils.get3DCSeries(cstack,0,i,frames,1,channels);
			Object mask=mstack.getPixels(i+1);
			float[] objects=fb.dofindblobs(mask,0.5f);
			float[][] stats=fb.get_all_object_stats(objects,substack,"Avg");
			for(int j=0;j<stats.length;j++){
				if(stats[j][0]>ch1thresh){
					if(stats[j][1]>ch2thresh) quads[i][3]++;
					else quads[i][1]++;
				} else {
					if(stats[j][1]>ch2thresh) quads[i][2]++;
					else quads[i][0]++;
				}
			}
			tw.append(""+(i+1)+"\t"+quads[i][0]+"\t"+quads[i][1]+"\t"+quads[i][2]+"\t"+quads[i][3]+"\n");
			IJ.showProgress(i,frames);
		}
	}

}
