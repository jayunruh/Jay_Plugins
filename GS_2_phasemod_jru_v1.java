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

public class GS_2_phasemod_jru_v1 implements PlugIn {

	public void run(String arg) {
		float[] sdata,gdata;
		float[][] mdata,pdata;
		int[] wlist = WindowManager.getIDList();
		String[] titles = new String[wlist.length];
		for(int i=0;i<wlist.length;i++){
			ImagePlus imp = WindowManager.getImage(wlist[i]);
			titles[i]=imp.getTitle();
		}
		GenericDialog gd = new GenericDialog("Calculate Phase and Mod");
		gd.addChoice("G image",titles,titles[0]);
		gd.addChoice("S image",titles,titles[1]);
		gd.showDialog();
		String choice = gd.getNextChoice();
		ImagePlus impg = WindowManager.getImage(choice);
		choice = gd.getNextChoice();
		ImagePlus imps = WindowManager.getImage(choice);
		FloatProcessor fpg =(FloatProcessor)impg.getProcessor();
		FloatProcessor fps = (FloatProcessor)imps.getProcessor();
		sdata = (float [])fps.getPixels();
		gdata = (float [])fpg.getPixels();
		int height = impg.getHeight();
		int width = impg.getWidth();
		mdata = new float[width][height];
		pdata = new float[width][height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				mdata[i][j]=(float)Math.sqrt(sdata[j+i*width]*sdata[j+i*width]+gdata[j+i*width]*gdata[j+i*width]);
				pdata[i][j]=(float)((Math.atan(sdata[j+i*width]/gdata[j+i*width])*180.0)/(Math.PI));
				if(gdata[j+i*width]<0.0f){pdata[i][j]+=180.0f;}
				if(sdata[j+i*width]<0.0f && gdata[j+i*width]>0.0f){pdata[i][j]+=360.0f;}
			}
		}
		FloatProcessor fpm = new FloatProcessor(mdata);
		FloatProcessor fpp = new FloatProcessor(pdata);
		ImagePlus impm = new ImagePlus("Mod Image",fpm);
		impm.show();
		ImagePlus impp = new ImagePlus("Phase Image",fpp);
		impp.show();
	}

}
