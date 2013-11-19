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
import ij.plugin.frame.*;
import jalgs.jseg.*;

public class gray_morphology_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] ops={"erode","dilate","open","close","tophat_black","tophat_white"};
		gd.addChoice("Operation",ops,ops[0]);
		gd.addNumericField("Iterations",1,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		int iterations=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		float[] pixels=(float[])imp.getProcessor().convertToFloat().getPixels();
		int width=imp.getWidth();
		int height=imp.getHeight();
		float[] result=pixels.clone();
		switch(index){
			case 0:
				for(int i=0;i<iterations;i++){
					result=gray_processing.erode(result,width,height);
				}
				break;
			case 1:
				for(int i=0;i<iterations;i++){
					result=gray_processing.dilate(result,width,height);
				}
				break;
			case 2:
				result=gray_processing.open(result,iterations,width,height);
				break;
			case 3:
				result=gray_processing.close(result,iterations,width,height);
				break;
			case 4:
				result=gray_processing.tophat_black(result,iterations,width,height);
				break;
			case 5:
				result=gray_processing.tophat_white(result,iterations,width,height);
				break;
		}
		new ImagePlus(imp.getTitle()+" "+ops[index],new FloatProcessor(width,height,result,null)).show();
	}

}
