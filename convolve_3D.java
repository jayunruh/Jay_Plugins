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
import jalgs.*;
import jalgs.jfft.*;

public class convolve_3D implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			ImageWindow iw=imp.getWindow();
			if(iw instanceof ImageWindow){titles[i]=imp.getTitle();}
			else{titles[i]="NA";}
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Image 1",titles,titles[0]);
		gd.addChoice("Image 2",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		ImageStack stack=imp.getStack();
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		ImageStack stack2=imp2.getStack();
		int width=imp.getWidth(); int height=imp.getHeight(); int slices=stack.getSize();
		int[] windex=fftutils.get_best_index(width,false,19);
		int[] hindex=fftutils.get_best_index(height,false,19);
		int[] sindex=fftutils.get_best_index(slices,false,19);

		convolution3D convfunc=new convolution3D(windex[1],hindex[1],sindex[1],windex[0],hindex[0],sindex[0]);
		Object[] conv=convfunc.convolve3D(stack.getImageArray(),stack2.getImageArray());
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<slices;i++){retstack.addSlice("",conv[i]);}
		new ImagePlus("Convolution",retstack).show();
	}

}
