/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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
//import ij.plugin.filter.*;

public class flatten_shifted_labels_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Font Size",14,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int fsize=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		RoiManager rman=RoiManager.getInstance();
		//rman.runCommand("show none");
		rman.runCommand("show all without labels");
		IJ.wait(100);
		ImagePlus imp2=imp.flatten();
		//(imp.flatten()).show();
		//IJ.selectWindow(imp.getTitle());
		IJ.wait(100);
		//rman.runCommand("show all");
		rman.runCommand("show all with labels");
		Roi[] rois=rman.getRoisAsArray();
		//Filler filler=new Filler();
		ImageProcessor ip2=imp2.getProcessor();
		for(int i=0;i<rois.length;i++){
			//Rectangle r=new Rectangle(rois[i].getBounds().x-10,rois[i].getBounds().y-10,20,20);
			//filler.drawLabel(imp2,imp2.getProcessor(),i+1,r);
			drawLabel(ip2,i+1,rois[i].getBounds().x,rois[i].getBounds().y,fsize);
		}
		imp2.show();
	}

	public void drawLabel(ImageProcessor ip,int count,int x, int y,int size){
		ip.setFont(new Font("SansSerif",Font.PLAIN,size));
		String label=""+count;
		int w=ip.getStringWidth(label);
		int x1=x-w/2;
		int y1=y+size/2;
		FontMetrics metrics=ip.getFontMetrics();
		int h=metrics.getHeight();
		ip.setColor(Color.black);
		ip.setRoi(x1-1,y1-h+2,w+1,h-3);
		ip.fill();
		ip.resetRoi();
		ip.setColor(Color.white);
		ip.drawString(label,x1,y1);
		return;
	}

}
