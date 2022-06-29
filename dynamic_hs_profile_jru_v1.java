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
import jguis.*;

public class dynamic_hs_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		Roi roi=imp.getRoi();
		if(roi==null && RoiManager.getInstance()==null){
			imp.setRoi(imp.getWidth()/4,imp.getHeight()/4,imp.getWidth()/2,imp.getHeight()/2);
		}
		dynamic_profile_panel dp=new dynamic_profile_panel();
		boolean mdim=false;
		int channels=imp.getNChannels();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int currchan=imp.getChannel();
		int currslice=imp.getSlice();
		int currframe=imp.getFrame();
		if(channels>1 && slices>1) mdim=true;
		if(channels>1 && frames>1) mdim=true;
		if(slices>1 && frames>1) mdim=true;
		int sdim=0;
		if(channels==1) sdim=1;
		if(channels==1 && slices==1) sdim=2;
		String[] sdims={"Channel","Slice","Frame"};
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Profile_Dimension",sdims,sdims[sdim]);
		gd.addCheckbox("Subtract_Background?",false);
		if(channels>1){
			gd.addCheckbox("Color_Ratio?",false);
			gd.addNumericField("Numerator_Channel",1,0);
			gd.addNumericField("Demoninator_Channel",2,0);
		}
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sdim=gd.getNextChoiceIndex();
		boolean back=gd.getNextBoolean();
		//imp.saveRoi();
		Polygon backpoly=null;
		if(back){
			//imp.saveRoi();
			imp.killRoi();
			new WaitForUserDialog("Select Background with Roi").show();
			Roi backroi=imp.getRoi();
			backpoly=copypoly(backroi.getPolygon());
			//Rectangle r=backpoly.getBounds();
			//IJ.log(""+r.x+" , "+r.y);
			//imp.killRoi();
			imp.restoreRoi();
		}
		int[] nindices={currslice-1,currframe-1};
		if(sdim==1) nindices=new int[]{currchan-1,currframe-1};
		if(sdim==2) nindices=new int[]{currslice-1,currchan-1};
		boolean ratio=false;
		int[] dindices=new int[2];
		if(channels>1 && sdim!=0){
			ratio=gd.getNextBoolean();
			if(ratio){
				if(sdim==1){
					nindices[0]=(int)gd.getNextNumber()-1;
					nindices[1]=currframe-1;
					dindices[0]=(int)gd.getNextNumber()-1;
					dindices[1]=currframe-1;
				} else {
					nindices[1]=(int)gd.getNextNumber()-1;
					nindices[0]=currslice-1;
					dindices[1]=(int)gd.getNextNumber()-1;
					dindices[0]=currslice-1;
				}
			}
		}
		dp.init(imp,"Avg",backpoly,sdim,nindices,ratio,dindices);
		dynamic_profile_panel.launch_frame(dp);
	}

	public Polygon copypoly(Polygon poly){
		return new Polygon(poly.xpoints.clone(),poly.ypoints.clone(),poly.npoints);
	}

}
