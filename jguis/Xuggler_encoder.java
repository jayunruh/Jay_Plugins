/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.ImagePlus;
import ij.process.*;
import java.awt.image.*;
import java.util.concurrent.TimeUnit;
import com.xuggle.mediatool.*;
import com.xuggle.xuggler.*;

public class Xuggler_encoder{
	int width,height,stackframes;
	private FrameInterface finterface;
	private double frameRate;
	private int nframes;

	public Xuggler_encoder(ImagePlus imp,int nframes,FrameInterface finterface,double frameRate){
		this.finterface=finterface;
		this.width=imp.getWidth();
		this.height=imp.getHeight();
		this.stackframes=nframes;
		this.frameRate=frameRate;
		this.nframes=nframes;
	}

	public void saveAsMovie(String path){
		final IMediaWriter writer=ToolFactory.makeWriter(path);
		writer.addVideoStream(0,0,IRational.make(frameRate),width,height);
		long frameTime=(long)(1000000000.0/frameRate);
		for(int i=0;i<nframes;i++){
			ColorProcessor cp=new ColorProcessor(width,height,(int[])finterface.getNextFrame());
			BufferedImage bi=new BufferedImage(width,height,BufferedImage.TYPE_3BYTE_BGR);
			bi.getGraphics().drawImage(cp.createImage(),0,0,null);
			writer.encodeVideo(0,bi,frameTime*(long)i,TimeUnit.NANOSECONDS);
			IJ.showProgress(i,nframes);
		}
		writer.close();
	}

}
