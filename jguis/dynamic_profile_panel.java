/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import jalgs.*;

import java.awt.*;
import java.awt.event.*;

import ij.*;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;

public class dynamic_profile_panel extends Panel implements ActionListener,MouseMotionListener,ImageListener{

	private Button quit_button,update_button;
	public ImagePlus imp;
	public PlotWindow4 pw;
	public String stat;
	public int[] backcoords,nindices,dindices;
	public Polygon backroi;
	public int sdim;
	boolean ratio;

	public static Frame launch_frame(dynamic_profile_panel panel){
		final Frame f=new Frame("Dynamic Profile");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				Component[] comps=f.getComponents();
				for(int i=0;i<comps.length;i++){
					comps[i].setVisible(false);
				}
				f.dispose();
			}
		});

		f.setLayout(null);
		panel.setBounds(10,40,180,150);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(200,200));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(ImagePlus imp,String stat){
		init(imp,stat,null,0,new int[]{0,0},false,null);
	}

	public void init(ImagePlus imp,String stat,Polygon back,int sdim,int[] nindices,boolean ratio,int[] dindices){
		setLayout(null);
		this.imp=imp;
		quit_button=new Button("Quit");
		quit_button.setBounds(10,10,100,30);
		quit_button.addActionListener(this);
		update_button=new Button("Update");
		update_button.setBounds(10,50,100,30);
		update_button.addActionListener(this);
		add(quit_button);
		add(update_button);
		this.stat=stat;
		backroi=back;
		if(backroi!=null)
			backcoords=get_roi_coords(backroi);
		this.sdim=sdim;
		this.nindices=nindices;
		this.dindices=dindices;
		this.ratio=ratio;
		update_plot();
		imp.getCanvas().addMouseMotionListener(this);
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void endprog(){
		if(imp!=null)
			imp.getCanvas().removeMouseMotionListener(this);
		ImagePlus.removeImageListener(this);
		this.getParent().setVisible(false);
		setVisible(false);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==quit_button){
			endprog();
			return;
		}
		if(e.getSource()==update_button){
			update_plot();
		}
	}

	public int[] get_roi_coords(Polygon poly){
		Rectangle r=poly.getBounds();
		int[] tcoords=new int[r.width*r.height];
		int counter=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					tcoords[counter]=j+i*imp.getWidth();
					counter++;
				}
			}
		}
		int[] coords=new int[counter];
		System.arraycopy(tcoords,0,coords,0,counter);
		return coords;
	}

	public Object get_roi_data(Object image,int[] coords){
		if(image instanceof float[]){
			float[] roidata=new float[coords.length];
			for(int i=0;i<coords.length;i++){
				roidata[i]=((float[])image)[coords[i]];
			}
			return roidata;
		}else if(image instanceof short[]){
			short[] roidata=new short[coords.length];
			for(int i=0;i<coords.length;i++){
				roidata[i]=((short[])image)[coords[i]];
			}
			return roidata;
		}else{
			byte[] roidata=new byte[coords.length];
			for(int i=0;i<coords.length;i++){
				roidata[i]=((byte[])image)[coords[i]];
			}
			return roidata;
		}
	}

	public float[] get_spectrum(Polygon poly,int[] indices){
		int[] roicoords=get_roi_coords(poly);
		ImageStack stack=imp.getStack();
		int channels=imp.getNChannels();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int totsize=stack.getSize();
		if(totsize==slices){
			channels=slices;
			slices=1;
		}
		if(totsize==frames){
			channels=frames;
			frames=1;
		}
		Object[] series=null;
		if(sdim==0)
			series=jutils.get3DCSeries(stack,indices[0],indices[1],frames,slices,channels);
		else if(sdim==1)
			series=jutils.get3DZSeries(stack,indices[0],indices[1],frames,slices,channels);
		else
			series=jutils.get3DTSeries(stack,indices[0],indices[1],frames,slices,channels);
		float[] spectrum=new float[series.length];
		for(int i=0;i<series.length;i++){
			Object pixels=get_roi_data(series[i],roicoords);
			spectrum[i]=jstatistics.getstatistic(stat,pixels,null);
			if(backroi!=null){
				Object back=get_roi_data(series[i],backcoords);
				float val=jstatistics.getstatistic(stat,back,null);
				// IJ.log(""+val+" , "+spectrum[i]);
				if(!Float.isNaN(val))
					spectrum[i]-=val;
			}
		}
		return spectrum;
	}

	public float[] get_ratio_spectrum(Polygon poly,int[] indices1,int[] indices2){
		float[] spectrum1=get_spectrum(poly,indices1);
		float[] spectrum2=get_spectrum(poly,indices2);
		for(int i=0;i<spectrum1.length;i++){
			if(spectrum2[i]>0.0f)
				spectrum1[i]/=spectrum2[i];
			else
				spectrum1[i]=0.0f;
		}
		return spectrum2;
	}

	public void update_plot(){
		if(imp.isLocked())
			return;
		// first get the spectrum
		RoiManager rman=RoiManager.getInstance();
		float[][] spectrum=null;
		if(rman==null||rman.getCount()==0){
			Polygon poly=imp.getRoi().getPolygon();
			spectrum=new float[1][];
			if(ratio)
				spectrum[0]=get_ratio_spectrum(poly,nindices,dindices);
			else
				spectrum[0]=get_spectrum(poly,nindices);
		}else{
			Roi[] rois=rman.getRoisAsArray();
			spectrum=new float[rois.length][];
			for(int i=0;i<rois.length;i++){
				Polygon poly=rois[i].getPolygon();
				if(ratio)
					spectrum[i]=get_ratio_spectrum(poly,nindices,dindices);
				else
					spectrum[i]=get_spectrum(poly,nindices);
			}
		}
		if(pw==null){
			pw=new PlotWindow4("Dynamic Spectrum","Position","Intensity",spectrum,null);
			ImagePlus.addImageListener(this);
			pw.draw();
		}else{
			int newlength=spectrum[0].length;
			if(newlength!=pw.getPlot().maxpts){
				/*
				 * Plot4 p4=new Plot4("Position","Intensity",spectrum,null);
				 * pw.p3=p4; pw.updatePlot();
				 */
				float[] xvals=new float[newlength];
				for(int i=0;i<xvals.length;i++)
					xvals[i]=(float)(i+1);
				for(int i=0;i<spectrum.length;i++){
					pw.updateSeries(xvals,spectrum[i],i,true);
				}
			}else{
				for(int i=0;i<spectrum.length;i++){
					pw.updateSeries(spectrum[i],i,true);
				}
			}
		}
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		update_plot();
	}

	public void mouseClicked(MouseEvent arg0){
		update_plot();
	}

	public void mouseEntered(MouseEvent arg0){
	}

	public void mouseExited(MouseEvent arg0){
	}

	public void mousePressed(MouseEvent arg0){
	}

	public void mouseReleased(MouseEvent arg0){
		update_plot();
	}

	public void imageClosed(ImagePlus arg0){
		if(arg0.equals(pw.getImagePlus())||arg0.equals(imp))
			endprog();
	}

	public void imageOpened(ImagePlus arg0){
	}

	public void imageUpdated(ImagePlus arg0){
		if(arg0.equals(imp)){
			update_plot();
		}
	}

}
