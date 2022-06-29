/*******************************************************************************
 * Copyright (c) 2020 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;

import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class contrast_saver_panel extends Panel implements ActionListener,ImageListener{

	private Button quit_button,reset_button,save_button;
	public ImagePlus imp;
	public double[][] currranges;

	public static Frame launch_frame(contrast_saver_panel panel){
		final Frame f=new Frame("Contrast Saver");
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
		panel.setBounds(10,40,180,200);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(200,250));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(ImagePlus imp){
		setLayout(null);
		this.imp=imp;
		ImagePlus.addImageListener(this);
		quit_button=new Button("Quit");
		quit_button.setBounds(10,10,100,30);
		quit_button.addActionListener(this);
		save_button=new Button("Save_Contrast");
		save_button.setBounds(10,50,100,30);
		save_button.addActionListener(this);
		reset_button=new Button("Reset_Contrast");
		reset_button.setBounds(10,90,100,30);
		reset_button.addActionListener(this);
		add(quit_button);
		add(save_button);
		add(reset_button);
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void endprog(){
		ImagePlus.removeImageListener(this);
		this.getParent().setVisible(false);
		setVisible(false);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==quit_button){
			endprog();
			return;
		}
		if(e.getSource()==reset_button){
			reset_contrast();
		}
		if(e.getSource()==save_button){
			save_contrast();
		}
	}
	
	public void save_contrast(){
		int nch=imp.getNChannels();
		int currframe=imp.getFrame();
		int currslice=imp.getSlice();	
		int currch=imp.getChannel();
		this.currranges=new double[nch][2];
		for(int i=0;i<nch;i++){
			imp.setPositionWithoutUpdate(i+1,currslice,currframe);
			currranges[i][0]=imp.getProcessor().getMin();
			currranges[i][1]=imp.getProcessor().getMax();
			IJ.log("ch "+(i+1)+" range: "+currranges[i][0]+" , "+currranges[i][1]);
		}
		imp.setPosition(currch,currslice,currframe);
	}
	
	public void reset_contrast(){
		if(currranges!=null){
			int nch=imp.getNChannels();
			int currframe=imp.getFrame();
			int currslice=imp.getSlice();	
			int currch=imp.getChannel();
			for(int i=0;i<nch;i++){
				imp.setPositionWithoutUpdate(i+1,currslice,currframe);
				imp.setDisplayRange(currranges[i][0],currranges[i][1]);
			}
			//imp.reset();
			imp.setPosition(currch,currslice,currframe);
			imp.updateChannelAndDraw();
		}
	}

	public void imageClosed(ImagePlus arg0){
		if(arg0.equals(imp))
			endprog();
	}

	public void imageOpened(ImagePlus arg0){
	}

	public void imageUpdated(ImagePlus arg0){
		/*if(arg0.equals(imp)){
			update_plot();
		}*/
	}

}
