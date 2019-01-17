/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.ImagePlus;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.text.TextPanel;
import ij.text.TextWindow;

import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class IJTableListener extends Panel implements ActionListener,MouseListener{

	private Button quit_button,update_button;
	public ImagePlus imp;
	public TextWindow tw;
	public int xcol,ycol,zcol;
	public float zratio,psize;

	public static Frame launch_frame(IJTableListener panel){
		final Frame f=new Frame("Table Listener");
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

	public void init(TextWindow tw,ImagePlus imp,int xcol,int ycol,int zcol,float zratio){
		setLayout(null);
		this.imp=imp;
		this.tw=tw;
		quit_button=new Button("Quit");
		quit_button.setBounds(10,10,100,30);
		quit_button.addActionListener(this);
		update_button=new Button("Update");
		update_button.setBounds(10,50,100,30);
		update_button.addActionListener(this);
		add(quit_button);
		add(update_button);
		tw.getTextPanel().addMouseListener(this);
		this.xcol=xcol; this.ycol=ycol; this.zcol=zcol; this.zratio=zratio; this.psize=1.0f;
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void endprog(){
		tw.getTextPanel().removeMouseListener(this);
		this.getParent().setVisible(false);
		setVisible(false);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==quit_button){
			endprog();
			return;
		}
		if(e.getSource()==update_button){
			update_table();
		}
	}
	
	public void update_table(){
		TextPanel tp=tw.getTextPanel();
		int selStart=tp.getSelectionStart();
		String line=tp.getLine(selStart);
		String[] line2=table_tools.split_string_tab(line);
		int x=(int)(Float.parseFloat(line2[xcol])/psize);
		int y=(int)(Float.parseFloat(line2[ycol])/psize);
		int z=(int)((Float.parseFloat(line2[zcol])/zratio)/psize);
		Roi roi=new PointRoi(x,y);
		roi.setPosition(1,z,1);
		imp.setRoi(roi);
		int currchan=imp.getChannel();
		int currframe=imp.getFrame();
		imp.setPosition(currchan,z,currframe);
		//imp.setSlice(z);
		imp.updateAndDraw();
	}

	public void mouseClicked(MouseEvent arg0){
		update_table();
	}

	public void mouseEntered(MouseEvent arg0){
	}

	public void mouseExited(MouseEvent arg0){
	}

	public void mousePressed(MouseEvent arg0){
		update_table();
	}

	public void mouseReleased(MouseEvent arg0){
	}

}
