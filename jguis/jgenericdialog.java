/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class jgenericdialog extends JDialog implements ActionListener{
	private static jgenericdialog dialog;
	private int length;
	private Object[] tabledata1;
	private Object[] tableitems;
	private static Object[] outdata;

	public static Object[] showDialog(Component frameComp,Component locationComp,String title,String[] labels,Object[] tabledata,int[] options){
		Frame frame=JOptionPane.getFrameForComponent(frameComp);
		dialog=new jgenericdialog(frame,locationComp,title,labels,tabledata,options);
		dialog.setVisible(true);
		return outdata;
	}

	private jgenericdialog(Frame frame,Component locationComp,String title,String[] labels,Object[] tabledata,int[] options){
		super(frame,title,true);

		// Create and initialize the buttons.
		final JButton cancelButton=new JButton("Cancel");
		cancelButton.setActionCommand("Cancel");
		cancelButton.addActionListener(this);
		//
		final JButton okButton=new JButton("Continue");
		okButton.setActionCommand("Continue");
		okButton.addActionListener(this);
		getRootPane().setDefaultButton(okButton);

		tabledata1=tabledata;
		length=tabledata1.length;

		JPanel tablepane=new JPanel();
		tablepane.setLayout(new GridLayout(0,2));
		tableitems=new Object[length];
		for(int i=0;i<length;i++){
			tablepane.add(new Label(labels[i]));
			if(tabledata1[i] instanceof Number){
				tableitems[i]=new JTextField(tabledata1[i].toString(),15);
				tablepane.add((JTextField)tableitems[i]);

			}else{
				if(tabledata1[i] instanceof Boolean){
					tableitems[i]=new JCheckBox("",((Boolean)tabledata1[i]).booleanValue());
					tablepane.add((JCheckBox)tableitems[i]);
				}else{
					if(tabledata1[i] instanceof String[]){
						String[] tempitemlist=(String[])tabledata1[i];
						Object[] tempitems=new Object[tempitemlist.length];
						for(int j=0;j<tempitems.length;j++){
							tempitems[j]=tempitemlist[j];
						}
						tableitems[i]=new JComboBox(tempitems);
						tablepane.add((JComboBox)tableitems[i]);
					}else{
						tableitems[i]=new JTextField(tabledata1[i].toString());
						tablepane.add((JTextField)tableitems[i]);
					}
				}
			}

		}

		// Lay out the buttons from left to right.
		JPanel buttonPane=new JPanel();
		buttonPane.setLayout(new BoxLayout(buttonPane,BoxLayout.LINE_AXIS));
		buttonPane.setBorder(BorderFactory.createEmptyBorder(0,10,10,10));
		buttonPane.add(Box.createHorizontalGlue());
		buttonPane.add(cancelButton);
		buttonPane.add(Box.createRigidArea(new Dimension(10,0)));
		buttonPane.add(okButton);
		// Put everything together, using the content pane's BorderLayout.
		Container contentPane=getContentPane();
		contentPane.add(tablepane,BorderLayout.CENTER);
		contentPane.add(buttonPane,BorderLayout.PAGE_END);

		// Initialize values.
		pack();
		setLocationRelativeTo(locationComp);
	}

	// Handle clicks on the Continue and Cancel buttons.
	public void actionPerformed(ActionEvent e){
		if("Continue".equals(e.getActionCommand())){
			jgenericdialog.outdata=new Object[length];
			for(int i=0;i<length;i++){
				if(tabledata1[i] instanceof Number){
					jgenericdialog.outdata[i]=Double.valueOf(((JTextField)tableitems[i]).getText());
				}else{
					if(tabledata1[i] instanceof Boolean){
						jgenericdialog.outdata[i]=Boolean.valueOf(((JCheckBox)tableitems[i]).isSelected());
					}else{
						if(tabledata1[i] instanceof String[]){
							jgenericdialog.outdata[i]=Integer.valueOf(((JComboBox)tableitems[i]).getSelectedIndex());
						}else{
							jgenericdialog.outdata[i]=((JTextField)tableitems[i]).getText();
						}
					}
				}
			}
		}
		if("Cancel".equals(e.getActionCommand())){
			jgenericdialog.outdata=null;
		}
		jgenericdialog.dialog.setVisible(false);
	}
}
