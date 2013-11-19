/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import java.awt.Dimension;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.JTable;
import javax.swing.table.*;

public class TablePanel extends JPanel{
	private int height,width;
	public JTable table;
	private Object[][] tabledata1;

	public void init(String[] column_labels,Object[][] tabledata,int[][] options){
		height=tabledata.length;
		width=tabledata[0].length;
		tabledata1=tabledata;
		Object[][] celldata=new Object[height][width];
		// data are arranged [row][column]
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(tabledata1[i][j]==null){
					tabledata1[i][j]=new String("");
				}
				if(tabledata1[i][j] instanceof String[]){
					celldata[i][j]=((String[])tabledata1[i][j])[options[i][j]];
				}else{
					if(tabledata1[i][j] instanceof Number){
						celldata[i][j]=tabledata1[i][j].toString();
					}else{
						celldata[i][j]=tabledata1[i][j];
					}
				}
			}
		}
		DefaultTableModel dm=new DefaultTableModel();
		dm.setDataVector(celldata,column_labels);
		table=new JTable(dm);
		table.setPreferredScrollableViewportSize(new Dimension(500,100));
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		JScrollPane scrollpane=new JScrollPane(table,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		custom_cell_editor2 tce=new custom_cell_editor2(tabledata1);
		// note that columns must have unique names for the following to work
		for(int i=0;i<width;i++){
			table.getColumn(table.getColumnName(i)).setCellRenderer(new custom_cell_renderer());
			table.getColumn(table.getColumnName(i)).setCellEditor(tce);
			table.getColumn(table.getColumnName(i)).setPreferredWidth(75);
		}
		table.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		table.setCellSelectionEnabled(true);
		table.changeSelection(0,1,false,false);
		scrollpane.setPreferredSize(new Dimension(100,100));
		setLayout(new BoxLayout(this,BoxLayout.PAGE_AXIS));
		add(Box.createRigidArea(new Dimension(0,5)));
		add(scrollpane);
	}

}
