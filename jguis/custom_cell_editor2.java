/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import java.awt.Component;
import java.util.EventObject;

import javax.swing.DefaultCellEditor;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.CellEditorListener;
import javax.swing.table.TableCellEditor;

public class custom_cell_editor2 implements TableCellEditor{
	protected Object[][] editors;
	int width,height;
	int row,column;

	public custom_cell_editor2(Object[][] tabledata){
		height=tabledata.length;
		width=tabledata[0].length;
		editors=new Object[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(tabledata[i][j] instanceof String){
					editors[i][j]=new DefaultCellEditor(new JTextField());
				}else{
					if(tabledata[i][j] instanceof Number){
						editors[i][j]=new DefaultCellEditor(new JTextField());
					}else{
						if(tabledata[i][j] instanceof Boolean){
							editors[i][j]=new DefaultCellEditor(new JCheckBox());
						}else{
							if(tabledata[i][j] instanceof String[]){
								JComboBox temp=new JComboBox();
								for(int k=0;k<((String[])tabledata[i][j]).length;k++){
									temp.addItem(((String[])tabledata[i][j])[k]);
								}
								editors[i][j]=new DefaultCellEditor(temp);
							}
						}
					}
				}
			}
		}
	}

	public Component getTableCellEditorComponent(JTable table,Object value,boolean isSelected,int row1,int column1){
		// currenteditor=(TableCellEditor)editors.get(new Integer(row));
		// currenteditor=checkboxeditor;
		row=row1;
		column=column1;
		return ((TableCellEditor)editors[row][column]).getTableCellEditorComponent(table,value,isSelected,row,column);
	}

	public Object getCellEditorValue(){
		return ((TableCellEditor)editors[row][column]).getCellEditorValue();
	}

	public boolean stopCellEditing(){
		return ((TableCellEditor)editors[row][column]).stopCellEditing();
	}

	public void cancelCellEditing(){
		((TableCellEditor)editors[row][column]).cancelCellEditing();
	}

	public boolean isCellEditable(EventObject anEvent){
		return ((TableCellEditor)editors[row][column]).isCellEditable(anEvent);
	}

	public void addCellEditorListener(CellEditorListener l){
		((TableCellEditor)editors[row][column]).addCellEditorListener(l);
	}

	public void removeCellEditorListener(CellEditorListener l){
		((TableCellEditor)editors[row][column]).removeCellEditorListener(l);
	}

	public boolean shouldSelectCell(EventObject anEvent){
		return ((TableCellEditor)editors[row][column]).shouldSelectCell(anEvent);
	}

}
