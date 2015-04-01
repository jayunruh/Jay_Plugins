/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import jalgs.delimit_string;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;

public class TableDialog2 extends JDialog implements ActionListener,KeyListener,TableModelListener{
	public static TableDialog2 dialog;
	private int height,width;
	private JTable table;
	private Object[][] tabledata1;
	private static Object[][] outdata;
	private List<TableDialogListener> listeners;

	public static Object[][] showDialog(Component frameComp,Component locationComp,String title,String[] column_labels,Object[][] tabledata,int[][] options){
		Frame frame=JOptionPane.getFrameForComponent(frameComp);
		dialog=new TableDialog2(frame,locationComp,title,column_labels,tabledata,options);
		dialog.setVisible(true);
		return outdata;
	}
	
	public static Object[][] showDialog(TableDialog2 dialog){
		TableDialog2.dialog=dialog;
		dialog.setVisible(true);
		return outdata;
	}

	public TableDialog2(Frame frame,Component locationComp,String title,String[] column_labels,Object[][] tabledata,int[][] options){
		//if((frame==null)?(frame=JOptionPane.getFrameForComponent(frameComp)):frame);
		super((frame==null)?(frame=JOptionPane.getFrameForComponent(null)):frame,title,true);

		// Create and initialize the buttons.
		final JButton cancelButton=new JButton("Quit");
		cancelButton.setActionCommand("Quit");
		cancelButton.addActionListener(this);
		//
		final JButton okButton=new JButton("Continue");
		okButton.setActionCommand("Continue");
		okButton.addActionListener(this);
		getRootPane().setDefaultButton(okButton);
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
		table.addKeyListener(this);
		table.setPreferredScrollableViewportSize(new Dimension(750,80));
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		// table.setFillsViewportHeight(true);
		JScrollPane scrollpane=new JScrollPane(table,ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
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
		table.getModel().addTableModelListener(this);
		scrollpane.setPreferredSize(new Dimension(700,150));
		// scrollpane.setAlignmentX(LEFT_ALIGNMENT);
		JPanel tablepane=new JPanel();
		tablepane.setLayout(new BoxLayout(tablepane,BoxLayout.PAGE_AXIS));
		tablepane.add(Box.createRigidArea(new Dimension(0,5)));
		tablepane.add(scrollpane);
		tablepane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));

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
			int colsel=table.getSelectedColumn();
			int rowsel=table.getSelectedRow();
			TableCellEditor tce2=table.getCellEditor(rowsel,colsel);
			tce2.stopCellEditing();
			TableDialog2.outdata=new Object[height][width];
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(tabledata1[i][j] instanceof Number){
						TableDialog2.outdata[i][j]=new Double((String)table.getValueAt(i,j));
					}else{
						TableDialog2.outdata[i][j]=table.getValueAt(i,j);
					}
				}
			}
		}
		if("Quit".equals(e.getActionCommand())){
			TableDialog2.outdata=null;
		}
		TableDialog2.dialog.setVisible(false);
	}
	
	public void addTableDialogListener(TableDialogListener listener){
		if(listeners==null){
			listeners=new ArrayList<TableDialogListener>();
		}
		listeners.add(listener);
	}
	
	public void removeTableDialogListener(TableDialogListener listener){
		if(listeners==null) return;
		if(listeners.size()==0) return;
		listeners.remove(listener);
	}
	
	public void tableChanged(TableModelEvent e){
		if(listeners==null) return;
		//int colsel=table.getSelectedColumn();
		//int rowsel=table.getSelectedRow();
		//TableCellEditor tce2=table.getCellEditor(rowsel,colsel);
		//tce2.stopCellEditing();
		Object[][] tempoutdata=new Object[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(tabledata1[i][j] instanceof Number){
					tempoutdata[i][j]=new Double((String)table.getValueAt(i,j));
				}else{
					tempoutdata[i][j]=table.getValueAt(i,j);
				}
			}
		}
		for(int i=0;i<listeners.size();i++){
			(listeners.get(i)).tableDataChanged(table,e,tempoutdata);
		}
		/*int row=e.getFirstRow();
		int column=e.getColumn();
		if(column==0){
			boolean valid=((Boolean)datapanel.table.getValueAt(row,0)).booleanValue();
			if(valid!=include[row]){
				include[row]=valid;
				updateavg();
				datapanel.table.setValueAt(""+(float)intensity1[ncurves],ncurves,2);
				datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
				pwavg.updateSeries(avg,0,true);
				float[][] temperrs=new float[2][];
				temperrs[0]=errs;
				temperrs[1]=new float[errs.length];
				pwavg.addErrors(temperrs);
				if(trajectories!=null){
					pwavgtraj.updateSeries(avgtraj,0,true);
				}
			}
		}*/
	}


	public void keyTyped(KeyEvent e){
	}

	public void keyReleased(KeyEvent e){
	}

	public void keyPressed(KeyEvent e){
		if(((e.getKeyCode()==KeyEvent.VK_V)&&e.isControlDown())){
			int[] selectedRows=table.getSelectedRows();
			int[] selectedColumns=table.getSelectedColumns();
			Clipboard systemClipboard=null;
			try{
				systemClipboard=getToolkit().getSystemClipboard();
			}catch(Exception ex){
				systemClipboard=null;
			}
			if(systemClipboard==null){
				System.out.println("error opening clipboard");
				return;
			}
			int selcolumns=selectedColumns.length;
			int selrows=selectedRows.length;
			String text;
			try{
				text=(String)systemClipboard.getData(DataFlavor.stringFlavor);
			}catch(UnsupportedFlavorException ex){
				text=null;
			}catch(IOException ex){
				text=null;
			}
			if(text!=null){
				delimit_string ds=new delimit_string('\t');
				String[] lines=ds.getrows(text);
				int inlines=lines.length;
				if(selrows<inlines){
					inlines=selrows;
				}
				int incolumns=ds.getnumcolumns(lines[0]);
				if(selcolumns<incolumns){
					incolumns=selcolumns;
				}
				for(int i=0;i<inlines;i++){
					String[] linedelim=ds.delim2string(lines[i],incolumns);
					for(int j=0;j<incolumns;j++){
						table.setValueAt(linedelim[j],selectedRows[i],selectedColumns[j]);
					}
				}
			}
		}
	}
}
