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
import javax.swing.event.*;
import ij.*;
import java.io.*;

import jhcsfcs.*;

public class fcs_results_panel extends Panel implements ActionListener,ItemListener,ListSelectionListener,TableModelListener{
	private TablePanel platetable;
	private TablePanel welltable;
	private TablePanel celltable;
	private PlotStack4 corrplots;
	private PlotStack4 trajplots;
	private Button savebutton,filterbutton,recalcbutton,allvalidbutton,getcellavgbutton,getwellavgbutton;
	private Choice wellstatchoice,cellstatchoice;
	private Label wellstatlabel,cellstatlabel;
	private int currwell,currcell,currrun,maxcells,maxruns,nparams,wellstat,cellstat;
	private hcs_fcs_plate plate;
	public String[] paramsnames;
	public final static int H=775;
	public final static int WR=650;
	public Dimension totalSize=new Dimension();
	private Object[][] filter;

	public static void launch_frame(fcs_results_panel panel){
		final Frame f=new Frame("FCS Results Panel");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				f.dispose();
			}
		});

		f.setLayout(null);
		panel.setBounds(10,40,750,700);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(800,800));
		f.setVisible(true);
		panel.requestFocus();
	}

	public void init(String filename){
		// here we initialize the panel from a saved file
		try{
			InputStream os=new BufferedInputStream(new FileInputStream(filename));
			jdataio jdio=new jdataio();
			int nfiles=jdio.readintelint(os);
			Plot4[] corrplots2=new Plot4[nfiles];
			Plot4[] trajplots2=new Plot4[nfiles];
			int wellstat1=jdio.readintelint(os);
			int cellstat1=jdio.readintelint(os);
			int nparams1=jdio.readintelint(os);
			String[] paramsnames1=new String[nparams1];
			for(int i=0;i<nparams1;i++){
				paramsnames1[i]=jdio.readstring(os);
			}
			String tempdir=jdio.readstring(os);
			int nwells=jdio.readintelint(os); // number of wells
			hcs_fcs_well[] wells=new hcs_fcs_well[nwells];
			hcs_fcs_plate plate1=new hcs_fcs_plate(wells,tempdir);
			int id=0;
			for(int i=0;i<nwells;i++){
				String fname=jdio.readstring(os);
				int ncells=jdio.readintelint(os); // number of cells
				hcs_fcs_cell[] cells=new hcs_fcs_cell[ncells];
				wells[i]=new hcs_fcs_well(cells,null,fname,0,plate1);
				for(int j=0;j<ncells;j++){
					String cellname=jdio.readstring(os);
					int nruns=jdio.readintelint(os); // number of runs
					hcs_fcs_run[] runs=new hcs_fcs_run[nruns];
					cells[j]=new hcs_fcs_cell(runs,null,0,0,0,0,cellname,wells[i]);
					for(int k=0;k<nruns;k++){
						int runnumber=jdio.readintelint(os);
						double[] params=new double[nparams1];
						jdio.readintelfloatfile(os,nparams1,params);
						corrplots2[id]=new Plot4(os);
						trajplots2[id]=new Plot4(os);
						boolean valid1=jdio.readboolean(os);
						int nfnames=jdio.readintelint(os);
						String[] rfnames=new String[nfnames];
						for(int l=0;l<nfnames;l++){
							rfnames[l]=jdio.readstring(os);
						}
						runs[k]=new hcs_fcs_run(rfnames,runnumber,cells[j]);
						runs[k].id=id;
						runs[k].params=params;
						runs[k].valid=valid1;
						id++;
					}
				}
			}
			PlotStack4 corrplots1=new PlotStack4("FCS Batch Correlation",corrplots2);
			corrplots1.draw();
			PlotStack4 trajplots1=new PlotStack4("FCS Batch Trajectories",trajplots2);
			trajplots1.draw();
			os.close();
			init(plate1,paramsnames1,corrplots1,trajplots1,wellstat1,cellstat1);

		}catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
	}

	public void init(hcs_fcs_plate plate,String[] paramsnames,PlotStack4 corrplots,PlotStack4 trajplots){
		init(plate,paramsnames,corrplots,trajplots,0,0);
	}

	public void init(hcs_fcs_plate plate,String[] paramsnames,PlotStack4 corrplots,PlotStack4 trajplots,int wellstat1,int cellstat1){
		setLayout(null);
		this.plate=plate;
		this.paramsnames=paramsnames;
		this.corrplots=corrplots;
		this.trajplots=trajplots;
		currwell=0;
		currcell=0;
		currrun=0;
		nparams=paramsnames.length;
		plate.generate_ids();
		maxcells=plate.wells[0].cells.length;
		maxruns=plate.wells[0].cells[0].runs.length;
		hcs_fcs_well[] wells=plate.wells;
		wellstat=wellstat1;
		cellstat=cellstat1;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			if(cells.length>maxcells){
				maxcells=cells.length;
			}
			for(int j=0;j<cells.length;j++){
				if(cells[j].runs.length>maxruns){
					maxruns=cells[j].runs.length;
				}
			}
		}
		update_params();
		String[] col_labels=new String[nparams+1];
		col_labels[0]="Name";
		String[] runcol_labels=new String[nparams+2];
		runcol_labels[0]="Name";
		runcol_labels[1]="Valid";
		for(int i=0;i<nparams;i++){
			col_labels[i+1]=paramsnames[i];
			runcol_labels[i+2]=paramsnames[i];
		}
		Object[][] plate_data=new Object[plate.wells.length][nparams+1];
		for(int i=0;i<plate.wells.length;i++){
			plate_data[i][0]=plate.wells[i].foldername;
			for(int j=0;j<nparams;j++){
				plate_data[i][j+1]=""+(float)plate.wells[i].params[j];
			}
		}

		platetable=new TablePanel();
		platetable.init(col_labels,plate_data,null);
		platetable.setBounds(10,30,730,300);
		platetable.table.getSelectionModel().addListSelectionListener(this);
		welltable=new TablePanel();
		Object[][] well_data=new Object[maxcells][nparams+1];
		for(int i=0;i<plate.wells[0].cells.length;i++){
			well_data[i][0]=plate.wells[0].cells[i].cellname;
			for(int j=0;j<nparams;j++){
				well_data[i][j+1]=""+(float)plate.wells[0].cells[i].params[j];
			}
		}
		for(int i=plate.wells[0].cells.length;i<maxcells;i++){
			for(int j=0;j<nparams+1;j++){
				well_data[i][j]="";
			}
		}
		welltable.init(col_labels,well_data,null);
		welltable.setBounds(10,340,730,200);
		welltable.table.getSelectionModel().addListSelectionListener(this);
		celltable=new TablePanel();
		Object[][] cell_data=new Object[maxruns][nparams+2];
		for(int i=0;i<plate.wells[0].cells[0].runs.length;i++){
			cell_data[i][0]=new Integer(plate.wells[0].cells[0].runs[i].runnumber);
			cell_data[i][1]=new Boolean(plate.wells[0].cells[0].runs[i].valid);
			for(int j=0;j<nparams;j++){
				cell_data[i][j+2]=""+(float)plate.wells[0].cells[0].runs[i].params[j];
			}
		}
		for(int i=plate.wells[0].cells[0].runs.length;i<maxruns;i++){
			cell_data[i][0]=new String();
			cell_data[i][1]=new Boolean(false);
			for(int j=0;j<nparams;j++){
				cell_data[i][j+2]="";
			}
		}
		celltable.init(runcol_labels,cell_data,null);
		celltable.setBounds(10,550,730,110);
		celltable.table.getSelectionModel().addListSelectionListener(this);
		celltable.table.getModel().addTableModelListener(this);
		int tempx=5;
		savebutton=new Button("Save");
		savebutton.setBounds(10,tempx,40,20);
		savebutton.addActionListener(this);
		filterbutton=new Button("Filter");
		filterbutton.setBounds(60,tempx,40,20);
		filterbutton.addActionListener(this);
		wellstatlabel=new Label("Well Stat");
		wellstatlabel.setBounds(110,tempx,60,20);
		add(wellstatlabel);
		wellstatchoice=new Choice();
		for(int i=0;i<jstatistics.stats.length;i++){
			wellstatchoice.add(jstatistics.stats[i]);
		}
		wellstatchoice.setBounds(180,tempx,60,20);
		wellstatchoice.addItemListener(this);
		add(wellstatchoice);
		cellstatlabel=new Label("Cell Stat");
		cellstatlabel.setBounds(250,tempx,60,20);
		add(cellstatlabel);
		cellstatchoice=new Choice();
		for(int i=0;i<jstatistics.stats.length;i++){
			cellstatchoice.add(jstatistics.stats[i]);
		}
		cellstatchoice.setBounds(320,tempx,60,20);
		cellstatchoice.addItemListener(this);
		recalcbutton=new Button("ReCalc");
		recalcbutton.setBounds(390,tempx,50,20);
		recalcbutton.addActionListener(this);
		allvalidbutton=new Button("All Valid");
		allvalidbutton.setBounds(450,tempx,60,20);
		allvalidbutton.addActionListener(this);
		getcellavgbutton=new Button("Cell Avg Corr");
		getcellavgbutton.setBounds(520,tempx,80,20);
		getcellavgbutton.addActionListener(this);
		getwellavgbutton=new Button("Well Avg Corr");
		getwellavgbutton.setBounds(610,tempx,80,20);
		getwellavgbutton.addActionListener(this);
		add(cellstatchoice);
		add(platetable);
		add(welltable);
		add(celltable);
		add(savebutton);
		add(filterbutton);
		add(recalcbutton);
		add(allvalidbutton);
		add(getcellavgbutton);
		add(getwellavgbutton);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==savebutton){
			saveAsObject();
		}
		if(e.getSource()==filterbutton){
			run_filter();
		}
		if(e.getSource()==recalcbutton){
			update_params();
			updatetables();
		}
		if(e.getSource()==allvalidbutton){
			validate_all();
			update_params();
			updatetables();
		}
		if(e.getSource()==getcellavgbutton){
			Object[] temp=getcellavgcorr(plate.wells[currwell].cells[currcell]);
			float[][] cellavg=(float[][])temp[1];
			float[][] cellavgxvals=(float[][])temp[0];
			PlotWindow4 pw=new PlotWindow4(plate.wells[currwell].cells[currcell].cellname+" Avg","G(tau)","tau (s)",cellavgxvals,cellavg,null);
			pw.draw();
			pw.setLogAxes(true,false);
		}
		if(e.getSource()==getwellavgbutton){
			Object[] temp=getwellavgcorr(plate.wells[currwell]);
			float[][] wellavg=(float[][])temp[1];
			float[][] wellavgxvals=(float[][])temp[0];
			PlotWindow4 pw=new PlotWindow4(plate.wells[currwell].foldername+" Avg","G(tau)","tau (s)",wellavgxvals,wellavg,null);
			pw.draw();
			pw.setLogAxes(true,false);
		}
	}

	public Object[] getcellavgcorr(hcs_fcs_cell cell){
		hcs_fcs_run[] temp=cell.runs;
		int corrsize=corrplots.p3[temp[0].id].getmaxpts();
		int ncorrs=corrplots.p3[temp[0].id].getNSeries();
		float[][] cellavg=new float[ncorrs][corrsize];
		float[][] xvals=corrplots.p3[temp[0].id].getXValues();
		float[][] cellavgxvals=xvals.clone();
		int nvalid=0;
		for(int i=0;i<temp.length;i++){
			if(temp[i].valid){
				float[][] yvals=corrplots.p3[temp[i].id].getYValues();
				for(int j=0;j<ncorrs;j++){
					for(int k=0;k<corrsize;k++){
						cellavg[j][k]+=yvals[j][k];
					}
				}
				nvalid++;
			}
		}
		for(int j=0;j<ncorrs;j++){
			for(int k=0;k<corrsize;k++){
				cellavg[j][k]/=(float)nvalid;
			}
		}
		Object[] corrs={cellavgxvals,cellavg,new Integer(nvalid)};
		return corrs;
	}

	public Object[] getwellavgcorr(hcs_fcs_well well){
		hcs_fcs_cell[] temp=well.cells;
		int corrsize=corrplots.p3[temp[0].runs[0].id].getmaxpts();
		int ncorrs=corrplots.p3[temp[0].runs[0].id].getNSeries();
		float[][] cellavg=new float[ncorrs][corrsize];
		float[][] xvals=corrplots.p3[temp[0].runs[0].id].getXValues();
		float[][] cellavgxvals=xvals.clone();
		int wellnvalid=0;
		for(int i=0;i<temp.length;i++){
			Object[] temp2=getcellavgcorr(temp[i]);
			int nvalid=((Integer)temp2[2]).intValue();
			if(nvalid>0){
				float[][] yvals=(float[][])temp2[1];
				for(int j=0;j<ncorrs;j++){
					for(int k=0;k<corrsize;k++){
						cellavg[j][k]+=(float)nvalid*yvals[j][k];
					}
				}
				wellnvalid+=nvalid;
			}
		}
		for(int j=0;j<ncorrs;j++){
			for(int k=0;k<corrsize;k++){
				cellavg[j][k]/=(float)wellnvalid;
			}
		}
		Object[] corrs={cellavgxvals,cellavg,new Integer(wellnvalid)};
		return corrs;
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==wellstatchoice){
			wellstat=wellstatchoice.getSelectedIndex();
			update_params();
			updatetables();
		}
		if(e.getSource()==cellstatchoice){
			cellstat=cellstatchoice.getSelectedIndex();
			update_params();
			updatetables();
		}
	}

	public void valueChanged(ListSelectionEvent e){
		if(e.getSource()==platetable.table.getSelectionModel()){
			int selrow=platetable.table.getSelectedRow();
			setcurrwell(selrow);
		}
		if(e.getSource()==welltable.table.getSelectionModel()){
			int selrow=welltable.table.getSelectedRow();
			setcurrcell(selrow);
		}
		if(e.getSource()==celltable.table.getSelectionModel()){
			int selrow=celltable.table.getSelectedRow();
			setcurrrun(selrow);
		}
	}

	public void tableChanged(TableModelEvent e){
		int row=e.getFirstRow();
		int column=e.getColumn();
		if(column==1){
			boolean valid=((Boolean)celltable.table.getValueAt(row,column)).booleanValue();
			plate.wells[currwell].cells[currcell].runs[row].valid=valid;
		}
	}

	private void setcurrwell(int currwell){
		int nwells=plate.wells.length;
		if(currwell>=nwells){
			currwell=nwells-1;
		}
		this.currwell=currwell;
		updatewelltable();
		setcurrcell(currcell);
	}

	private void setcurrcell(int currcell){
		int ncells=plate.wells[currwell].cells.length;
		if(currcell>=ncells){
			currcell=ncells-1;
		}
		this.currcell=currcell;
		updatecelltable();
		setcurrrun(currrun);
	}

	private void setcurrrun(int currrun){
		int nruns=plate.wells[currwell].cells[currcell].runs.length;
		if(currrun>=nruns){
			currrun=nruns-1;
		}
		this.currrun=currrun;
		hcs_fcs_run runs[]=plate.wells[currwell].cells[currcell].runs;
		int id=runs[currrun].id;
		corrplots.setSlice(id);
		trajplots.setSlice(id);
	}

	private void updateplatetable(){
		hcs_fcs_well[] wells=plate.wells;
		int nwells=wells.length;
		for(int i=0;i<nwells;i++){
			platetable.table.setValueAt(wells[i].foldername,i,0);
			for(int j=0;j<nparams;j++){
				platetable.table.setValueAt(""+(float)wells[i].params[j],i,j+1);
			}
		}
	}

	private void updatewelltable(){
		// now update the welltable
		hcs_fcs_cell[] cells=plate.wells[currwell].cells;
		int ncells=cells.length;
		for(int i=0;i<ncells;i++){
			welltable.table.setValueAt(cells[i].cellname,i,0);
			for(int j=0;j<nparams;j++){
				welltable.table.setValueAt(""+(float)cells[i].params[j],i,j+1);
			}
		}
		for(int i=ncells;i<maxcells;i++){
			welltable.table.setValueAt("",i,0);
			for(int j=0;j<nparams;j++){
				welltable.table.setValueAt("",i,j+1);
			}
		}
	}

	private void updatecelltable(){
		hcs_fcs_run[] runs=plate.wells[currwell].cells[currcell].runs;
		int nruns=runs.length;
		// now update the cell table
		for(int i=0;i<nruns;i++){
			celltable.table.setValueAt(""+runs[i].runnumber,i,0);
			celltable.table.setValueAt(new Boolean(runs[i].valid),i,1);
			for(int j=0;j<nparams;j++){
				celltable.table.setValueAt(""+(float)runs[i].params[j],i,j+2);
			}
		}
		for(int i=nruns;i<maxruns;i++){
			celltable.table.setValueAt("",i,0);
			celltable.table.setValueAt(new Boolean(false),i,1);
			for(int j=0;j<nparams;j++){
				celltable.table.setValueAt("",i,j+2);
			}
		}
	}

	private void updatetables(){
		updatecelltable();
		updatewelltable();
		updateplatetable();
	}

	void saveAsObject(){
		Frame parent=(Frame)this.getParent();
		FileDialog fd=new FileDialog(parent,"Save as Results Object...",FileDialog.SAVE);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		fd.dispose();
		save_plate(directory+File.separator+name);
	}

	private void save_plate(String outfile){
		// here we save everything including the plots to one file
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(outfile));
			jdataio jdio=new jdataio();
			jdio.writeintelint(os,plate.getnfiles());
			jdio.writeintelint(os,wellstat);
			jdio.writeintelint(os,cellstat);
			jdio.writeintelint(os,nparams);
			for(int i=0;i<nparams;i++){
				jdio.writestring(os,paramsnames[i]);
			}
			jdio.writestring(os,plate.dir);
			hcs_fcs_well[] wells=plate.wells;
			jdio.writeintelint(os,wells.length); // number of wells
			for(int i=0;i<wells.length;i++){
				jdio.writestring(os,wells[i].foldername);
				hcs_fcs_cell[] cells=wells[i].cells;
				jdio.writeintelint(os,cells.length); // number of cells
				for(int j=0;j<cells.length;j++){
					jdio.writestring(os,cells[j].cellname);
					hcs_fcs_run[] runs=cells[j].runs;
					jdio.writeintelint(os,runs.length); // number of runs
					for(int k=0;k<runs.length;k++){
						jdio.writeintelint(os,runs[k].runnumber);
						jdio.writeintelfloatarray(os,runs[k].params);
						int id=runs[k].id;
						corrplots.p3[id].saveplot2os(os);
						trajplots.p3[id].saveplot2os(os);
						jdio.writeboolean(os,runs[k].valid);
						jdio.writeintelint(os,runs[k].filenames.length);
						for(int l=0;l<runs[k].filenames.length;l++){
							jdio.writestring(os,runs[k].filenames[l]);
						}
					}
				}
			}
			os.close();
		}catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
	}

	private void run_filter(){
		boolean[] oldfiltervalid=new boolean[8];
		if(filter!=null){
			for(int j=0;j<8;j++){
				oldfiltervalid[j]=((Boolean)filter[j][0]).booleanValue();
			}
		}
		if(!filterdialog()){
			return;
		}
		boolean[] filtervalid=new boolean[8];
		int[] filterparam=new int[8];
		int[] filtertype=new int[8];
		float[] filtervalue=new float[8];
		for(int j=0;j<8;j++){
			filtervalid[j]=((Boolean)filter[j][0]).booleanValue();
			filterparam[j]=((Integer)filter[j][1]).intValue();
			filtertype[j]=((Integer)filter[j][2]).intValue();
			filtervalue[j]=((Double)filter[j][3]).floatValue();
		}
		// loop over all the filters
		for(int j=0;j<8;j++){
			if(filtervalid[j]){
				if(!oldfiltervalid[j]){
					// we have a brand new filter
					plate.filter_data(filterparam[j],filtertype[j],filtervalue[j]);
				}
			}else{
				if(oldfiltervalid[j]){
					// we undid an old filter
					plate.unfilter_data(filterparam[j],filtertype[j],filtervalue[j]);
				}
			}
		}
		update_params();
		updatetables();
	}

	private void validate_all(){
		hcs_fcs_well[] wells=plate.wells;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				hcs_fcs_run[] runs=cells[j].runs;
				for(int k=0;k<runs.length;k++){
					runs[k].valid=true;
				}
			}
		}
	}

	private boolean filterdialog(){
		String[] columnlabels={"Use Filter?","Filter Parameter","Filter Type","Filter Value"};
		if(filter==null){
			filter=new Object[8][4];
			for(int i=0;i<8;i++){
				filter[i][0]=new Boolean(false);
				filter[i][1]=new Integer(0);
				filter[i][2]=new Integer(0);
				filter[i][3]=new Double(0.0);
			}
		}
		Object[][] tabledata=new Object[8][4];
		int[][] options=new int[8][4];
		String[] types={">","<"};
		for(int i=0;i<8;i++){
			tabledata[i][0]=filter[i][0];
			tabledata[i][1]=paramsnames;
			tabledata[i][2]=types;
			tabledata[i][3]=""+((Double)filter[i][3]).floatValue();
			options[i][1]=((Integer)filter[i][1]).intValue();
			options[i][2]=((Integer)filter[i][2]).intValue();
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Filters",columnlabels,tabledata,options);
		if(retvals==null){
			return false;
		}
		for(int i=0;i<8;i++){
			filter[i][0]=retvals[i][0];
			String retparam=(String)retvals[i][1];
			for(int j=0;j<paramsnames.length;j++){
				if(paramsnames[j].equals(retparam)){
					filter[i][1]=new Integer(j);
					break;
				}
			}
			String rettype=(String)retvals[i][2];
			for(int j=0;j<types.length;j++){
				if(types[j].equals(rettype)){
					filter[i][2]=new Integer(j);
					break;
				}
			}
			filter[i][3]=new Double((String)retvals[i][3]);
		}
		return true;
	}

	private void update_params(){
		hcs_fcs_well[] wells=plate.wells;
		for(int i=0;i<wells.length;i++){
			hcs_fcs_cell[] cells=wells[i].cells;
			for(int j=0;j<cells.length;j++){
				cells[j].update_params(cellstat);
			}
			wells[i].update_params(wellstat);
		}
	}

}
