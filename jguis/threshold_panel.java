/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import jalgs.*;
import jalgs.jseg.*;

import java.awt.*;
import java.awt.event.*;
import ij.*;
import ij.gui.*;
import ij.plugin.frame.RoiManager;
import ij.process.*;
import ij.text.TextWindow;
import ij.io.*;

public class threshold_panel extends Panel implements ActionListener,ItemListener,MouseMotionListener{

	private Button combine_objects,separate_objects,add_object,delete_object,undo_button;
	private Button enlarge_object,edit_object,launch_sky,dilate_objects,erode_objects,fill_holes;
	private Button save_objects,log_number,obj_stats;
	private Checkbox mask_check,rank_check;
	private Label idlabel,nobjects,arealabel;
	public ImagePlus imp;
	public float[] objects,oldobjects;
	private int[] areas;
	private int[] arearank;
	public findblobs3 fb;
	public boolean showmask,showrank;

	public static Frame launch_frame(threshold_panel panel){
		final Frame f=new Frame("Object Selection");
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
		panel.setBounds(10,40,180,560);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(200,600));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(ImagePlus imp,float[] objects,boolean showmask){
		setLayout(null);
		showrank=false;
		this.showmask=showmask;
		this.imp=imp;
		fb=new findblobs3(imp.getWidth(),imp.getHeight());
		this.objects=objects;
		combine_objects=new Button("Combine Objects");
		combine_objects.setBounds(10,10,100,30);
		combine_objects.addActionListener(this);
		add(combine_objects);
		separate_objects=new Button("Separate Objects");
		separate_objects.setBounds(10,40,100,30);
		separate_objects.addActionListener(this);
		add(separate_objects);
		add_object=new Button("Add Object");
		add_object.setBounds(10,70,100,30);
		add_object.addActionListener(this);
		add(add_object);
		delete_object=new Button("Delete Object");
		delete_object.setBounds(10,100,100,30);
		delete_object.addActionListener(this);
		add(delete_object);
		enlarge_object=new Button("Enlarge Object");
		enlarge_object.setBounds(10,130,100,30);
		enlarge_object.addActionListener(this);
		add(enlarge_object);
		edit_object=new Button("Edit Object");
		edit_object.setBounds(10,160,100,30);
		edit_object.addActionListener(this);
		add(edit_object);
		mask_check=new Checkbox("Show Mask",showmask);
		mask_check.setBounds(10,190,100,20);
		mask_check.addItemListener(this);
		add(mask_check);
		rank_check=new Checkbox("Show Rank",showrank);
		rank_check.setBounds(10,210,100,20);
		rank_check.addItemListener(this);
		add(rank_check);
		idlabel=new Label("id = n/a");
		idlabel.setBounds(10,230,100,20);
		add(idlabel);
		nobjects=new Label("# of objects = n/a");
		nobjects.setBounds(10,250,100,20);
		add(nobjects);
		arealabel=new Label("object area = n/a");
		arealabel.setBounds(10,270,100,20);
		add(arealabel);
		launch_sky=new Button("Launch SKY");
		launch_sky.setBounds(10,290,100,30);
		launch_sky.addActionListener(this);
		add(launch_sky);
		dilate_objects=new Button("Dilate Objects");
		dilate_objects.setBounds(10,320,100,30);
		dilate_objects.addActionListener(this);
		add(dilate_objects);
		erode_objects=new Button("Erode Objects");
		erode_objects.setBounds(10,350,100,30);
		erode_objects.addActionListener(this);
		add(erode_objects);
		fill_holes=new Button("Fill Holes");
		fill_holes.setBounds(10,380,100,30);
		fill_holes.addActionListener(this);
		add(fill_holes);
		save_objects=new Button("Save");
		save_objects.setBounds(10,410,100,30);
		save_objects.addActionListener(this);
		add(save_objects);
		log_number=new Button("Log");
		log_number.setBounds(10,440,100,30);
		log_number.addActionListener(this);
		add(log_number);
		obj_stats=new Button("Obj Stats");
		obj_stats.setBounds(10,470,100,30);
		obj_stats.addActionListener(this);
		add(obj_stats);
		undo_button=new Button("Undo");
		undo_button.setBounds(10,500,100,30);
		undo_button.addActionListener(this);
		add(undo_button);
		update_image();
		imp.getCanvas().addMouseMotionListener(this);
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==undo_button){
			objects=oldobjects.clone();
			fb.set_objects(objects);
			update_image();
			return;
		}
		oldobjects=objects.clone();
		Roi roi=imp.getRoi();
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=new Roi[0];
		if(rman!=null){
			rois=rman.getRoisAsArray();
		}
		if(roi!=null||rois.length>0){
			if(e.getSource()==combine_objects){
				Polygon poly=roi.getPolygon();
				fb.combine_objects(objects,poly);
			}
			if(e.getSource()==separate_objects){
				if(rois.length>0){
					for(int i=0;i<rois.length;i++){
						separate_objects(rois[i]);
					}
				}else{
					separate_objects(roi);
				}
			}
			if(e.getSource()==add_object){
				Polygon poly=roi.getPolygon();
				fb.add_object(objects,poly);
			}
			if(e.getSource()==delete_object){
				Rectangle rect=roi.getBounds();
				float id=objects[rect.x+rect.y*imp.getWidth()];
				if(id>0.0f){
					fb.delete_object(objects,id);
				}
			}
			if(e.getSource()==enlarge_object){
				Polygon poly=roi.getPolygon();
				fb.combine_objects(objects,poly);
			}
			if(e.getSource()==edit_object){
				Polygon poly=roi.getPolygon();
				fb.edit_object(objects,poly);
			}
		}
		if(e.getSource()==launch_sky){
			String[] labels={"SKY_Image","Spectral_Image","Spectra"};
			ImagePlus[] imps=jutils.selectImages(true,3,labels);
			if(imps==null){
				return;
			}
			// if the spectra are available, get them
			float[][][] spectra=null;
			Object[] data=null;
			// IJ.log("test1");
			if(imps[0]!=null&&imps[1]!=null&&imps[2]!=null){
				ImageWindow iw=imps[2].getWindow();
				if(iw.getClass().getName().equals("jguis.PlotWindow4")){
					float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
					data=jutils.stack2array(imps[1].getStack());
					Object[] coef=jutils.stack2array(imps[0].getStack());
					spectra=new float[fb.nobjects][2][];
					for(int i=0;i<fb.nobjects;i++){
						spectra[i][0]=fb.get_object_spectrum(objects,(i+1),data,"Sum");
						spectra[i][1]=new float[yvals[0].length];
						float[] tempcoef=fb.get_object_spectrum(objects,(i+1),coef,"Sum");
						for(int j=0;j<yvals[0].length;j++){
							for(int k=0;k<5;k++){
								spectra[i][1][j]+=tempcoef[k]*yvals[k][j];
							}
						}
					}
				}
			}
			// IJ.log("test2");
			SkyPanel sp=new SkyPanel();
			sp.init(imp,imps[0],objects,areas,arearank,fb,showmask,spectra,data);
			SkyPanel.launch_frame(sp);
			this.getParent().setVisible(false);
			return;
		}
		if(e.getSource()==save_objects){
			boolean oldshowmask=showmask;
			showmask=true;
			boolean oldshowrank=showrank;
			showrank=false;
			update_image();
			new FileSaver(imp).saveAsTiff();
			showmask=oldshowmask;
			showrank=oldshowrank;
			update_image();
		}
		if(e.getSource()==log_number){
			IJ.log(""+fb.nobjects);
		}
		if(e.getSource()==dilate_objects){
			fb.dilateobjects(objects,false);
		}
		if(e.getSource()==erode_objects){
			fb.erodeobjects(objects);
		}
		if(e.getSource()==fill_holes){
			fb.fill_holes(objects);
		}
		if(e.getSource()==obj_stats){
			get_stats();
		}
		update_image();
	}

	public void get_stats(){
		Object[] windowlist=jutils.getImageWindowList(true);
		String[] titles=(String[])windowlist[1];
		int[] ids=(int[])windowlist[0];
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Measurement Image",titles,titles[0]);
		gd.addChoice("Measurement Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd.addCheckbox("Outside Edge Measure?",false);
		gd.addNumericField("Edge Thickness (if edge measure)",4,0);
		gd.addCheckbox("Show Edge Image (if edge measure",false);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		int index=gd.getNextChoiceIndex();
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		boolean circ=gd.getNextBoolean();
		int circrad=(int)gd.getNextNumber();
		boolean showedge=gd.getNextBoolean();
		Object[] measurement=null;
		if(index==ids.length){
			measurement=new Object[1];
			measurement[0]=imp.getStack().getPixels(1);
		}else{
			ImagePlus imp=WindowManager.getImage(ids[index]);
			ImageStack mstack=imp.getStack();
			measurement=new Object[mstack.getSize()];
			for(int i=0;i<measurement.length;i++)
				measurement[i]=mstack.getPixels(i+1);
		}
		StringBuffer labels=new StringBuffer();
		labels.append("objid");
		for(int i=0;i<measurement.length;i++)
			labels.append("\tslice"+(i+1));
		TextWindow tw=new TextWindow("Object Stats",labels.toString(),"",400,400);
		float[] tempobj=null;
		if(circ){
			tempobj=objects.clone();
			for(int i=0;i<circrad;i++)
				fb.dilateobjects(objects,false);
			for(int i=0;i<objects.length;i++){
				if(tempobj[i]>0.0f)
					objects[i]=0.0f;
			}
			if(showedge) new ImagePlus("Edge Objects",new FloatProcessor(fb.width,fb.height,objects.clone(),null)).show();
		}
		// new ImagePlus("Circ Objects",new
		// FloatProcessor(fb.width,fb.height,objects,null)).show();
		for(int i=0;i<fb.nobjects;i++){
			float[] spectrum=fb.get_object_spectrum(objects,i+1,measurement,stat);
			StringBuffer sb2=new StringBuffer();
			sb2.append(""+(i+1));
			for(int j=0;j<spectrum.length;j++)
				sb2.append("\t"+spectrum[j]);
			tw.append(sb2.toString()+"\n");
		}
		if(circ){
			objects=tempobj;
		}
	}

	private void separate_objects(Roi roi){
		Polygon poly;
		if(roi instanceof Line){
			int[] xpts={((Line)roi).x1,((Line)roi).x2};
			int[] ypts={((Line)roi).y1,((Line)roi).y2};
			// IJ.log(""+xpts[0]+" , "+ypts[0]+" , "+xpts[1]+" , "+ypts[1]);
			poly=new Polygon(xpts,ypts,2);
		}else{
			poly=roi.getPolygon();
		}
		fb.separateobjects(objects,poly,false);
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==mask_check){
			showmask=mask_check.getState();
			update_image();
		}
		if(e.getSource()==rank_check){
			showrank=rank_check.getState();
			update_image();
		}
	}

	public void mouseMoved(MouseEvent e){
		int x=e.getX();
		int y=e.getY();
		ImageCanvas ic=imp.getCanvas();
		int ox=(int)ic.offScreenX(x);
		int oy=(int)ic.offScreenY(y);
		if(ox<imp.getWidth()&&oy<imp.getHeight()){
			float tempid=objects[ox+imp.getWidth()*oy];
			if(tempid>0.0f){
				int sortid=arearank[(int)tempid-1]+1;
				idlabel.setText("id = "+(int)tempid+" rank = "+sortid);
				arealabel.setText("area = "+areas[(int)tempid-1]);
			}
		}
	}

	public void mouseDragged(MouseEvent e){

	}

	public void update_image(){
		areas=fb.get_areas(objects);
		int[] temprank=jsort.get_javasort_order(areas);
		arearank=jsort.get_javasort_order(temprank);
		for(int i=0;i<fb.nobjects;i++){
			arearank[i]=fb.nobjects-arearank[i]-1;
		}
		jutils.update_threshold_image(imp,showmask,fb,objects,areas,showrank);
		nobjects.setText("# of objects = "+fb.nobjects);
	}

}
