/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import jalgs.*;
import jalgs.jseg.findblobs3;
import jalgs.jseg.measure_object;

import java.awt.*;
import java.awt.event.*;

import ij.*;

public class interactive_counter_panel extends Panel implements ActionListener,MouseMotionListener,MouseListener{

	private Button inc_obj_button,dec_obj_button,del_obj_button,add_obj_button,out_stats_button;
	private Label countlabel,objlabel;
	public ImagePlus imp;
	public int nobjects,currobj,currx,curry;
	public findblobs3 fb;
	public int[] objcount;
	public float[] objects;
	public Polygon[] outlines;
	public float[][] coords,objstats;
	public float medarea,percentile;

	public static Frame launch_frame(interactive_counter_panel panel){
		final Frame f=new Frame("Interactive Object Counter");
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
		panel.setBounds(10,40,180,250);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(200,300));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(float[] objects,findblobs3 fb,float[] image,float percentile){
		setLayout(null);
		this.objects=objects;
		this.fb=fb;
		this.percentile=percentile;
		init_objects();
		currobj=-1;
		create_image(image);
		inc_obj_button=new Button("+ obj");
		inc_obj_button.setBounds(10,10,50,30);
		inc_obj_button.addActionListener(this);
		add(inc_obj_button);
		dec_obj_button=new Button("- obj");
		dec_obj_button.setBounds(10,50,50,30);
		dec_obj_button.addActionListener(this);
		add(dec_obj_button);
		del_obj_button=new Button("del obj");
		del_obj_button.setBounds(10,90,50,30);
		del_obj_button.addActionListener(this);
		add(del_obj_button);
		add_obj_button=new Button("add obj");
		add_obj_button.setBounds(10,130,50,30);
		add_obj_button.addActionListener(this);
		add(add_obj_button);
		out_stats_button=new Button("out stats");
		out_stats_button.setBounds(70,10,60,30);
		out_stats_button.addActionListener(this);
		add(out_stats_button);
		countlabel=new Label("# objects = "+nobjects);
		countlabel.setBounds(10,160,100,20);
		add(countlabel);
		objlabel=new Label("Object ct = 0");
		objlabel.setBounds(10,190,100,20);
		add(objlabel);
		imp.getCanvas().addMouseMotionListener(this);
		imp.getCanvas().addMouseListener(this);
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==out_stats_button){
			table_tools.create_table("Object Stats",objstats,new String[]{"Area","Circ"});
			return;
		}
		if(currobj>=0){
			if(e.getSource()==inc_obj_button){
				objcount[currobj]++;
				nobjects++;
			}
			if(e.getSource()==dec_obj_button){
				if(objcount[currobj]>1){
					objcount[currobj]--;
					nobjects--;
				}
			}
			if(e.getSource()==del_obj_button){
				delete_object();
			}
		}else{
			if(e.getSource()==add_obj_button){
				add_object();
			}
		}
		update_position();
		update_image();
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
	}

	public void mouseClicked(MouseEvent arg0){
	}

	public void mouseEntered(MouseEvent arg0){
	}

	public void mouseExited(MouseEvent arg0){
	}

	public void mousePressed(MouseEvent arg0){
	}

	public void mouseReleased(MouseEvent arg0){
		update_position();
	}

	public void init_objects(){
		int[] temp=fb.get_areas(objects);
		objstats=new float[temp.length][2];
		float[] areas=new float[temp.length];
		for(int i=0;i<temp.length;i++){
			areas[i]=(float)temp[i];
			objstats[i][0]=areas[i];
		}
		// medarea=jstatistics.getstatistic("Median",areas,null);
		float[] extras={percentile};
		medarea=jstatistics.getstatistic("Percentile",areas,extras);
		// objects greater than twice median area are counted more than once
		nobjects=fb.nobjects;
		objcount=new int[fb.nobjects];
		for(int i=0;i<fb.nobjects;i++){
			objcount[i]=1;
			if(areas[i]>2.0f*medarea){
				int extraobjects=(int)(areas[i]/medarea)-1;
				nobjects+=extraobjects;
				objcount[i]+=extraobjects;
			}
		}
	}

	private float[][] delete_table_row(float[][] table,int row){
		float[][] newtable=new float[table.length-1][];
		for(int i=0;i<row;i++)
			newtable[i]=table[i];
		for(int i=row+1;i<table.length;i++)
			newtable[i-1]=table[i];
		return newtable;
	}

	private float[][] add_table_row(float[][] table,float[] newrow){
		float[][] newtable=new float[table.length+1][];
		for(int i=0;i<table.length;i++)
			newtable[i]=table[i];
		newtable[table.length]=newrow;
		return newtable;
	}

	public void delete_object(){
		// here we delete and contract the appropriate arrays
		int oldobjcount=objcount[currobj];
		coords=delete_table_row(coords,currobj);
		objstats=delete_table_row(objstats,currobj);
		int[] newobjcount=new int[objcount.length-1];
		for(int i=0;i<currobj;i++){
			newobjcount[i]=objcount[i];
		}
		for(int i=currobj+1;i<objcount.length;i++){
			newobjcount[i-1]=objcount[i];
		}
		objcount=newobjcount;
		Polygon[] newoutlines=new Polygon[outlines.length-1];
		for(int i=0;i<currobj;i++){
			newoutlines[i]=outlines[i];
		}
		for(int i=currobj+1;i<outlines.length;i++){
			newoutlines[i-1]=outlines[i];
		}
		outlines=newoutlines;
		nobjects-=oldobjcount;
		currobj=-1;
	}

	public void add_object(){
		// here we add to the appropriate arrays
		coords=add_table_row(coords,new float[]{currx,curry});
		objstats=add_table_row(objstats,new float[]{medarea,1.0f});
		int[] newobjcount=new int[objcount.length+1];
		System.arraycopy(objcount,0,newobjcount,0,objcount.length);
		newobjcount[objcount.length]=1;
		objcount=newobjcount;
		Polygon[] newoutlines=new Polygon[outlines.length+1];
		for(int i=0;i<outlines.length;i++){
			newoutlines[i]=outlines[i];
		}
		newoutlines[outlines.length]=make_circle(currx,curry);
		outlines=newoutlines;
		nobjects++;
	}

	public void drawobject(Polygon poly,int id){
		Rectangle r=poly.getBounds();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					objects[j+i*fb.width]=(float)(id+1);
				}
			}
		}
	}

	public Polygon make_circle(int x,int y){
		float radius=(float)Math.sqrt(medarea/Math.PI);
		return jutils.circle2poly(x,y,radius);
	}

	public void update_position(){
		if(imp.getRoi()==null)
			return;
		Rectangle r=imp.getRoi().getBounds();
		currx=r.x;
		curry=r.y;
		get_curr_obj();
		int currcount=0;
		if(currobj>=0)
			currcount=objcount[currobj];
		objlabel.setText("Object ct = "+currcount);
		countlabel.setText("# objects = "+nobjects);
	}

	public void get_curr_obj(){
		// do hit testing to find out which object we are in if any
		currobj=-1;
		for(int i=0;i<outlines.length;i++){
			if(outlines[i].contains(currx,curry)){
				currobj=i;
				return;
			}
		}
	}

	public void create_image(float[] image){
		outlines=fb.get_object_outlines(objects);
		for(int i=0;i<outlines.length;i++){
			float perim=(float)fb.get_perimeter(outlines[i]);
			objstats[i][1]=4.0f*(float)Math.PI*(objstats[i][0]/(perim*perim));
		}
		coords=measure_object.centroids(objects,fb.width,fb.height);
		imp=jutils.create_threshold_image(fb,objects,image,false,objcount,true,outlines,coords);
	}

	public void update_image(){
		jutils.update_threshold_image(imp,false,fb,objects,objcount,true,outlines,coords);
	}

}
