/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.CompositeImage;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import jalgs.algutils;
import jalgs.jdataio;
import jalgs.jstatistics;

import java.awt.Button;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Label;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class plate_viewer_panel extends Panel implements ActionListener,ItemListener,MouseListener,MouseMotionListener,ImageListener{

	private Button quit_button;
	private Choice rowsel,colsel,fisel,frsel,zsel;
	private Label rowlab,collab,filab,frlab,zlab;
	public ImagePlus imp;
	public Image plateimg;
	public String dir;
	public int[][] lims;
	public String[][][][][][] nmatrix;
	public int width,height,dtype;

	public static Frame launch_frame(plate_viewer_panel panel){
		final Frame f=new Frame("Plate Viewer");
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
		panel.setBounds(10,40,400,500);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(420,520));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}
	
	/********************************
	 * this assembles a file list from a set of PE files with format rrrccc-f-tttzzzccc.tif
	 * format could also be r##c##f##p##-ch#sk#fk#fl#.tiff--not supported for now (row,col,field,plane,ch,t?,?,?)
	 * @param dir
	 */
	public static List<List<String>> makePEFlist(String dir){
		jdataio jdio=new jdataio();
		String[] flist1=jdio.get_numeric_sorted_string_list(dir,null);
		List<List<String>> flist=new ArrayList<List<String>>();
		for(int i=0;i<flist1.length;i++){
			if(!(flist1[i].endsWith(".tif") || flist1[i].endsWith(".tiff"))) continue;
			String[] trow=null;
			if(flist1[i].substring(0,1).equals("r")){
				String[] split1=flist1[i].split("-");
				String[] first4=new String[]{split1[0].substring(1,3),split1[0].substring(4,6),split1[0].substring(7,9),split1[0].substring(10,12)};
				int skpos=split1[1].indexOf("sk");
				int fkpos=split1[1].indexOf("fk");
				int flpos=split1[1].indexOf("fl");
				String[] last4=new String[]{split1[1].substring(2,skpos),split1[1].substring(skpos+2,fkpos),split1[1].substring(fkpos+2,flpos),split1[1].substring(flpos+2,split1[1].length()-4)};
				trow=new String[]{flist1[i],first4[0],first4[1],first4[2],last4[1],first4[3],last4[0]};
			} else {
				trow=new String[]{flist1[i],flist1[i].substring(0,3),flist1[i].substring(3,6),flist1[i].substring(7,8),flist1[i].substring(9,12),flist1[i].substring(12,15),flist1[i].substring(15,18)};
			}
			flist.add(Arrays.asList(trow));
		}
		return flist;
	}
	
	public static String[][][][][][] getNameMatrix(List<List<String>> flist){
		int[] rlist=algutils.convert_arr_int(table_tools.get_column_array(flist,1));
		int[] clist=algutils.convert_arr_int(table_tools.get_column_array(flist,2));
		int[] filist=algutils.convert_arr_int(table_tools.get_column_array(flist,3));
		int[] tlist=algutils.convert_arr_int(table_tools.get_column_array(flist,4));
		int[] zlist=algutils.convert_arr_int(table_tools.get_column_array(flist,5));
		int[] chlist=algutils.convert_arr_int(table_tools.get_column_array(flist,6));
		int[][] lims=getLims(flist);
		int[] dims=new int[lims.length];
		for(int i=0;i<dims.length;i++){
			dims[i]=lims[i][1]-lims[i][0]+1;
		}
		String[][][][][][] nmatrix=new String[dims[0]][dims[1]][dims[2]][dims[3]][dims[4]][dims[5]];
		for(int i=0;i<rlist.length;i++){
			nmatrix[rlist[i]-lims[0][0]][clist[i]-lims[1][0]][filist[i]-lims[2][0]][tlist[i]-lims[3][0]][zlist[i]-lims[4][0]][chlist[i]-lims[5][0]]=flist.get(i).get(0);
		}
		return nmatrix;
	}
	
	public static int[][] getLims(List<List<String>> flist){
		int[] rlist=algutils.convert_arr_int(table_tools.get_column_array(flist,1));
		int[] clist=algutils.convert_arr_int(table_tools.get_column_array(flist,2));
		int[] filist=algutils.convert_arr_int(table_tools.get_column_array(flist,3));
		int[] tlist=algutils.convert_arr_int(table_tools.get_column_array(flist,4));
		int[] zlist=algutils.convert_arr_int(table_tools.get_column_array(flist,5));
		int[] chlist=algutils.convert_arr_int(table_tools.get_column_array(flist,6));
		int[][] lims=new int[6][2];
		lims[0]=new int[]{(int)jstatistics.getstatistic("Min",rlist,null),(int)jstatistics.getstatistic("Max",rlist,null)};
		lims[1]=new int[]{(int)jstatistics.getstatistic("Min",clist,null),(int)jstatistics.getstatistic("Max",clist,null)};
		lims[2]=new int[]{(int)jstatistics.getstatistic("Min",filist,null),(int)jstatistics.getstatistic("Max",filist,null)};
		lims[3]=new int[]{(int)jstatistics.getstatistic("Min",tlist,null),(int)jstatistics.getstatistic("Max",tlist,null)};
		lims[4]=new int[]{(int)jstatistics.getstatistic("Min",zlist,null),(int)jstatistics.getstatistic("Max",zlist,null)};
		lims[5]=new int[]{(int)jstatistics.getstatistic("Min",chlist,null),(int)jstatistics.getstatistic("Max",chlist,null)};
		return lims;
	}
	
	public static String[] makeIntList(int start,int end){
		String[] outlist=new String[end-start+1];
		for(int i=start;i<=end;i++){
			outlist[i-start]=""+i;
		}
		return outlist;
	}

	/***********************
	 * the main initializer, flist should have columns: name,row,col,field,time,z,ch
	 * @param dir
	 * @param flist
	 */
	public void init(String dir,List<List<String>> flist){
		setLayout(null);
		this.dir=dir;
		this.lims=getLims(flist);
		this.nmatrix=getNameMatrix(flist);
		this.width=-1;
		//this.imp=IJ.openImage(dir+flist.get(0).get(0));
		int[][] tlims=this.lims;
		Object[] tstack=getStack(tlims[0][0],tlims[1][0],tlims[2][0],tlims[3][0],tlims[4][0]);
		ImagePlus timp=new ImagePlus("Plate_Stack",jutils.array2stack(tstack,width,height));
		this.imp=new CompositeImage(timp,CompositeImage.COMPOSITE);
		this.imp.show();
		ImagePlus.addImageListener(this);
		quit_button=new Button("Quit");
		quit_button.setBounds(10,10,100,30);
		quit_button.addActionListener(this);
		rowlab=new Label("Row:");
		rowlab.setBounds(10,10+40,50,15);
		rowsel=new Choice();
		rowsel.setBounds(10+60,10+40,100,15);
		rowsel.addItemListener(this);
		String[] rowlist=makeIntList(lims[0][0],lims[0][1]);
		for(int i=0;i<rowlist.length;i++) rowsel.add(rowlist[i]);
		collab=new Label("Col:");
		collab.setBounds(10,10+40+30,50,15);
		colsel=new Choice();
		colsel.setBounds(10+60,10+40+30,100,15);
		colsel.addItemListener(this);
		String[] collist=makeIntList(lims[1][0],lims[1][1]);
		for(int i=0;i<collist.length;i++) colsel.add(collist[i]);
		filab=new Label("Field:");
		filab.setBounds(10,10+40+60,50,15);
		fisel=new Choice();
		fisel.setBounds(10+60,10+40+60,100,15);
		fisel.addItemListener(this);
		String[] filist=makeIntList(lims[2][0],lims[2][1]);
		for(int i=0;i<filist.length;i++) fisel.add(filist[i]);
		frlab=new Label("Frame:");
		frlab.setBounds(10,10+40+90,50,15);
		frsel=new Choice();
		frsel.setBounds(10+60,10+40+90,100,15);
		frsel.addItemListener(this);
		String[] frlist=makeIntList(lims[3][0],lims[3][1]);
		for(int i=0;i<frlist.length;i++) frsel.add(frlist[i]);
		zlab=new Label("Slice:");
		zlab.setBounds(10,10+40+120,50,15);
		zsel=new Choice();
		zsel.setBounds(10+60,10+40+120,100,15);
		zsel.addItemListener(this);
		String[] zlist=makeIntList(lims[4][0],lims[4][1]);
		for(int i=0;i<zlist.length;i++) zsel.add(zlist[i]);
		
		add(quit_button);
		add(rowlab); add(rowsel);
		add(collab); add(colsel);
		add(filab); add(fisel);
		add(frlab); add(frsel);
		add(zlab); add(zsel);
		addMouseListener(this);
		updateImage(tlims[0][0],tlims[1][0],tlims[2][0],tlims[3][0],tlims[4][0]);
	}
	
	public void paint(Graphics g){
		g.drawImage(plateimg,20,220,this);
		Font f=g.getFont();
		g.setFont(new Font("Arial",Font.PLAIN,8));
		boolean is384=lims[0][1]>8 || lims[1][1]>12;
		String[] letts=new String[]{"A","B","C","D","E","F","G","H","I","L","M","N","O","P"};
		String[] nums=new String[24];
		for(int i=0;i<24;i++) nums[i]=""+(i+1);
		if(is384){
			for(int i=0;i<16;i++){
				g.drawString(letts[i],5,230+i*12);
			}
			for(int i=0;i<24;i++){
				g.drawString(nums[i],25+i*12,210);
			}
		} else {
			for(int i=0;i<8;i++){
				g.drawString(letts[i],5,235+i*25);
			}
			for(int i=0;i<12;i++){
				g.drawString(nums[i],30+i*25,210);
			}
		}
		g.setFont(f);
	}

	public void update(Graphics g){
		paint(g);
	}
	
	public Image drawPlate(int row,int col){
		//draw a grid of squares representing wells
		//shade the occupied ones
		//select one of the wells
		int imgwidth=300;
		int imgheight=200;
		int bgcolor=jutils.argb2intval(255,255,255,255);
		int wellcolor=jutils.argb2intval(255,200,200,200);
		int fullcolor=jutils.argb2intval(255,0,0,255);
		int selcolor=jutils.argb2intval(255,255,0,0);
		int[] pix=new int[imgwidth*imgheight];
		for(int i=0;i<pix.length;i++) pix[i]=bgcolor;
		boolean is384=lims[0][1]>8 || lims[1][1]>12;
		if(!is384){
			int wwidth=25;
			int wborder=5;
			int remwidth=wwidth-wborder;
			for(int i=0;i<12;i++){
				int xpos=i*wwidth;
				for(int j=0;j<8;j++){
					int ypos=j*wwidth;
					boolean issel=(((i+1)==col)&&((j+1)==row));
					boolean isfull=((i+1)>=lims[1][0]&&(i+1)<=lims[1][1]&&(j+1)>=lims[0][0]&&(j+1)<=lims[0][1]);
					for(int k=0;k<remwidth;k++){
						int xpos2=xpos+k;
						for(int l=0;l<remwidth;l++){
							int ypos2=ypos+l;
							if(isfull){
								if(issel) pix[xpos2+ypos2*imgwidth]=selcolor;
								else pix[xpos2+ypos2*imgwidth]=fullcolor;
							} else {
								pix[xpos2+ypos2*imgwidth]=wellcolor;
							}
						}
					}
				}
			}
		} else {
			int wwidth=12;
			int wborder=3;
			int remwidth=wwidth-wborder;
			for(int i=0;i<24;i++){
				int xpos=i*wwidth;
				for(int j=0;j<16;j++){
					int ypos=j*wwidth;
					boolean issel=(i==col)&&(j==row);
					boolean isfull=(i<lims[1][0]||i>lims[1][1])||(j<lims[0][0]||j>lims[0][1]);
					for(int k=0;k<remwidth;k++){
						int xpos2=xpos+k;
						for(int l=0;l<remwidth;l++){
							int ypos2=ypos+l;
							if(isfull){
								if(issel) pix[xpos2+ypos2*imgwidth]=selcolor;
								else pix[xpos2+ypos2*imgwidth]=fullcolor;
							} else {
								pix[xpos2+ypos2*imgwidth]=wellcolor;
							}
						}
					}
				}
			}
		}
		return (new ColorProcessor(imgwidth,imgheight,pix)).createImage();
	}
	
	public void updateImage(int row,int col,int field,int frame,int slice){
		Object[] stack=getStack(row,col,field,frame,slice);
		this.imp.setStack(jutils.array2stack(stack,this.width,this.height));
		this.imp.updateAndDraw();
		//now update the plate diagram
		plateimg=drawPlate(row,col);
		repaint();
	}
	
	public Object[] getStack(int row,int col,int field,int frame,int slice){
		int nch=lims[5][1]-lims[5][0]+1;
		Object[] stack=new Object[nch];
		int trow=row-lims[0][0];
		int tcol=col-lims[1][0];
		int tfield=field-lims[2][0];
		int tframe=frame-lims[3][0];
		int tslice=slice-lims[4][0];
		for(int i=0;i<nch;i++){
			String tname=this.nmatrix[trow][tcol][tfield][tframe][tslice][i];
			if(tname!=null){
				ImagePlus timp=IJ.openImage(this.dir+tname);
				stack[i]=timp.getProcessor().getPixels();
				if(this.width<0){
					this.width=timp.getWidth();
					this.height=timp.getHeight();
					this.dtype=algutils.get_array_type(stack[i]);
				}
			} else {
				stack[i]=algutils.create_array(this.width*this.height,this.dtype);
			}
		}
		return stack;
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
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
	}

	public void mouseClicked(MouseEvent arg0){
		int xpos=arg0.getX()-20;
		int ypos=arg0.getY()-220;
		if(xpos<0 || ypos<0) return;
		if(xpos>300 || ypos>200) return;
		boolean is384=lims[0][1]>8 || lims[1][1]>12;
		if(!is384){
			//IJ.log(""+xpos+" , "+ypos);
			float newx=(float)xpos/25.0f;
			float newy=(float)ypos/25.0f;
			int row=(int)newy+1;
			int col=(int)newx+1;
			float xrem=newx-(col-1);
			float yrem=newy-(row-1);
			if(xrem>0.8f || yrem>0.8f) return;
			if(row<lims[0][0]||row>lims[0][1]||col<lims[1][0]||col>lims[1][1]) return;
			rowsel.select(row-lims[0][0]);
			colsel.select(col-lims[1][0]);
		} else {
			//IJ.log(""+xpos+" , "+ypos);
			float newx=(float)xpos/12.0f;
			float newy=(float)ypos/12.0f;
			int row=(int)newy+1;
			int col=(int)newx+1;
			float xrem=newx-(col-1);
			float yrem=newy-(row-1);
			if(xrem>0.8f || yrem>0.8f) return;
			if(row<lims[0][0]||row>lims[0][1]||col<lims[1][0]||col>lims[1][1]) return;
			rowsel.select(row-lims[0][0]);
			colsel.select(col-lims[1][0]);
		}
		//update all of the positions and the image
		int row=rowsel.getSelectedIndex()+lims[0][0];
		int col=colsel.getSelectedIndex()+lims[1][0];
		int field=fisel.getSelectedIndex()+lims[2][0];
		int frame=frsel.getSelectedIndex()+lims[3][0];
		int slice=zsel.getSelectedIndex()+lims[4][0];
		updateImage(row,col,field,frame,slice);
	}

	public void mouseEntered(MouseEvent arg0){
	}

	public void mouseExited(MouseEvent arg0){
	}

	public void mousePressed(MouseEvent arg0){
	}

	public void mouseReleased(MouseEvent arg0){
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

	public void itemStateChanged(ItemEvent arg0){
		//update all of the positions and the image
		int row=rowsel.getSelectedIndex()+lims[0][0];
		int col=colsel.getSelectedIndex()+lims[1][0];
		int field=fisel.getSelectedIndex()+lims[2][0];
		int frame=frsel.getSelectedIndex()+lims[3][0];
		int slice=zsel.getSelectedIndex()+lims[4][0];
		updateImage(row,col,field,frame,slice);
	}

}
