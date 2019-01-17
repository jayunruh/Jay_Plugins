/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.gui.GenericDialog;
import ij.plugin.frame.Editor;
import jalgs.delimit_string;

import java.awt.Button;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class Layout_Generator extends Panel implements MouseListener,MouseMotionListener,ActionListener{
	private Rectangle[] rects;
	private Rectangle oldrect;
	private int[] recttypes;
	private String[] rectnames;
	private Object[] recttext;
	private boolean[] haslistener;
	private Button gen_code_button;
	private int inrect,inrectID,oldxpos,oldypos,nrects;
	private FontMetrics fm;
	private static final int maxrects=1000;
	// public static int width,height;
	public static int width;
	public static int height;
	public static String[] types={"TextField","Label","Button","Checkbox","Choice","Scrollbar"};

	public static void launch_frame(Layout_Generator panel){
		final Frame f=new Frame("Layout Generator");
		f.setLocation(10,10);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				/*
				 * Component[] comps=f.getComponents(); for(int
				 * i=0;i<comps.length;i++){ comps[i].setVisible(false); }
				 */
				f.dispose();
			}
		});

		f.setLayout(null);
		Insets ins=f.getInsets();
		int tempheight=Layout_Generator.height+ins.bottom+ins.top+100;
		int tempwidth=Layout_Generator.width+ins.left+ins.right;
		panel.setBounds(ins.left+5,ins.top+5,Layout_Generator.width,Layout_Generator.height+100);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(tempwidth,tempheight));
		f.setVisible(true);
		panel.requestFocus();
	}

	public void paint(Graphics g){
		if(fm==null){
			fm=g.getFontMetrics();
		}
		g.clearRect(0,0,width-1,height-1);
		g.drawRect(0,0,width-1,height-1);
		for(int i=0;i<nrects;i++){
			g.setClip(0,0,width,height);
			g.drawRect(rects[i].x,rects[i].y,rects[i].width,rects[i].height);
			// draw the text for each item
			g.setClip(rects[i].x,rects[i].y,rects[i].width,rects[i].height);
			if(recttypes[i]==0){
				draw_left_centered_string(g,(String)recttext[i],rects[i].x+4,rects[i].y+rects[i].height/2);
			}
			if(recttypes[i]==1){
				draw_left_centered_string(g,(String)recttext[i],rects[i].x+1,rects[i].y+rects[i].height/2);
			}
			if(recttypes[i]==2){
				draw_centered_string(g,(String)recttext[i],rects[i].x+rects[i].width/2,rects[i].y+rects[i].height/2);
			}
			if(recttypes[i]==3){
				draw_left_centered_string(g,(String)recttext[i],rects[i].x+1,rects[i].y+rects[i].height/2);
			}
			if(recttypes[i]==4){
				draw_left_centered_string(g,(String)recttext[i],rects[i].x+4,rects[i].y+rects[i].height/2);
			}
		}
		g.setClip(null);
	}

	public void update(Graphics g){
		paint(g);
	}

	public void draw_centered_string(Graphics g,String out,int xcenter,int ycenter){
		int width=fm.stringWidth(out);
		int height=fm.getMaxAscent();// +fm.getMaxDescent();
		g.drawString(out,xcenter-width/2,ycenter+height/2);
	}

	public void draw_left_centered_string(Graphics g,String out,int xleft,int ycenter){
		int height=fm.getMaxAscent();// +fm.getMaxDescent();
		g.drawString(out,xleft,ycenter+height/2);
	}

	public void init(int width,int height){
		// here we start an empty layout
		Layout_Generator.width=width;
		Layout_Generator.height=height;
		rects=new Rectangle[maxrects];
		recttypes=new int[maxrects];
		rectnames=new String[maxrects];
		recttext=new Object[maxrects];
		haslistener=new boolean[maxrects];
		setLayout(null);
		gen_code_button=new Button("Gen Code");
		gen_code_button.setBounds(10,height+10,60,30);
		gen_code_button.addActionListener(this);
		add(gen_code_button);
		addMouseListener(this);
		addMouseMotionListener(this);
		inrect=0;
		nrects=0;
		fm=null;
		repaint();
	}

	/*
	 * public void init(int width,int height,String code){
	 * 
	 * }
	 */

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==gen_code_button){
			generate_code();
		}
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
		int xpos=e.getX();
		int ypos=e.getY();
		int[] hit=hit_test(xpos,ypos);
		if(e.getButton()==MouseEvent.BUTTON3){
			// right mouse button clicked
			if(hit!=null){
				// here we just modify an existing object
				inrectID=hit[0];
				edit_item(inrectID);
				inrect=0;
			}
		}else{
			if(xpos<width&&ypos<height){
				if(hit==null){
					// start a new item
					inrect=3;
					inrectID=nrects;
					rects[nrects]=new Rectangle(xpos,ypos,1,1);
					nrects++;
				}else{
					// modify an existing item
					if(hit[1]==0){
						inrect=1;
					}else{
						inrect=2;
					}
					inrectID=hit[0];
				}
				oldxpos=xpos;
				oldypos=ypos;
				oldrect=new Rectangle(rects[inrectID]);
			}
		}
	}

	public void mouseReleased(MouseEvent e){
		if(inrect==3){
			recttypes[nrects-1]=0;
			rectnames[nrects-1]="";
			recttext[nrects-1]="";
			haslistener[nrects-1]=false;
			if(!edit_item(nrects-1)){
				inrect=0;
				nrects--;
				repaint();
				return;
			}
		}
		inrect=0;
		repaint();
	}

	public boolean edit_item(int itemID){
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("delete item",false);
		gd.addChoice("Component Type",types,types[recttypes[itemID]]);
		gd.addStringField("Component Name",rectnames[itemID]);
		String temp;
		boolean wasscrollbar=(recttypes[itemID]==5);
		if(recttext[itemID] instanceof String){
			temp=recttext[itemID].toString();
		}else{
			if(recttext[itemID] instanceof String[]){
				String[] temp2=(String[])recttext[itemID];
				StringBuffer temp3=new StringBuffer();
				temp3.append(temp2[0]);
				for(int i=1;i<temp2.length;i++){
					temp3.append(","+temp2[i]);
				}
				temp=temp3.toString();
			}else{
				temp="";

			}
		}
		gd.addStringField("Component Text (or csv list)",temp);
		gd.addNumericField("X",rects[itemID].x,0);
		gd.addNumericField("Y",rects[itemID].y,0);
		gd.addNumericField("Width",rects[itemID].width,0);
		gd.addNumericField("Height",rects[itemID].height,0);
		gd.addCheckbox("Has Listener?",haslistener[itemID]);
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
		if(gd.getNextBoolean()){
			// here we need to delete the item
			for(int i=itemID;i<nrects;i++){
				rects[i]=rects[i+1];
				rectnames[i]=rectnames[i+1];
				recttext[i]=recttext[i+1];
				recttypes[i]=recttypes[i+1];
				haslistener[i]=haslistener[i+1];
				nrects--;
				repaint();
				return true;
			}
		}
		recttypes[itemID]=gd.getNextChoiceIndex();
		rectnames[itemID]=gd.getNextString();
		String temp4=gd.getNextString();
		rects[itemID].x=(int)gd.getNextNumber();
		rects[itemID].y=(int)gd.getNextNumber();
		rects[itemID].width=(int)gd.getNextNumber();
		rects[itemID].height=(int)gd.getNextNumber();
		if(recttypes[itemID]<4){
			recttext[itemID]=temp4;
		}else{
			if(recttypes[itemID]==4){
				delimit_string ds=new delimit_string(',');
				int cols=ds.getnumcolumns(temp4);
				String[] options=ds.delim2string(temp4,cols);
				recttext[itemID]=options;
			}else{
				int[] scrolloptions={0,0,10,1};
				if(wasscrollbar){
					scrolloptions=(int[])recttext[itemID];
				}
				GenericDialog gd2=new GenericDialog("Scrollbar Options");
				gd2.addCheckbox("Vertical?",scrolloptions[0]==1);
				gd2.addNumericField("Minimum Value",scrolloptions[1],0);
				gd2.addNumericField("Maximum Value",scrolloptions[2],0);
				gd2.addNumericField("Bar Thickness",scrolloptions[3],0);
				gd2.showDialog();
				if(gd2.wasCanceled()){
					inrect=0;
					nrects--;
					repaint();
					return false;
				}
				if(gd2.getNextBoolean()){
					scrolloptions[0]=1;
				}
				scrolloptions[1]=(int)gd2.getNextNumber();
				scrolloptions[2]=(int)gd2.getNextNumber();
				scrolloptions[3]=(int)gd2.getNextNumber();
				recttext[itemID]=scrolloptions;
			}
		}
		haslistener[itemID]=gd.getNextBoolean();
		return true;
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		int xpos=e.getX();
		int ypos=e.getY();
		if(inrect==1){
			// here we are simply moving the rectangle
			rects[inrectID].x=oldrect.x+(xpos-oldxpos);
			rects[inrectID].y=oldrect.y+(ypos-oldypos);
			repaint();
		}else{
			if(inrect==2||inrect==3){
				// here we are dragging an edge of the rectangle
				int oldxcorn=oldrect.x;
				int oldycorn=oldrect.y;
				int newwidth=Math.abs(xpos-oldxcorn);
				int newheight=Math.abs(ypos-oldycorn);
				if(xpos<oldxcorn){
					oldxcorn=xpos;
				}
				if(ypos<oldycorn){
					oldycorn=ypos;
				}
				rects[inrectID]=new Rectangle(oldxcorn,oldycorn,newwidth,newheight);
				repaint();
			}
		}
	}

	private int[] hit_test(int x,int y){
		int th=2;
		int[] retvals=new int[2];
		for(int i=0;i<nrects;i++){
			if(x>=(rects[i].x-th)&&x<=(rects[i].x+rects[i].width+th)&&y>=(rects[i].y-th)&&y<=(rects[i].y+rects[i].height+th)){
				// we are inside this rectangle
				retvals[0]=i;
				if(x>(rects[i].x+th)&&x<(rects[i].x+rects[i].width-th)&&y>(rects[i].y+th)&&y<(rects[i].y+rects[i].height-th)){
					// we are inside the inside of this rectangle
					retvals[1]=0;
				}else{
					// we are inside the border of this rectangle
					retvals[1]=1;
				}
				return retvals;
			}
		}
		// we are not inside any rectangles
		return null;
	}

	private void generate_code(){
		StringBuffer decl=new StringBuffer();
		StringBuffer inits=new StringBuffer();
		decl.append("//starting gui declarations\n");
		inits.append("//starting gui definitions\n");
		inits.append("setLayout(null);\n");
		for(int i=0;i<nrects;i++){
			decl.append("private "+types[recttypes[i]]+" "+rectnames[i]+";\n");
			if(recttypes[i]<3){
				inits.append(rectnames[i]+"=new "+types[recttypes[i]]+"(\""+(String)recttext[i]+"\");\n");
				if(haslistener[i]&&recttypes[i]!=1){
					inits.append(rectnames[i]+".addActionListener(this);\n");
				}
			}else{
				if(recttypes[i]==3){
					inits.append(rectnames[i]+"=new "+types[recttypes[i]]+"(\""+(String)recttext[i]+"\","+"true);\n");
					if(haslistener[i]){
						inits.append(rectnames[i]+".addItemListener(this);\n");
					}
				}else{
					if(recttypes[i]==4){
						inits.append(rectnames[i]+"=new "+types[recttypes[i]]+"();\n");
						String[] items=(String[])recttext[i];
						for(int j=0;j<items.length;j++){
							inits.append(rectnames[i]+".add(\""+items[j]+"\");\n");
						}
						if(haslistener[i]){
							inits.append(rectnames[i]+".addItemListener(this);\n");
						}
					}else{
						int[] so=(int[])recttext[i];
						if(so[0]==0){
							inits.append(rectnames[i]+"=new "+types[recttypes[i]]+"(Scrollbar.HORIZONTAL"+","+so[1]+","+so[3]+","+so[1]+","+(so[2]-1+so[3])+");\n");
						}else{
							inits.append(rectnames[i]+"=new "+types[recttypes[i]]+"(Scrollbar.VERTICAL"+","+so[1]+","+so[3]+","+so[1]+","+(so[2]-1+so[3])+");\n");
						}
						if(haslistener[i]){
							inits.append(rectnames[i]+".addAdjustmentListener(this);\n");
						}
					}
				}
			}
			inits.append(rectnames[i]+".setBounds("+rects[i].x+","+rects[i].y+","+rects[i].width+","+rects[i].height+");\n");
			inits.append("add("+rectnames[i]+");\n");
		}
		inits.append("//ending gui definitions\n");
		decl.append("//ending gui declarations\n");
		Editor ed=new Editor();
		ed.create("generated_code",decl.toString()+inits.toString());
	}

}
