/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.CompositeImage;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.StackWindow;
import ij.process.LUT;
import jalgs.algutils;
import jalgs.interpolation;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * This class is an extended StackWindow that displays carpets with different
 * binning.
 */
public class CarpetWindow extends StackWindow implements ActionListener{

	private Button shr,shl,pr,pl,home,end,fit,orig,scup,scdn,edit;
	private TextField scale;
	private Label scalelabel;
	public static ImageStack stack;
	public int width2,height2,channels; // width2 is the number of columns,
	// height2 is the number of rows
	// private float[][] carpet;
	public Object[] carpet;
	public float[] limits;
	public float binsize;
	public static int WIDTH=750;

	public CarpetWindow(String title1,Object carpet,int width,int height){
		this(title1,new Object[]{carpet},width,height,null,0);
	}

	public CarpetWindow(String title1,Object[] carpet,int width,int height){
		this(title1,carpet,width,height,null,0);
	}

	public CarpetWindow(String title1,Object[] carpet,int width,int height,LUT[] luts,int orientation){
		super(createImage(title1,carpet,width,luts));
		this.width2=width;
		this.height2=height;
		this.channels=carpet.length;
		// this.carpet=carpet;
		// need to transpose carpet if orientation is incorrect
		// careful--memory may be limited
		if(orientation==0)
			this.carpet=transposecarpet(carpet);
		else
			this.carpet=carpet;
		limits=new float[]{0.0f,WIDTH};
		binsize=1.0f;
		updateStack();
		imp.resetDisplayRange();
	}

	public Object[] transposecarpet(Object[] carpet){
		Object[] temp=new Object[carpet.length];
		for(int i=0;i<carpet.length;i++){
			temp[i]=algutils.transpose_image(carpet[i],width2,height2);
		}
		return temp;
	}

	static ImagePlus createImage(String title1,Object[] carpet,int height,LUT[] luts){
		stack=new ImageStack(WIDTH,height);
		for(int j=0;j<carpet.length;j++){
			stack.addSlice("",new float[WIDTH*height]);
		}
		ImagePlus imp=new ImagePlus(title1,stack);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(carpet.length,1,1);
		if(carpet.length>1){
			CompositeImage ci=new CompositeImage(imp,CompositeImage.COMPOSITE);
			if(luts!=null)
				ci.setLuts(luts);
			return ci;
		}else{
			if(luts!=null)
				imp.getProcessor().setLut(luts[0]);
			return imp;
		}
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(float[] limits){
		this.limits=limits;
		binsize=(limits[1]-limits[0])/WIDTH;
		scale.setText(""+1.0f/binsize);
		updateStack();
	}

	public void setLimits(float xmin){
		float[] newlimits={xmin,xmin+binsize*WIDTH};
		setLimits(newlimits);
	}

	public void scaleroi(){
		Rectangle rect=imp.getProcessor().getRoi();
		if(rect!=null){
			float[] newlimits={rect.x*binsize+limits[0],(rect.x+rect.width)*binsize+limits[0]};
			imp.killRoi();
			setLimits(newlimits);
		}
	}

	public void updateStack(){
		// in order to bin an arbitrary amount, we have to integrate from limit
		// 0 to limit 1 in binsize chunks
		float[][] limits2=new float[channels][2];
		for(int i=0;i<channels;i++){
			limits2[i][0]=(float)stack.getProcessor(i+1).getMin();
			limits2[i][1]=(float)stack.getProcessor(i+1).getMax();
		}
		stack=new ImageStack(WIDTH,width2);
		for(int i=0;i<channels;i++){
			float[] temp3=new float[WIDTH*width2]; // here is the binned
			// slice
			for(int j=0;j<width2;j++){
				Object temp=algutils.get_image_row(carpet[i],height2,width2,j);
				float[] temp2=integrate_line(temp,height2);
				System.arraycopy(temp2,0,temp3,j*WIDTH,WIDTH);
			}
			stack.addSlice("",temp3);
		}
		imp.setStack(null,stack);
		if(channels>1){
			LUT[] lut=((CompositeImage)imp).getLuts();
			for(int i=0;i<channels;i++){
				stack.getProcessor(i+1).setMinAndMax(limits2[i][0],limits2[i][1]);
				lut[i].min=limits2[i][0];
				lut[i].max=limits2[i][1];
			}
		}else{
			imp.setDisplayRange(limits2[0][0],limits2[0][1]);
		}
		imp.updateAndDraw();
	}

	public float[] integrate_line(Object oldline,int len){
		float[] newline=new float[WIDTH];
		float pos=limits[0];
		for(int i=0;i<WIDTH;i++){
			newline[i]=interpolation.integrate(oldline,len,pos,pos+binsize,false);
			pos+=binsize;
		}
		return newline;
	}

	/** Displays the plot. */
	public void draw(){
		Panel buttons=new Panel();
		buttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
		home=new Button("|<<");
		home.addActionListener(this);
		buttons.add(home);
		pl=new Button("<<");
		pl.addActionListener(this);
		buttons.add(pl);
		shl=new Button("<");
		shl.addActionListener(this);
		buttons.add(shl);
		shr=new Button(">");
		shr.addActionListener(this);
		buttons.add(shr);
		pr=new Button(">>");
		pr.addActionListener(this);
		buttons.add(pr);
		end=new Button(">>|");
		end.addActionListener(this);
		buttons.add(end);
		add(buttons);
		Panel buttons2=new Panel();
		buttons2.setLayout(new FlowLayout(FlowLayout.RIGHT));
		orig=new Button("reset");
		orig.addActionListener(this);
		buttons2.add(orig);
		fit=new Button("fit");
		fit.addActionListener(this);
		buttons2.add(fit);
		scup=new Button("zoom+");
		scup.addActionListener(this);
		buttons2.add(scup);
		scdn=new Button("zoom-");
		scdn.addActionListener(this);
		buttons2.add(scdn);
		scalelabel=new Label("Scaling");
		buttons2.add(scalelabel);
		scale=new TextField("1.00000");
		scale.addActionListener(this);
		buttons2.add(scale);
		edit=new Button("Edit");
		edit.addActionListener(this);
		buttons2.add(edit);
		add(buttons2);
		pack();
		// coordinates.setText("");
		updateStack();
	}

	void editStack(){
		GenericDialog gd=new GenericDialog("Stack Options");
		gd.addNumericField("x min",limits[0],5,10,null);
		gd.addCheckbox("Scale Roi",false);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		float xmin=(float)gd.getNextNumber();
		if(limits[0]!=xmin)
			setLimits(xmin);
		else{
			if(gd.getNextBoolean())
				scaleroi();
		}
	}

	public void actionPerformed(ActionEvent e){
		Object b=e.getSource();
		if(b==shr){
			limits[0]+=1.0f;
			limits[1]+=1.0f;
		}else if(b==shl){
			limits[0]-=1.0f;
			limits[1]-=1.0f;
		}else if(b==pr){
			limits[0]=limits[1];
			limits[1]=limits[0]+binsize*WIDTH;
		}else if(b==pl){
			limits[1]=limits[0];
			limits[0]=limits[1]-binsize*WIDTH;
		}else if(b==home){
			limits[0]=0.0f;
			limits[1]=limits[0]+binsize*WIDTH;
		}else if(b==end){
			limits[1]=height2-1;
			limits[0]=limits[1]-binsize*WIDTH;
		}else if(b==fit){
			limits[1]=height2-1;
			limits[0]=0.0f;
			binsize=(float)height2/(float)WIDTH;
		}else if(b==orig){
			binsize=1.0f;
			float center=0.5f*(limits[1]+limits[0]);
			limits[1]=center+0.5f*binsize*WIDTH;
			limits[0]=limits[1]-binsize*WIDTH;
		}else if(b==scup){
			float center=0.5f*(limits[1]+limits[0]);
			binsize*=0.5f;
			limits[1]=center+0.5f*binsize*WIDTH;
			limits[0]=limits[1]-binsize*WIDTH;
		}else if(b==scdn){
			float center=0.5f*(limits[1]+limits[0]);
			binsize*=1.5f;
			limits[1]=center+0.5f*binsize*WIDTH;
			limits[0]=limits[1]-binsize*WIDTH;
		}else if(b==edit){
			editStack();
			return;
		}else{
			float text=1.0f/Float.parseFloat(scale.getText());
			if(text!=binsize)
				binsize=text;
			float center=0.5f*(limits[1]+limits[0]);
			limits[1]=center+0.5f*binsize*WIDTH;
			limits[0]=limits[1]-binsize*WIDTH;
		}
		// check to make sure we haven't gone past the bounds
		if(limits[0]<0.0f){
			float diff=-limits[0];
			limits[0]=0.0f;
			limits[1]+=diff;
		}else{
			if(limits[1]>height2-1){
				float diff=limits[1]-(height2-1);
				limits[0]-=diff;
				limits[1]-=diff;
				// check again on the lower bound since it takes precedence
				if(limits[0]<0.0f){
					diff=-limits[0];
					limits[0]=0.0f;
					limits[1]+=diff;
				}
			}

		}
		scale.setText(""+1.0f/binsize);
		updateStack();
	}

	public float[] getLimits(){
		return limits;
	}

}
