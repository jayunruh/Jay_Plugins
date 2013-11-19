/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import j3D.*;

import jalgs.jsim.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import ij.*;
import ij.gui.*;
import ij.process.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class maxproj3D_panel extends JPanel implements ActionListener,ItemListener,ImageListener{

	private JButton rotright,rotleft,rotup,rotdown,rotclock,rotcounter,reset,makemovie;
	private Choice methodchoice;
	private Label threshlabel,methodlabel,xysteplabel;
	private TextField threshval,xystepval;
	public ImagePlus imp3D,srcimp;
	public ImageStack stack,stack3D;
	public int method,width,height,slices,channels,frames,newzslices,maxsize,imagetype;
	public int currframe,currchan,xystep,nthreads,xstart1,xend1,ystart1,yend1,zstart1,zend1;
	public float zratio,xrot,yrot,zrot,invzr,halfzr;
	private float[] distances,thresh;
	private float[][] images;
	private rngs random;
	private element3D[] elements;
	private line3D[] lines;
	private renderer r,r2;

	public static Frame launch_frame(maxproj3D_panel panel){
		final Frame f=new Frame("3D Viewer");
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
		panel.setBounds(10,40,440,200);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(450,250));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(ImagePlus imp,float zratio,float thresh){
		init(imp,zratio,thresh,1);
	}

	public void init(ImagePlus imp,float zratio,float thresh,int nthreads){
		setLayout(null);
		this.zratio=zratio;
		this.nthreads=nthreads;
		xystep=1;
		srcimp=imp;
		int bsize=55;
		// 55
		int bheight=30;
		Font bfont=new Font("Lucida Sans Unicode",Font.BOLD,12);
		rotup=new JButton("\u2191");
		rotup.setFont(bfont);
		rotup.setBounds(10,10,bsize,bheight);
		rotup.addActionListener(this);
		add(rotup);
		rotdown=new JButton("\u2193");
		rotdown.setFont(bfont);
		rotdown.setBounds(10+bsize,10,bsize,bheight);
		rotdown.addActionListener(this);
		add(rotdown);
		rotleft=new JButton("\u2190");
		rotleft.setFont(bfont);
		rotleft.setBounds(10+2*bsize,10,bsize,bheight);
		rotleft.addActionListener(this);
		add(rotleft);
		rotright=new JButton("\u2192");
		rotright.setFont(bfont);
		rotright.setBounds(10+3*bsize,10,bsize,bheight);
		rotright.addActionListener(this);
		add(rotright);
		rotclock=new JButton("\u21b6");
		rotclock.setFont(bfont);
		rotclock.setBounds(10+4*bsize,10,bsize,bheight);
		rotclock.addActionListener(this);
		add(rotclock);
		rotcounter=new JButton("\u21b7");
		rotcounter.setFont(bfont);
		rotcounter.setBounds(10+5*bsize,10,bsize,bheight);
		rotcounter.addActionListener(this);
		add(rotcounter);
		reset=new JButton("Reset");
		reset.setFont(bfont);
		reset.setBounds(10+6*bsize,10,75,bheight);
		reset.addActionListener(this);
		add(reset);
		methodlabel=new Label("Interpolation");
		methodlabel.setBounds(10,bheight+20,80,20);
		add(methodlabel);
		methodchoice=new Choice();
		methodchoice.add("None");
		methodchoice.add("Random");
		methodchoice.add("Step");
		methodchoice.add("Linear");
		methodchoice.setBounds(90,bheight+20,100,20);
		methodchoice.addItemListener(this);
		add(methodchoice);
		threshlabel=new Label("Threshold");
		threshlabel.setBounds(10,bheight+50,80,20);
		add(threshlabel);
		threshval=new TextField(""+thresh);
		threshval.setBounds(90,bheight+50,60,20);
		threshval.addActionListener(this);
		add(threshval);
		xysteplabel=new Label("Resample");
		xysteplabel.setBounds(160,bheight+50,80,20);
		add(xysteplabel);
		xystepval=new TextField(""+xystep);
		xystepval.setBounds(250,bheight+50,60,20);
		xystepval.addActionListener(this);
		add(xystepval);
		makemovie=new JButton("Make Movie");
		makemovie.setFont(bfont);
		makemovie.setBounds(10,bheight+80,150,bheight);
		makemovie.addActionListener(this);
		add(makemovie);
		width=imp.getWidth();
		height=imp.getHeight();
		stack=imp.getStack();
		slices=imp.getNSlices();
		channels=imp.getNChannels();
		frames=imp.getNFrames();
		currframe=imp.getFrame()-1;
		currchan=imp.getChannel()-1;
		xstart1=0;
		xend1=width-1;
		ystart1=0;
		yend1=height-1;
		zstart1=0;
		zend1=slices;
		this.thresh=new float[channels];
		for(int i=0;i<channels;i++)
			this.thresh[i]=thresh;
		newzslices=(int)(zratio*(float)slices);
		maxsize=(int)Math.sqrt(width*width+height*height+newzslices*newzslices);
		int xrem=(maxsize-width)/2;
		int yrem=(maxsize-height)/2;
		int zrem=(maxsize-newzslices)/2;
		lines=new line3D[width*height];
		elements=new element3D[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				lines[j+i*width]=new line3D(j+xrem,i+yrem,zrem,j+xrem,i+yrem,newzslices+zrem,Color.black);
				elements[j+i*width]=lines[j+i*width];
			}
		}
		r=new renderer(maxsize,maxsize);
		r.centerz=newzslices/2+zrem;
		r2=new renderer(maxsize,maxsize);
		r2.centerz=r.centerz;
		element3D[] cube=new element3D[1];
		cube[0]=new rect3D(r.centerx,r.centery,r.centerz,width,height,newzslices,Color.white);
		r2.setelementarray(cube);
		r2.background=Color.black;
		r.setelementarray(elements);
		stack3D=new ImageStack(maxsize,maxsize);
		images=new float[channels+1][maxsize*maxsize];
		for(int i=0;i<(channels+1);i++)
			stack3D.addSlice("",images[i]);
		imp3D=new ImagePlus("3D Rendering",stack3D);
		imp3D.setOpenAsHyperStack(true);
		imp3D.setDimensions(channels+1,1,1);
		imp3D=new CompositeImage(imp3D,CompositeImage.COMPOSITE);
		imp3D.copyScale(imp);
		if(imp.isComposite()){
			LUT[] luts=((CompositeImage)imp).getLuts();
			LUT[] luts2=new LUT[channels+1];
			System.arraycopy(luts,0,luts2,0,channels);
			luts2[channels]=LUT.createLutFromColor(Color.white);
			((CompositeImage)imp3D).setLuts(luts2);
		}else{
			((CompositeImage)imp3D).setChannelLut(LUT.createLutFromColor(Color.white),channels+1);
		}
		imp3D.show();
		r.set_perspective(maxsize*2.0);
		xrot=0.0f;
		yrot=0.0f;
		zrot=0.0f;
		init_distances();
		method=0;
		random=new rngs();
		imagetype=0;
		if(stack.getPixels(1) instanceof short[])
			imagetype=1;
		if(stack.getPixels(1) instanceof float[])
			imagetype=2;
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		r2.setrotation((int)xrot,(int)yrot,(int)zrot);
		ImagePlus.addImageListener(this);
		update_image();
		((CompositeImage)imp3D).resetDisplayRanges();
	}

	public void setVisible(boolean b){
		ImagePlus.removeImageListener(this);
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==makemovie){
			make_movie();
			return;
		}
		if(e.getSource()==rotup)
			xrot+=10.0f;
		if(e.getSource()==rotdown)
			xrot-=10.0f;
		if(e.getSource()==rotright)
			yrot+=10.0f;
		if(e.getSource()==rotleft)
			yrot-=10.0f;
		if(e.getSource()==rotclock)
			zrot+=10.0f;
		if(e.getSource()==rotcounter)
			zrot-=10.0f;
		if(e.getSource()==reset){
			xrot=0.0f;
			yrot=0.0f;
			zrot=0.0f;
		}
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		r2.setrotation((int)xrot,(int)yrot,(int)zrot);
		float temp=Float.parseFloat(threshval.getText());
		thresh[currchan]=temp;
		xystep=Integer.parseInt(xystepval.getText());
		update_image();
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==methodchoice){
			method=methodchoice.getSelectedIndex();
			update_image();
		}
	}

	public void imageClosed(ImagePlus imp){
	}

	public void imageOpened(ImagePlus imp){
	}

	public void imageUpdated(ImagePlus imp){
		if(imp.equals(imp3D)){
			int tempcurrchan=imp3D.getChannel()-1;
			if(tempcurrchan!=currchan&&tempcurrchan<channels){
				currchan=tempcurrchan;
				srcimp.setPosition(currchan+1,srcimp.getSlice(),currframe+1);
				threshval.setText(""+thresh[currchan]);
			}
		}else{
			if(imp.equals(srcimp)){
				int tempcurrchan=srcimp.getChannel()-1;
				if(tempcurrchan!=currchan){
					currchan=tempcurrchan;
					imp3D.setPosition(currchan+1,1,1);
					threshval.setText(""+thresh[currchan]);
				}
				int tempcurrframe=srcimp.getFrame()-1;
				if(tempcurrframe!=currframe){
					currframe=tempcurrframe;
					update_image();
				}
			}
		}
	}

	public void init_distances(){
		distances=new float[maxsize*maxsize];
		for(int i=0;i<width*height;i++){
			float distx=(lines[i].pt2.rx-lines[i].pt1.rx);
			float disty=(lines[i].pt2.ry-lines[i].pt1.ry);
			float distz=(lines[i].pt2.rz-lines[i].pt1.rz);
			distances[i]=(float)Math.sqrt(distx*distx+disty*disty+distz*distz);
		}
	}

	public void update_image(){
		float[][] images2=maxprojimage();
		for(int i=0;i<channels;i++){
			System.arraycopy(images2[i],0,images[i],0,maxsize*maxsize);
		}
		float[] grid=r2.renderfloat(255.0f);
		System.arraycopy(grid,0,images[channels],0,maxsize*maxsize);
		imp3D.updateAndDraw();
	}

	public float interpolatez(int index,float z,int channel){
		float scaledz=z/zratio;
		int zprev=(int)scaledz;
		int znext=zprev+1;
		float zrem=scaledz-(float)zprev;
		float z1=get_value(index,zprev,channel);
		float z2=get_value(index,znext,channel);
		return z1+zrem*(z2-z1);
	}

	public float get_value(int index,int z,int channel){
		int index2=currframe*slices*channels+z*channels+channel;
		return get_value(index,index2);
	}

	public float get_value(int index,int index2){
		if(imagetype==0){
			return (float)(((byte[])stack.getPixels(index2+1))[index]&0xff);
		}else{
			if(imagetype==1){
				return (float)(((short[])stack.getPixels(index2+1))[index]&0xffff);
			}else{
				return ((float[])stack.getPixels(index2+1))[index];
			}
		}
	}

	public float[] get_line(int index,int index21){
		float[] temp=new float[slices];
		int index2=index21;
		for(int i=0;i<slices;i++){
			temp[i]=get_value(index,index2);
			index2+=channels;
		}
		return temp;
	}

	public void put_value(float value,float[] image,int xpoint,int ypoint){
		int index=xpoint+ypoint*maxsize;
		if(image[index]<value){
			if(xystep==1)
				image[index]=value;
			else{
				for(int i=0;i<xystep;i++){
					int temp=i*maxsize;
					for(int j=0;j<xystep;j++){
						image[index+j+temp]=value;
					}
				}
			}
		}
	}

	public void put_line(float[] inline,float thresh2,float[] outimage,float xstart,float ystart,float xinc,float yinc){
		float x=xstart;
		float y=ystart;
		float z=0.0f;
		// here xinc and yinc are the x and y steps corresponding to a z pixel
		// unit
		if(method==0){
			// here we just draw a single value per slice
			for(int k=0;k<slices;k++){
				float val=inline[k];
				if(val>thresh2){
					int xpoint=(int)(x+halfzr*xinc);
					int ypoint=(int)(y+halfzr*yinc);
					if(xpoint>=0&&xpoint<maxsize&&ypoint>=0&&ypoint<maxsize){
						put_value(val,outimage,xpoint,ypoint);
					}
				}
				x+=xinc*zratio;
				y+=yinc*zratio;
			}
		}
		if(method==1){
			// here we randomly choose a z value per slice (takes twice as long
			// as method 0)
			for(int k=0;k<slices;k++){
				float val=inline[k];
				if(val>thresh2){
					int xpoint=(int)random.unidev(x+zratio*xinc,x);
					int ypoint=(int)random.unidev(y+zratio*yinc,y);
					if(xpoint>=0&&xpoint<maxsize&&ypoint>=0&&ypoint<maxsize){
						put_value(val,outimage,xpoint,ypoint);
					}
				}
				x+=zratio*xinc;
				y+=zratio*yinc;
			}
		}
		if(method==2){
			// here we draw a line for each z segment (takes 4 times longer than
			// method 0)
			float temp=0.0f;
			for(int k=0;k<slices;k++){
				float val=inline[k];
				if(val>thresh2){
					for(int j=0;j<(int)zratio;j++){
						z=temp+(float)j;
						x=xstart+xinc*z;
						y=ystart+yinc*z;
						// int xpoint=Math.round(x);
						// int ypoint=Math.round(y);
						int xpoint=random.random_integerize(x);
						int ypoint=random.random_integerize(y);
						if(xpoint>=0&&xpoint<maxsize&&ypoint>=0&&ypoint<maxsize){
							put_value(val,outimage,xpoint,ypoint);
						}
					}
				}
				temp+=zratio;
			}
		}
		if(method==3){
			// here we linear interpolate for each z segment (takes 10 times
			// longer than method 0)
			float temp=0.0f;
			float oldval=inline[0];
			for(int k=0;k<(slices-1);k++){
				float newval=inline[k+1];
				if(oldval>thresh2||newval>thresh2){
					float val=oldval;
					float inc=invzr*(newval-oldval);
					for(int j=0;j<(int)zratio;j++){
						z=temp+(float)j;
						x=xstart+xinc*z;
						y=ystart+yinc*z;
						if(val>thresh2){
							int xpoint=random.random_integerize(x);
							int ypoint=random.random_integerize(y);
							// int xpoint=Math.round(x);
							// int ypoint=Math.round(y);
							if(xpoint>=0&&xpoint<maxsize&&ypoint>=0&&ypoint<maxsize){
								put_value(val,outimage,xpoint,ypoint);
							}
						}
						val+=inc;
					}
				}
				oldval=newval;
				temp+=zratio;
			}
		}
	}

	public void make_movie(){
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_steps",36,0);
		gd.addNumericField("X_Angle_Step_(deg)",10.0,5,15,null);
		gd.addNumericField("Y_Angle_Step_(deg)",10.0,5,15,null);
		gd.addNumericField("Z_Angle_Step_(deg)",0.0,5,15,null);
		gd.addNumericField("Frame_Step",0,0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		int nsteps=(int)gd.getNextNumber();
		float dx=(float)gd.getNextNumber();
		float dy=(float)gd.getNextNumber();
		float dz=(float)gd.getNextNumber();
		int dt=(int)gd.getNextNumber();
		ImageStack movstack=new ImageStack(maxsize,maxsize);
		float oldxrot=xrot;
		float oldyrot=yrot;
		float oldzrot=zrot;
		int oldt=currframe;
		for(int i=0;i<nsteps;i++){
			float[][] temp=maxprojimage();
			for(int j=0;j<temp.length;j++){
				movstack.addSlice("",temp[j]);
			}
			movstack.addSlice("",r2.renderfloat(255.0f));
			xrot+=dx;
			yrot+=dy;
			zrot+=dz;
			currframe+=dt;
			if(currframe>=frames)
				currframe=0;
			r.setrotation(xrot,yrot,zrot);
			r2.setrotation(xrot,yrot,zrot);
			IJ.showProgress(i,nsteps);
			IJ.showStatus("frame "+(i+1)+" of "+nsteps+" rendered");
		}
		xrot=oldxrot;
		yrot=oldyrot;
		zrot=oldzrot;
		currframe=oldt;
		r.setrotation((int)xrot,(int)yrot,(int)zrot);
		r2.setrotation((int)xrot,(int)yrot,(int)zrot);
		ImagePlus imp=new ImagePlus("Max Proj Movie",movstack);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(channels+1,1,nsteps);
		imp=new CompositeImage(imp,CompositeImage.COMPOSITE);
		imp.copyScale(imp3D);
		((CompositeImage)imp).copyLuts(imp3D);
		imp.show();
	}

	public float[][] maxprojimage(){
		float[][] outimage=new float[channels][maxsize*maxsize];
		int[] intparams={method,slices,maxsize,xystep};
		halfzr=0.5f*zratio;
		invzr=1.0f/zratio;
		float[] fltparams={zratio,halfzr,invzr};
		ExecutorService executor=null;
		if(nthreads>1)
			executor=Executors.newFixedThreadPool(nthreads);
		for(int ch=0;ch<channels;ch++){
			int tempindex=currframe*slices*channels+ch;
			for(int i=0;i<height;i+=xystep){
				for(int l=0;l<width;l+=xystep){
					int index=i*width+l;
					float distx=(lines[index].pt2.rx-lines[index].pt1.rx);
					float disty=(lines[index].pt2.ry-lines[index].pt1.ry);
					float dist=distances[index];
					float xunit2=distx/dist;
					float yunit2=disty/dist;
					// float xunit2=xunit*zratio;
					// float yunit2=yunit*zratio;
					float xstart=lines[index].pt1.rx;
					float ystart=lines[index].pt1.ry;
					float[] line=get_line(index,tempindex);
					if(nthreads>1){
						Runnable worker=new interp_line(intparams,fltparams,outimage[ch],line,xstart,ystart,xunit2,yunit2,thresh[ch]);
						executor.execute(worker);
					}else{
						put_line(line,thresh[ch],outimage[ch],xstart,ystart,xunit2,yunit2);
					}
				}
			}
		}
		if(nthreads>1){
			if(IJ.escapePressed()) executor.shutdownNow();
			else{
				executor.shutdown();
				// Wait until all threads are finished
				while(!executor.isTerminated()){}
			}
		}
		return outimage;
	}

}

class interp_line implements Runnable{

	private final float xstart,ystart,xinc,yinc,thresh;
	private final int[] intparams; // intparams is method,slices,maxsize,xystep
	private final float[] fltparams; // fltparams is zratio,halfzr,invzr
	private final float[] line;
	rngs random;
	public float[] outimage;

	public interp_line(int[] intparams,float[] fltparams,float[] outimage,float[] line,float xstart,float ystart,float xinc,float yinc,float thresh){
		this.intparams=intparams;
		this.fltparams=fltparams;
		this.xstart=xstart;
		this.ystart=ystart;
		this.xinc=xinc;
		this.yinc=yinc;
		this.thresh=thresh;
		this.line=line;
		random=new rngs();
		this.outimage=outimage;
	}

	public void run(){
		put_line();
	}

	public void put_line(){
		float x=xstart;
		float y=ystart;
		float z=0.0f;
		if(intparams[0]==0){
			// here we just draw a single value per slice
			for(int k=0;k<intparams[1];k++){
				float val=line[k];
				if(val>thresh){
					int xpoint=(int)(x+xinc*fltparams[1]);
					int ypoint=(int)(y+yinc*fltparams[1]);
					if(xpoint>=0&&xpoint<intparams[2]&&ypoint>=0&&ypoint<intparams[2]){
						put_value(val,outimage,xpoint,ypoint);
					}
				}
				x+=xinc*fltparams[0];
				y+=yinc*fltparams[0];
			}
		}
		if(intparams[0]==1){
			// here we randomly choose a z value per slice (takes twice as long
			// as method 0)
			for(int k=0;k<intparams[1];k++){
				float val=line[k];
				if(val>thresh){
					int xpoint=(int)random.unidev(x+xinc*fltparams[0],x);
					int ypoint=(int)random.unidev(y+yinc*fltparams[0],y);
					if(xpoint>=0&&xpoint<intparams[2]&&ypoint>=0&&ypoint<intparams[2]){
						put_value(val,outimage,xpoint,ypoint);
					}
				}
				x+=xinc*fltparams[0];
				y+=yinc*fltparams[0];
			}
		}
		if(intparams[0]==2){
			// here we draw a line for each z segment (takes 4 times longer than
			// method 0)
			float temp=0.0f;
			for(int k=0;k<intparams[1];k++){
				float val=line[k];
				if(val>thresh){
					for(int j=0;j<(int)fltparams[0];j++){
						z=temp+(float)j;
						x=xstart+xinc*z;
						y=ystart+yinc*z;
						// int xpoint=Math.round(x);
						// int ypoint=Math.round(y);
						int xpoint=random.random_integerize(x);
						int ypoint=random.random_integerize(y);
						if(xpoint>=0&&xpoint<intparams[2]&&ypoint>=0&&ypoint<intparams[2]){
							put_value(val,outimage,xpoint,ypoint);
						}
					}
				}
				temp+=fltparams[0];
			}
		}
		if(intparams[0]==3){
			// here we linear interpolate for each z segment (takes 10 times
			// longer than method 0)
			float temp=0.0f;
			float oldval=line[0];
			for(int k=0;k<(intparams[1]-1);k++){
				float newval=line[k+1];
				if(oldval>thresh||newval>thresh){
					float val=oldval;
					float inc=fltparams[2]*(newval-oldval);
					for(int j=0;j<(int)fltparams[0];j++){
						z=temp+(float)j;
						x=xstart+xinc*z;
						y=ystart+yinc*z;
						if(val>thresh){
							int xpoint=random.random_integerize(x);
							int ypoint=random.random_integerize(y);
							// int xpoint=Math.round(x);
							// int ypoint=Math.round(y);
							if(xpoint>=0&&xpoint<intparams[2]&&ypoint>=0&&ypoint<intparams[2]){
								put_value(val,outimage,xpoint,ypoint);
							}
						}
						val+=inc;
					}
				}
				oldval=newval;
				temp+=fltparams[0];
			}
		}
	}

	public void put_value(float value,float[] image,int xpoint,int ypoint){
		int index=xpoint+ypoint*intparams[2];
		if(image[index]<value){
			if(intparams[3]==1)
				image[index]=value;
			else{
				for(int i=0;i<intparams[3];i++){
					int temp=i*intparams[2];
					for(int j=0;j<intparams[3];j++){
						image[index+j+temp]=value;
					}
				}
			}
		}
	}
}
