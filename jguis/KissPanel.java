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
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.text.TextWindow;
import jalgs.interpolation;
import jalgs.jdataio;
import jalgs.jsort;
import jalgs.jstatistics;
import jalgs.jseg.findblobs3;
import jalgs.jseg.jsmooth;
import jalgs.jseg.measure_object;
import jalgs.jsim.rngs;

import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

public class KissPanel extends Panel implements ActionListener,ItemListener,MouseMotionListener,MouseListener{
	//this version allows for 7 channels: one extra plus dapi as a "channel"
	//it also implements a more versatile coloring scheme

	private static String defaultDirectory=null;
	public int nch=7; //number of channels
	public int nchrom=24; //number of chromosomes
	public float objthresh=0.35f; //threshold for whether a chromosome is "positive" for a dye
	public float[] objthresh2; //individual thresholds
	public float areathresh=30.0f; //threshold for whether a chromosome is the right size
	public boolean dapilast=true; //if last channel is really just dapi
	public boolean coloravg=false;
	public boolean dosmooth=false;
	public int smoothd=10;
	public int xsizerank=8;
	public int ysizerank=23;
	public int cellwidth=150;
	public int karyowidth;
	public int cellheight=100;
	public int chromspacing=3;
	public int karyoheight;
	public int[][] code;
	public int[][] positions;
	public boolean[] flip;
	private Button create_karyotype,save_data,guess_assign,out_stats,out_roi_contr,flip_button;
	private Checkbox mask_check;
	private Label idlabel,assignlabel,smoothlabel,smoothdlabel;
	private Label[] statlabels;
	private Choice assignchoice,smoothchoice;
	private Checkbox[] contrcheck;
	private Checkbox[] showcheck;
	private TextField smoothdval;
	private int nobjects,currobj;
	public ImagePlus threshimp,colorkissimp,karyoimp;
	public PlotWindow4 pw;
	public float[] unclust,karyounclust,objects,karyoobj,maxs,mins;
	public boolean[] iscluster;
	public float[][] object_stats,unmixed;
	private float[][] spec_stack; // here is the multispectral image
	private float[][][] spectra; // here are the spectra [dye][2][values]
	private int[] areas;
	private int[] arearank;
	private ImagePlus oldcolorkissimp;
	public boolean showmask;
	private boolean objsel;
	private findblobs3 fb,karyofb;
	private Image barimg;
	public static final int framewidth=300;
	public static final int frameheight=650;
	private int[] assignments;
	private int[][] contr;
	private Color[] colors;
	private boolean[] showcolors;
	public String[] names;
	public int[] colorindices;
	
	public static int[][] mousecode={
		{1,1,0,1,0,0,1,85}, //1
    	{1,0,0,0,1,0,1,69}, //2
    	{1,1,0,0,1,0,1,85}, //3
    	{0,0,1,1,0,0,1,100}, //4
    	{0,0,0,1,0,0,1,74}, //5
    	{1,1,1,1,0,0,1,83}, //6
    	{1,0,0,1,0,0,1,79}, //7
    	{0,0,0,1,1,0,1,74}, //8
    	{0,1,0,0,0,0,1,76}, //9
    	{0,1,0,0,1,0,1,61}, //10
    	{1,0,1,1,0,0,1,56}, //11
    	{0,1,0,1,0,0,1,75}, //12
    	{0,1,1,0,0,0,1,87}, //13
    	{1,1,0,0,0,0,1,63}, //14
    	{1,0,0,0,0,0,1,49}, //15
    	{0,0,0,0,1,0,1,58}, //16
    	{1,1,1,0,0,0,1,55}, //17
    	{0,1,1,1,0,0,1,66}, //18
    	{0,0,1,0,0,0,1,41}, //19
    	{1,0,1,0,1,0,1,63}, //X
    	{1,0,1,0,0,0,1,58} //Y
	};

	public static int[][] humancode={
		{1,1,0,1,0,0,1,92}, // 1
		{1,0,0,0,1,0,1,100}, // 2
		{1,1,0,0,1,0,1,86}, // 3
		{0,0,1,1,0,0,1,88}, // 4
		{0,0,0,1,0,0,1,70}, // 5
		{1,1,1,1,0,0,1,71}, // 6
		{1,0,0,1,0,0,1,68}, // 7
		{1,0,0,0,0,1,1,64}, // 8
		{1,1,1,0,1,0,1,55}, // 9
		{1,0,0,1,1,0,1,60}, // 10
		{0,1,0,0,0,0,1,59}, // 11
		{0,1,0,0,1,0,1,54}, // 12
		{1,0,1,1,0,0,1,48}, // 13
		{0,1,0,1,0,0,1,48}, // 14
		{0,0,0,1,1,0,1,50}, // 15
		{1,1,0,0,0,0,1,43}, // 16
		{1,1,0,1,1,0,1,48}, // 17
		{1,0,1,0,1,0,1,45}, // 18
		{0,0,0,0,1,0,1,36}, // 19
		{0,0,1,0,1,0,1,36}, // 20
		{0,1,1,1,0,0,1,24}, // 21
		{0,0,1,0,0,0,1,32}, // 22
		{1,0,1,1,1,0,1,75}, // x
		{1,0,1,0,0,0,1,39} // y
	};
	
	public static String[] alphabet={"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
	
	public static Color[] allcolors={Color.white,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};

	public static String[] colornames={"white","blue","green","red","magenta","cyan","yellow","orange"};

	public static Frame launch_frame(KissPanel panel){
		final Frame f=new Frame("Karyotype Identification by Spectral Separation");
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
		panel.setBounds(10,40,KissPanel.framewidth-20,KissPanel.frameheight-50);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(KissPanel.framewidth,KissPanel.frameheight));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}
	
	public static Object[] getCustomCodes(){
		//this method gets custom codes from the ImageJ-defaults folder
		//the returned object array contains the name String[], the mouse boolean[], and then the code int[][]
		String dir=System.getProperty("user.home");
		File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"edit_KISS_codes_jru_v1.jrn");
		if(!b.exists()) return null;
		List<List<String>> listtable=table_tools.getTableFromFile(b,"\t",false);
		if(listtable==null) return null;
		//table_tools.create_table("test import",listtable);
		listtable.remove(0); //remove the labels
		//the first line of the table contains the number of entries
		//IJ.log(listtable.get(0).get(0));
		int ncodes=table_tools.get_integer(listtable,0,0); //IJ.log(""+ncodes);
		if(ncodes==0) return null;
		String[] names=new String[ncodes];
		boolean[] mouse=new boolean[ncodes];
		//the first few lines of the table will contain the names of the codes and whether or not they are mouse
		for(int i=0;i<ncodes;i++){
			names[i]=listtable.get(i+1).get(0); //IJ.log(names[i]);
			mouse[i]=(table_tools.get_integer(listtable,i+1,1)==1); //IJ.log(""+mouse[i]);
		}
		Object[] retvals=new Object[2+ncodes];
		retvals[0]=names; retvals[1]=mouse;
		int counter=ncodes+1;
		for(int i=0;i<ncodes;i++){
			int[][] code=null;
			if(mouse[i]) code=new int[21][8];
			else code=new int[24][8];
			for(int j=0;j<code.length;j++){
				float[] temp=table_tools.get_row_array(listtable,counter);
				for(int k=0;k<8;k++) code[j][k]=(int)temp[k];
				counter++;
			}
			retvals[2+i]=code;
		}
		return retvals;
	}
	
	public static void setCustomCodes(Object[] data){
		List<List<String>> listtable=new ArrayList<List<String>>();
		List<String> temp=new ArrayList<String>();
		String[] names=(String[])data[0];
		if(names.length>=1){
    		boolean[] mouse=(boolean[])data[1];
    		temp.add(""+names.length); for(int i=0;i<7;i++) temp.add("0");
    		listtable.add(temp);
    		for(int i=0;i<names.length;i++){
    			int tempmouse=mouse[i]?1:0;
    			List<String> temp2=new ArrayList<String>();
    			temp2.add(names[i]); temp2.add(""+tempmouse); for(int j=0;j<6;j++) temp2.add("0");
    			listtable.add(temp2);
    		}
    		for(int i=0;i<names.length;i++){
    			int[][] code=(int[][])data[2+i];
    			for(int j=0;j<code.length;j++){
    				List<String> temp2=new ArrayList<String>();
    				for(int k=0;k<8;k++) temp2.add(""+code[j][k]);
    				listtable.add(temp2);
    			}
    		}
		} else {
			for(int i=0;i<8;i++) temp.add("0");
			listtable.add(temp);
		}
		//table_tools.create_table("test",listtable,null);
		String dir=System.getProperty("user.home");
		File a=new File(dir+File.separator+"ImageJ_defaults");
		if(!a.exists()) a.mkdir();
		String path=dir+File.separator+"ImageJ_defaults"+File.separator+"edit_KISS_codes_jru_v1.jrn";
		table_tools.writeTableToFile(path,null,listtable,0);
	}
	
	public void init(InputStream is){
		init(is,false);
	}

	public void init(InputStream is,boolean mouse){
		init(is,mouse,null);
	}
	
	public void init(InputStream is,boolean mouse,int[][] custcode){
		if(custcode!=null) code=custcode;
		else code=humancode;
		if(mouse){
			if(custcode!=null) code=custcode;
			else code=mousecode;
			nchrom=21;
		}
		init_names();
		jdataio jdio=new jdataio();
		int width=jdio.readintelint(is);
		int height=jdio.readintelint(is);
		nobjects=jdio.readintelint(is);
		int hasspectra=jdio.readintelint(is);
		showmask=(jdio.readintelint(is)==1);
		float[] threshpix=new float[width*height];
		float threshmin=jdio.readintelfloat(is);
		float threshmax=jdio.readintelfloat(is);
		jdio.readintelfloatfile(is,threshpix);
		objects=new float[width*height];
		jdio.readintelfloatfile(is,objects);
		unmixed=new float[nch][width*height];
		jdio.readintelfloatfile(is,unmixed);
		int[] cpixels=new int[width*height];
		jdio.readintelintfile(is,cpixels);
		object_stats=new float[nobjects][nch+1];
		jdio.readintelfloatfile(is,object_stats);
		contr=new int[nobjects][nch];
		jdio.readintelintfile(is,contr);
		assignments=new int[nobjects];
		jdio.readintelintfile(is,assignments);
		areas=new int[nobjects];
		arearank=new int[nobjects];
		jdio.readintelintfile(is,areas);
		jdio.readintelintfile(is,arearank);
		if(hasspectra==1){
			int spectrallength=jdio.readintelint(is);
			spectra=new float[nobjects][2][spectrallength];
			jdio.readintelfloatfile(is,spectra);
			spec_stack=new float[spectrallength][width*height];
			jdio.readintelfloatfile(is,spec_stack);
		}
		ImageStack dispstack=new ImageStack(width,height);
		dispstack.addSlice("",threshpix);
		dispstack.addSlice("",new float[width*height]);
		threshimp=new ImagePlus("Outlined Objects",dispstack);
		threshimp.setOpenAsHyperStack(true);
		threshimp.setDimensions(2,1,1);
		threshimp=new CompositeImage(threshimp,CompositeImage.COMPOSITE);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		graylut.min=threshmin;
		graylut.max=threshmax;
		((CompositeImage)threshimp).setChannelLut(graylut,1);
		((CompositeImage)threshimp).setDisplayRange(threshmin,threshmax);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		redlut.min=0.0;
		redlut.max=255.0;
		((CompositeImage)threshimp).setChannelLut(redlut,2);
		threshimp.show();
		fb=new findblobs3(width,height);
		fb.nobjects=nobjects;
		initClusters();
		update_image();
		int[] oldcpixels=new int[width*height];
		System.arraycopy(cpixels,0,oldcpixels,0,width*height);
		oldcolorkissimp=new ImagePlus("",new ColorProcessor(width,height,oldcpixels));
		colorkissimp=new ImagePlus("Kiss Image",new ColorProcessor(width,height,cpixels));
		colorkissimp.show();
		initui();
		objsel=false;
		update_karyotype();
		this.threshimp.getCanvas().addMouseMotionListener(this);
		colorkissimp.getCanvas().addMouseMotionListener(this);
		karyoimp.getCanvas().addMouseMotionListener(this);
		karyoimp.getCanvas().addMouseListener(this);
	}
	
	public void init(ImagePlus threshimp,ImagePlus kissimp,float[] objects,int[] areas,int[] arearank,findblobs3 fb,boolean showmask,float[][][] spectra,Object[] spec_stack){
		init(threshimp,kissimp,objects,areas,arearank,fb,showmask,spectra,spec_stack,false);
	}

	public void init(ImagePlus threshimp,ImagePlus kissimp,float[] objects,int[] areas,int[] arearank,findblobs3 fb,boolean showmask,float[][][] spectra,Object[] spec_stack,boolean mouse){
		init(threshimp,kissimp,objects,areas,arearank,fb,showmask,spectra,spec_stack,mouse,null);
	}
	
	public void init(ImagePlus threshimp,ImagePlus kissimp,float[] objects,int[] areas,int[] arearank,findblobs3 fb,boolean showmask,float[][][] spectra,Object[] spec_stack,boolean mouse,int[][] custcode){
		if(custcode!=null) code=custcode;
		else code=humancode;
		if(mouse){
			if(custcode!=null) code=custcode;
			else code=mousecode;
			nchrom=21;
		}
		this.showmask=showmask;
		this.threshimp=threshimp;
		this.objects=objects;
		this.nobjects=areas.length;
		this.areas=areas;
		this.arearank=arearank;
		this.fb=fb;
		initClusters();
		this.spectra=spectra;
		init_names();
		if(spec_stack!=null){
			this.spec_stack=new float[spec_stack.length][];
			for(int i=0;i<spec_stack.length;i++){
				this.spec_stack[i]=jutils.convert_arr_float(spec_stack[i]);
			}
		}
		if(kissimp!=null){
			unmixed=new float[nch][];
			ImageStack stack=kissimp.getStack();
			for(int i=0;i<(nch-1);i++) unmixed[i]=(float[])stack.getProcessor(i+1).convertToFloat().getPixelsCopy();
			if(dapilast) unmixed[nch-1]=(float[])threshimp.getStack().getProcessor(1).convertToFloat().getPixelsCopy();
			else unmixed[nch-1]=(float[])stack.getProcessor(nch).convertToFloat().getPixelsCopy();
			//truncate the negatives
			for(int i=0;i<nch;i++){
				for(int j=0;j<unmixed[0].length;j++){
					if(unmixed[i][j]<0.0f){
						unmixed[i][j]=0.0f;
					}
				}
			}
			get_obj_contr();
		}else{
			sim_obj_contr();
		}
		initui();
		objsel=false;
		reset_contr();
		reset_assignments();
		update_karyotype();
		this.threshimp.getCanvas().addMouseMotionListener(this);
		colorkissimp.getCanvas().addMouseMotionListener(this);
		karyoimp.getCanvas().addMouseMotionListener(this);
		karyoimp.getCanvas().addMouseListener(this);
	}

	private void initui(){
		setLayout(null);
		create_karyotype=new Button("Reset Karyotype");
		create_karyotype.setBounds(10,10,100,30);
		create_karyotype.addActionListener(this);
		add(create_karyotype);
		mask_check=new Checkbox("Show Mask",showmask);
		mask_check.setBounds(10,50,100,20);
		mask_check.addItemListener(this);
		add(mask_check);
		idlabel=new Label("object id = n/a");
		idlabel.setBounds(10,80,100,20);
		add(idlabel);
		statlabels=new Label[nch+1];
		contrcheck=new Checkbox[nch];
		showcheck=new Checkbox[nch];
		String[] dyelabels=new String[nch];
		showcolors=new boolean[nch+1];

		for(int i=0;i<nch;i++) dyelabels[i]=alphabet[i];
		for(int i=0;i<nch;i++){
			statlabels[i]=new Label("dye"+dyelabels[i]);
			statlabels[i].setBounds(10,110+i*30,30,20);
			add(statlabels[i]);
			contrcheck[i]=new Checkbox("",false);
			contrcheck[i].setBounds(160,110+i*30,20,20);
			contrcheck[i].addItemListener(this);
			add(contrcheck[i]);
			showcolors[i]=true;
			showcheck[i]=new Checkbox("Show",showcolors[i]);
			showcheck[i].setBounds(190,110+i*30,50,20);
			showcheck[i].addItemListener(this);
			add(showcheck[i]);
		}
		statlabels[nch]=new Label("area value = n/a");
		statlabels[nch].setBounds(10,110+nch*30,150,20);
		add(statlabels[nch]);
		assignlabel=new Label("Assignment");
		assignlabel.setBounds(10,110+(nch+1)*30,75,20);
		add(assignlabel);
		assignchoice=new Choice();
		assignchoice.setBounds(90,110+(nch+1)*30,75,20);
		for(int i=0;i<nchrom;i++){
			assignchoice.add(names[i]);
		}
		assignchoice.add("NA");
		assignchoice.addItemListener(this);
		add(assignchoice);
		save_data=new Button("Save");
		save_data.setBounds(10,110+(nch+2)*30,75,30);
		add(save_data);
		save_data.addActionListener(this);
		guess_assign=new Button("Guess");
		guess_assign.setBounds(10,110+(nch+3)*30,75,30);
		add(guess_assign);
		guess_assign.addActionListener(this);
		out_stats=new Button("Print Stats");
		out_stats.setBounds(10,110+(nch+4)*30,75,30);
		add(out_stats);
		out_stats.addActionListener(this);
		out_roi_contr=new Button("Roi Contr");
		out_roi_contr.setBounds(10,110+(nch+5)*30,75,30);
		add(out_roi_contr);
		out_roi_contr.addActionListener(this);
		smoothlabel=new Label("Smoothing");
		smoothlabel.setBounds(10,110+(nch+6)*30,75,20);
		add(smoothlabel);
		smoothchoice=new Choice();
		smoothchoice.setBounds(90,110+(nch+6)*30,75,20);
		smoothchoice.add("None");
		smoothchoice.add("Boxcar");
		smoothchoice.add("Show_Avg");
		smoothchoice.addItemListener(this);
		add(smoothchoice);
		smoothdlabel=new Label("Diameter");
		smoothdlabel.setBounds(10,110+(nch+7)*30,75,20);
		add(smoothdlabel);
		smoothd=10;
		smoothdval=new TextField(""+smoothd);
		smoothdval.setBounds(90,110+(nch+7)*30,75,20);
		add(smoothdval);
		flip_button=new Button("Flip");
		flip_button.setBounds(10,150+(nch+7)*30,75,30);
		add(flip_button);
		flip_button.addActionListener(this);
	}
	
	public void initClusters(){
		//this finds clusters and ids them
		iscluster=fb.find_clusters(objects);
		unclust=getUnclust(objects,fb);
	}
	
	public float[] getUnclust(float[] objects1,findblobs3 fb1){
		return (new findblobs3(fb1.width,fb1.height)).dofindblobs(objects1,0.5f);
	}

	public void paint(Graphics g){
		g.drawImage(barimg,50,110,this);
	}

	public void update(Graphics g){
		paint(g);
	}

	public void setVisible(boolean b){
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==create_karyotype){
			objsel=false;
			reset_contr();
			reset_assignments();
			update_karyotype();
			update_kiss();
			update_image();
		}
		if(e.getSource()==save_data){
			// here we save all of the data required to recreate this analysis
			FileDialog fd=new FileDialog((Frame)getParent(),"Save as Kiss Object...",FileDialog.SAVE);
			if(defaultDirectory!=null)
				fd.setDirectory(defaultDirectory);
			fd.setFile(".kiss");
			fd.show();
			String name=fd.getFile();
			String directory=fd.getDirectory();
			defaultDirectory=directory;
			fd.dispose();
			saveAsObject(directory+File.separator+name);
		}
		if(e.getSource()==guess_assign){
			if(objsel){
				float[] c2=new float[code.length];
				for(int i=0;i<code.length;i++){
					//compare the channels to the code
					for(int j=0;j<nch;j++){
						c2[i]+=(object_stats[currobj-1][j]-code[i][j])*(object_stats[currobj-1][j]-code[i][j]);
					}
					//here is the area comparison
					float areadif=object_stats[currobj-1][nch]-0.01f*code[i][code[0].length-1];
					c2[i]+=areadif*areadif;
				}
				int[] order=jsort.javasort_order(c2);
				GenericDialog gd=new GenericDialog("Choose Chromosome");
				String[] best=new String[5];
				for(int i=0;i<5;i++){
					best[i]=names[order[i]];
				}
				gd.addChoice("In order of probabilty",best,best[0]);
				gd.showDialog();
				if(gd.wasCanceled()){
					return;
				}
				int index1=gd.getNextChoiceIndex();
				int index=order[index1];
				if(index!=assignments[currobj-1]){
					assignments[currobj-1]=index;
					objsel=false;
					update_karyotype();
					update_kiss();
					update_image();
				}
			}
		}
		if(e.getSource()==flip_button){
			if(objsel){
				//IJ.log("flipped");
				if(flip[currobj-1]) flip[currobj-1]=false;
				else flip[currobj-1]=true;
				update_karyotype();
				//update_kiss();
				//update_image();
			}
		}
		if(e.getSource()==out_stats){
			String headings="";
			for(int i=0;i<nch;i++) headings+=alphabet[i]+"\t";
			headings+="Area\tAssignment";
			StringBuffer output=new StringBuffer();
			int[] areaorder=jsort.get_javasort_order(arearank);
			for(int i=0;i<nobjects;i++){
				for(int j=0;j<(nch+1);j++){
					output.append(object_stats[areaorder[i]][j]+"\t");
				}
				int tempassign=assignments[areaorder[i]];
				if(tempassign>=0){
					output.append(names[tempassign]+"\n");
				}else{
					output.append("Unassigned"+"\n");
				}
			}
			new TextWindow("Chromosome Stats",headings,output.toString(),400,200);
		}
		if(e.getSource()==out_roi_contr){
			get_roi_contr();
		}
	}

	public void saveAsObject(String filename){
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(filename));
			save2os(os);
			os.close();
		}catch(IOException e){
			return;
		}
	}

	public void save2os(OutputStream os){
		// output the width,height,nobjects,spectra flag,and showmask flag
		jdataio jdio=new jdataio();
		jdio.writeintelint(os,threshimp.getWidth());
		jdio.writeintelint(os,threshimp.getHeight());
		jdio.writeintelint(os,nobjects);
		jdio.writeintelint(os,(spectra!=null)?1:0);
		jdio.writeintelint(os,showmask?1:0);
		// now write the transmitted light image and its scalings
		ImageProcessor threship=threshimp.getStack().getProcessor(1);
		jdio.writeintelfloat(os,(float)threship.getMin());
		jdio.writeintelfloat(os,(float)threship.getMax());
		jdio.writeintelfloatarray(os,(float[])threship.getPixels());
		// now the objects image
		jdio.writeintelfloatarray(os,objects);
		// now write the unmixed array
		jdio.writeintelfloatarray(os,unmixed);
		// now write the kiss spread image
		jdio.writeintelintarray(os,(int[])colorkissimp.getProcessor().getPixels());
		// now write the object stats array
		jdio.writeintelfloatarray(os,object_stats);
		// now write the contributions array
		jdio.writeintelintarray(os,contr);
		// now write the assignments array
		jdio.writeintelintarray(os,assignments);
		// now write the area and arearank arrays
		jdio.writeintelintarray(os,areas);
		jdio.writeintelintarray(os,arearank);
		// now write the spectral array if not null
		// first write the spectral length
		if(spectra!=null){
			jdio.writeintelint(os,spectra[0][0].length);
			jdio.writeintelfloatarray(os,spectra);
			jdio.writeintelfloatarray(os,spec_stack);
		}
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==mask_check){
			showmask=mask_check.getState();
			update_image();
		}
		if(e.getSource()==smoothchoice){
			int index=smoothchoice.getSelectedIndex();
			coloravg=false;
			dosmooth=false;
			if(index==1){
				dosmooth=true;
				smoothd=(int)Float.parseFloat(smoothdval.getText());
				//smoothdval.
			}
			if(index==2) coloravg=true;
			regenerate_kiss();
			update_kiss();
			update_karyotype();
		}
		for(int i=0;i<nch;i++){
			boolean wasthis=false;
			if(e.getSource()==showcheck[i]){
				showcolors[i]=showcheck[i].getState();
				wasthis=true;
			}
			if(wasthis){
				regenerate_kiss();
				update_kiss();
				update_karyotype();
			}
		}
		if(objsel){
			boolean changed=false;
			for(int i=0;i<nch;i++){
				if(e.getSource()==contrcheck[i]){
					int newstate=contrcheck[i].getState()?1:0;
					if(newstate!=contr[currobj-1][i]){
						changed=true;
						contr[currobj-1][i]=newstate;
					}
				}
			}
			if(changed){
				update_assignment(currobj);
				objsel=false;
				update_karyotype();
				update_kiss();
				update_image();
			}
			if(e.getSource()==assignchoice){
				int index=assignchoice.getSelectedIndex();
				if(index==(nchrom)){
					index=-1;
				}
				if(index!=assignments[currobj-1]){
					assignments[currobj-1]=index;
					objsel=false;
					update_karyotype();
					update_kiss();
					update_image();
				}
			}
		}
	}

	public void mouseMoved(MouseEvent e){
		if(!objsel){
			ImageCanvas ic=threshimp.getCanvas();
			ImageCanvas ic2=colorkissimp.getCanvas();
			ImageCanvas ic3=karyoimp.getCanvas();
			if(e.getSource()==ic){
				int x=e.getX();
				int y=e.getY();
				int ox=ic.offScreenX(x);
				int oy=ic.offScreenY(y);
				if(ox<threshimp.getWidth()&&oy<threshimp.getHeight()){
					currobj=(int)objects[ox+threshimp.getWidth()*oy];
					update_object_stats(currobj);
				}
			}
			if(e.getSource()==ic2){
				int x=e.getX();
				int y=e.getY();
				int ox=ic2.offScreenX(x);
				int oy=ic2.offScreenY(y);
				if(ox<threshimp.getWidth()&&oy<threshimp.getHeight()){
					currobj=(int)objects[ox+threshimp.getWidth()*oy];
					update_object_stats(currobj);
				}
			}
			if(e.getSource()==ic3){
				int x=e.getX();
				int y=e.getY();
				int ox=ic3.offScreenX(x);
				int oy=ic3.offScreenY(y);
				if(ox<karyoimp.getWidth()&&oy<karyoimp.getHeight()){
					currobj=(int)karyoobj[ox+karyoimp.getWidth()*oy];
					update_object_stats(currobj);
				}
			}
			if(currobj>0){
				int tempassign=assignments[currobj-1];
				if(tempassign==-1){
					tempassign=nchrom;
				}
				assignchoice.select(tempassign);
			}
		}
	}

	public void mouseDragged(MouseEvent e){
	}

	public void mouseClicked(MouseEvent e){
		ImageCanvas ic=karyoimp.getCanvas();
		if(e.getSource()==ic){
			if(e.getClickCount()==2){
				int x=e.getX();
				int y=e.getY();
				int ox=ic.offScreenX(x);
				int oy=ic.offScreenY(y);
				if(ox<karyoimp.getWidth()&&oy<karyoimp.getHeight()){
					int tempobj=(int)karyoobj[ox+karyoimp.getWidth()*oy];
					if(tempobj>0){
						select_object(tempobj);
					}
				}
			}
		}
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
	}

	public void mouseReleased(MouseEvent e){
	}

	public void select_object(int id){
		// this method acts as a toggle to select and unselect objects and swap
		// object selections
		//karyounclust=getUnclust(karyoobj,karyofb); //maybe could do this only when moving objects
		if(objsel&&id==currobj){
			// here we are unselecting the current object
			objsel=false;
			karyoimp.setColor(Color.black);
			if(iscluster[currobj-1]){
				Polygon poly=karyofb.get_cluster_outline(karyounclust,karyoobj,currobj);
				jutils.draw_polygon(karyoimp.getProcessor(),poly,true);
			} else {
				Polygon poly=karyofb.get_object_outline(karyoobj,currobj); //breaks multiobject
				jutils.draw_polygon(karyoimp.getProcessor(),poly,true); //breaks multiobject
			}
			update_kiss();
		}else{
			// here we are either swapping selections or creating a new one
			if(objsel=true){
				karyoimp.setColor(Color.black);
				if(iscluster[currobj-1]){
					Polygon poly=karyofb.get_cluster_outline(karyounclust,karyoobj,currobj);
					jutils.draw_polygon(karyoimp.getProcessor(),poly,true);
				} else {
					Polygon poly=karyofb.get_object_outline(karyoobj,currobj); //breaks multiobject
					jutils.draw_polygon(karyoimp.getProcessor(),poly,true); //breaks multiobject
				}
			}
			currobj=id;
			objsel=true;
			update_object_stats(currobj);
			karyoimp.setColor(Color.white);
			if(iscluster[currobj-1]){
				Polygon poly=karyofb.get_cluster_outline(karyounclust,karyoobj,currobj);
				jutils.draw_polygon(karyoimp.getProcessor(),poly,true);
			} else {
				Polygon poly=karyofb.get_object_outline(karyoobj,currobj); //breaks multiobject
				jutils.draw_polygon(karyoimp.getProcessor(),poly,true); //breaks multiobject
			}
		}
		karyoimp.updateAndDraw();
		update_kiss();
		update_image();
		update_spectra();
		int tempassign=assignments[currobj-1];
		if(tempassign==-1){
			tempassign=nchrom;
		}
		assignchoice.select(tempassign);
	}

	public void update_object_stats(int id){
		if(id>0){
			int sortid=arearank[id-1]+1;
			idlabel.setText("id = "+id+" rank = "+sortid);
			ColorProcessor cp=new ColorProcessor(100,nch*30-10);
			barimg=cp.createImage();
			Graphics g=barimg.getGraphics();
			g.setColor(Color.white);
			g.fillRect(0,0,100,nch*30-10);
			g.setColor(Color.red);
			for(int i=0;i<nch;i++){
				g.fillRect(0,i*30,(int)(object_stats[id-1][i]*100.0f),20);
				if(contr[id-1][i]==1)
					contrcheck[i].setState(true);
				else
					contrcheck[i].setState(false);
			}
			g.setColor(Color.black);
			g.drawRect(0,0,99,nch*30-10-1);
			statlabels[nch].setText("area value = "+(int)(object_stats[id-1][nch]*100.0));
			repaint();
		}
	}

	public void reset_contr(){
		contr=new int[nobjects][nch];
		if(objthresh2==null){
			for(int i=0;i<nch;i++) objthresh2[i]=objthresh;
		}
		for(int i=0;i<nobjects;i++){
			for(int j=0;j<nch;j++){
				contr[i][j]=(object_stats[i][j]>objthresh2[j])?1:0;
			}
		}
	}

	public void reset_assignments(){
		assignments=new int[nobjects];
		for(int i=0;i<nobjects;i++){
			update_assignment(i+1);
		}
	}

	public void update_assignment(int id){
		assignments[id-1]=-1;
		for(int j=0;j<code.length;j++){
			boolean matching=true;
			for(int k=0;k<nch;k++){
				if(contr[id-1][k]!=code[j][k]){
					matching=false;
					break;
				}
			}
			if(matching){ //check if the area is within 30% of the code area
				if(Math.abs(object_stats[id-1][nch]*100.0f-code[j][code[0].length-1])>areathresh){
					matching=false;
				}
			}
			if(matching){
				assignments[id-1]=j;
				break;
			}
		}
	}

	public void get_obj_contr(){
		object_stats=new float[nobjects][nch+1];
		maxs=new float[nch+1];
		mins=new float[nch+1];
		float[][] pixels=unmixed;
		for(int i=0;i<nch;i++){
			float[] temp2=new float[nobjects];
			for(int j=0;j<nobjects;j++){
				object_stats[j][i]=fb.get_object_stats(objects,j+1,pixels[i],"Avg");
				if(object_stats[j][i]>maxs[i]){
					maxs[i]=object_stats[j][i];
				}
				if(i==0){
					object_stats[j][nch]=areas[j];
					if(object_stats[j][nch]>maxs[nch]){
						maxs[nch]=object_stats[j][nch];
					}
				}
				temp2[j]=object_stats[j][i];
			}
			jsort.javasort_order(temp2);
			maxs[i]=0.25f*(temp2[nobjects-1]+temp2[nobjects-2]+temp2[nobjects-3]+temp2[nobjects-4]);
			mins[i]=0.25f*(temp2[0]+temp2[1]+temp2[2]+temp2[3]);
			if(dapilast && i==(nch-1)) mins[i]=0.0f; //this is a dapi specific construct (there are no dapi negative chromosomes)
		}
		// now normalize all values to the maximum value
		float[] objmax=new float[nobjects];
		for(int i=0;i<(nch+1);i++){
			for(int j=0;j<nobjects;j++){
				object_stats[j][i]=(object_stats[j][i]-mins[i])/(maxs[i]-mins[i]);
				if(i<nch && object_stats[j][i]>objmax[j]){
					objmax[j]=object_stats[j][i];
				}
			}
		}
		// renormalize colors within each chromosome
		for(int i=0;i<nobjects;i++){
			for(int j=0;j<nch;j++){
				object_stats[i][j]/=objmax[i];
			}
		}
		regenerate_kiss();
	}

	public void get_roi_contr(){
		// IJ.log("getting Roi");
		Roi roi=colorkissimp.getRoi();
		if(roi==null){
			roi=threshimp.getRoi();
			if(roi==null)
				return;
		}
		Polygon poly=roi.getPolygon();
		float[] roicontr=new float[nch];
		float tempmax=0.0f;
		// IJ.log("getting contributions");
		if(mins==null)
			get_mins_maxs();
		for(int i=0;i<nch;i++){
			roicontr[i]=jstatistics.getstatistic("Avg",unmixed[i],colorkissimp.getWidth(),colorkissimp.getHeight(),poly,null);
			roicontr[i]=(roicontr[i]-mins[i])/(maxs[i]-mins[i]);
			if(roicontr[i]>tempmax)
				tempmax=roicontr[i];
		}
		for(int i=0;i<nch;i++)
			roicontr[i]/=tempmax;
		// IJ.log("drawing plot");
		ColorProcessor cp=new ColorProcessor(100,nch*30-10);
		Image roicontrimg=cp.createImage();
		Graphics g=roicontrimg.getGraphics();
		g.setColor(Color.white);
		g.fillRect(0,0,100,nch*30-10);
		g.setColor(Color.red);
		for(int i=0;i<nch;i++){
			g.fillRect(0,i*30,(int)(roicontr[i]*100.0f),20);
		}
		g.setColor(Color.black);
		g.drawRect(0,0,99,nch*30-10-1);
		cp=new ColorProcessor(roicontrimg);
		(new ImagePlus("Roi Contributions",cp)).show();
	}

	private void get_mins_maxs(){
		maxs=new float[nch+1];
		mins=new float[nch+1];
		for(int i=0;i<nch;i++){
			float[] temp2=new float[nobjects];
			for(int j=0;j<nobjects;j++){
				temp2[j]=fb.get_object_stats(objects,j+1,unmixed[i],"Avg");
			}
			jsort.javasort_order(temp2);
			maxs[i]=0.25f*(temp2[nobjects-1]+temp2[nobjects-2]+temp2[nobjects-3]+temp2[nobjects-4]);
			mins[i]=0.25f*(temp2[0]+temp2[1]+temp2[2]+temp2[3]);
			if(dapilast && i==(nch-1)) mins[i]=0.0f;
		}
	}
	
	private float get_biggest_avg(float[] sorted,int nbiggest,boolean biggest){
		if(biggest){
			float avg=0.0f;
			for(int i=0;i<nbiggest;i++) avg+=sorted[sorted.length-1-i];
			return avg/nbiggest;
		} else {
			float avg=0.0f;
			for(int i=0;i<nbiggest;i++) avg+=sorted[i];
			return avg/nbiggest;
		}
	}

	public void sim_obj_contr(){
		//note that this is human specific
		int width=threshimp.getWidth();
		int height=threshimp.getHeight();
		float[] multiplier=(float[])threshimp.getStack().getPixels(1);
		float[][] pixels=new float[nch][width*height];
		double stdev=100.0;
		rngs random=new rngs();
		float fback=0.2f;
		ImageStack stack=new ImageStack(width,height);
		for(int i=0;i<nch;i++){
			for(int j=0;j<width*height;j++){
				if(objects[j]>0.0f){
					int rank=arearank[(int)objects[j]-1];
					int chromosome=0;
					if(rank<xsizerank){
						chromosome=(int)(0.5f*rank);
					}else{
						if(rank>=xsizerank){
							if(rank==xsizerank){
								//chromosome=22;
								chromosome=nchrom-2;
							}else{
								if(rank<ysizerank){
									chromosome=(int)(0.5f*(rank-1));
								}else{
									if(rank==ysizerank){
										//chromosome=23;
										chromosome=nchrom-1;
									}else{
										chromosome=(int)(0.5f*(rank-2));
									}
								}
							}
						}
					}
					if(code[chromosome][i]==0){
						pixels[i][j]=fback*(float)random.gasdev(multiplier[j],stdev);
					}else{
						pixels[i][j]=(float)random.gasdev(multiplier[j],stdev);
					}
				}
			}
			stack.addSlice("",pixels[i]);
		}
		unmixed=pixels;
		ImagePlus tempimp=new ImagePlus("Simulated Kiss Image",stack);
		tempimp.show();
		get_obj_contr();
	}

	public void update_image(){
		ImageProcessor ip=threshimp.getStack().getProcessor(2);
		float[] disppix=(float[])ip.getPixels();
		if(showmask){
			float[] temp=(float[])jutils.convert_array(fb.tobinary(objects,true),2);
			System.arraycopy(temp,0,disppix,0,disppix.length);
			if(objsel){	
				Polygon poly=null;
				if(iscluster[currobj-1]){
					poly=fb.get_cluster_outline(unclust,objects,currobj);//breaks multiobject
				} else {
					poly=fb.get_object_outline(objects,currobj);//breaks multiobject
				}
				ip.setColor(Color.black);
				jutils.fill_polygon(ip,poly);//breaks multiobject
				ip.setColor(Color.white);
				jutils.draw_polygon(ip,poly,true);//breaks multiobject
			}
		}else{
			for(int i=0;i<disppix.length;i++){
				disppix[i]=0.0f;
			}
			Polygon[] objects2=fb.get_object_outlines(objects);//breaks multiobject
			for(int i=0;i<nobjects;i++){
				if(iscluster[i]){
					objects2[i]=fb.get_cluster_outline(unclust,objects,i+1);
					//IJ.log("polygon output");
					//for(int j=0;j<objects2[i].npoints;j++) IJ.log(objects2[i].xpoints[j]+",\t"+objects2[i].ypoints[j]);
				}
			}
			ip.setColor(Color.white);
			for(int i=0;i<objects2.length;i++){
				jutils.draw_polygon(ip,objects2[i],true);//breaks multiobject
			}
			if(objsel){
				Polygon poly=null;
				if(iscluster[currobj-1]){
					poly=fb.get_cluster_outline(unclust,objects,currobj);//breaks multiobject
				} else {
					poly=fb.get_object_outline(objects,currobj);//breaks multiobject
				}
				ip.setColor(Color.black);
				jutils.draw_polygon(ip,poly,true);//breaks multiobject
				ip.setColor(Color.white);
				jutils.fill_polygon(ip,poly);//breaks multiobject
			}
		}
		threshimp.updateAndDraw();
		// nobjects.setText("# of objects = "+fb.nobjects);
	}

	public void regenerate_kiss(){
		//here we generate the false colored kiss images using DAPI as our coloring mask
		float[][] temp_obj_stats=new float[nobjects][nch];
		int width=threshimp.getWidth();
		int height=threshimp.getHeight();
		float[][] pixels=new float[nch+1][];
		for(int i=0;i<nch;i++)
			pixels[i]=unmixed[i].clone();
		pixels[nch]=((float[])threshimp.getStack().getPixels(1)).clone();
		for(int i=0;i<nch;i++){
			float[] temp=new float[nobjects];
			for(int j=0;j<nobjects;j++){
				temp_obj_stats[j][i]=fb.get_object_stats(objects,j+1,pixels[i],"Avg");
				temp[j]=temp_obj_stats[j][i];
			}
		}
		//here is an option to smooth the objects
		//the neighboring regions are not smoothed
		//IJ.log("smoothing");
		if(dosmooth){
    		for(int i=0;i<nch;i++){
    			for(int j=0;j<width*height;j++){
    				if(objects[j]==0.0f) pixels[i][j]=0.0f;
    			}
    			pixels[i]=jsmooth.smooth2Dobjects(pixels[i],width,height,smoothd,objects);
    		}
    		//IJ.log("smoothd = "+smoothd);
    		//new ImagePlus("test2",jutils.array2stack(pixels,width,height)).show();
		}
		if(coloravg){
		//here we can replace the objects with their average
    		for(int i=0;i<nch;i++){
    			for(int j=0;j<width*height;j++){
    				if(objects[j]==0.0f) pixels[i][j]=0.0f;
    				else pixels[i][j]=temp_obj_stats[(int)objects[j]-1][i];
    			}
    		}
		}
		//need to figure out a way to average the color but retain the shading
		// now normalize all values to the maximum value
		//IJ.log("find obj max's");
		if(mins==null || maxs==null) get_mins_maxs();
		float[] objmax=new float[nobjects];
		for(int i=0;i<nch;i++){
			for(int j=0;j<nobjects;j++){
				temp_obj_stats[j][i]=(temp_obj_stats[j][i]-mins[i])/(maxs[i]-mins[i]);
				if(temp_obj_stats[j][i]>objmax[j]){
					objmax[j]=temp_obj_stats[j][i];
				}
			}
		}
		// now create a false color image based on the normalized contributions
		int[] cpixels=new int[width*height];
		//by default show only the first 5 channels
		//IJ.log("set up colors");
		if(colors==null){
			colors=new Color[nch+1]; 
			for(int i=0;i<colors.length;i++) colors[i]=allcolors[0];
			if(colorindices==null){colorindices=new int[]{4,1,2,6,3};}
			for(int i=0;i<colorindices.length;i++) colors[i]=allcolors[colorindices[i]];
		}
		if(showcolors==null){
			showcolors=new boolean[nch+1];
			for(int i=0;i<nch;i++) showcolors[i]=true;
			showcolors[nch]=false;
		}
		int[][] lutmaxs=new int[3][nch+1];
		for(int i=0;i<(nch+1);i++){
			int[] rgb=jutils.get_color_RGB(colors[i]);
			lutmaxs[0][i]=rgb[0];
			lutmaxs[1][i]=rgb[1];
			lutmaxs[2][i]=rgb[2];
		}
		//IJ.log("test1");
		for(int i=0;i<width*height;i++){
			if(objects[i]>0.0f){
				int[] rgb=new int[3];
				float tempobjmax=objmax[(int)objects[i]-1];
				for(int j=0;j<(nch+1);j++){
					if(showcolors[j]){
						float percentmax=0.4f*(pixels[j][i]-mins[j])/(maxs[j]-mins[j]);
						//float percentmax=temp_obj_stats[(int)objects[i]-1][j]; //this would give uniform shading
						percentmax/=tempobjmax;
						for(int k=0;k<3;k++){
							rgb[k]+=(int)(lutmaxs[k][j]*percentmax);
						}
					}
				}
				if(rgb[0]>255) rgb[0]=255;
				if(rgb[1]>255) rgb[1]=255;
				if(rgb[2]>255) rgb[2]=255;
				if(rgb[0]<0) rgb[0]=0;
				if(rgb[1]<0) rgb[1]=0;
				if(rgb[2]<0) rgb[2]=0;
				cpixels[i]=jutils.rgb2intval(rgb[0],rgb[1],rgb[2]);
			}
		}
		//IJ.log("test2");
		int[] oldcpixels=cpixels.clone();
		//new ImagePlus("test image",new ColorProcessor(width,height,cpixels.clone())).show();
		oldcolorkissimp=new ImagePlus("",new ColorProcessor(width,height,oldcpixels));
		if(colorkissimp==null){
			colorkissimp=new ImagePlus("Kiss Image",new ColorProcessor(width,height,cpixels));
			colorkissimp.show();
		}else{
			colorkissimp.setProcessor(new ColorProcessor(width,height,cpixels));
		}
	}

	public void update_kiss(){
		int[] cpixels=(int[])oldcolorkissimp.getProcessor().getPixelsCopy();
		ImageProcessor cp=colorkissimp.getProcessor();
		cp.setPixels(cpixels);
		if(objsel){
			Polygon poly=null;
			if(iscluster[currobj-1]){
				poly=fb.get_cluster_outline(unclust,objects,currobj);
			} else {
				poly=fb.get_object_outline(objects,currobj);//breaks multiobject
			}
			cp.setColor(Color.white);
			jutils.draw_polygon(cp,poly,true);//breaks multiobject
		}
		colorkissimp.updateAndDraw();
	}

	public void update_spectra(){
		if(spectra!=null){
			if(pw==null){
				pw=new PlotWindow4("Chromosome Spectrum","Spectral Unit","Intensity",spectra[currobj-1][0]);
				pw.draw();
				pw.addPoints(spectra[currobj-1][1],true);
			}else{
				pw.updateSeries(spectra[currobj-1][0],0,true);
				pw.updateSeries(spectra[currobj-1][1],1,true);
			}
		}
	}

	public void init_names(){
		names=new String[code.length];
		for(int i=0;i<code.length;i++){
			StringBuffer sb=new StringBuffer();
			if(i<(code.length-2)){
				sb.append(""+(i+1));
			}else{
				if(i==(code.length-2))
					sb.append("X");
				else
					sb.append("Y");
			}
			sb.append("(");
			for(int j=0;j<nch;j++){
				if(code[i][j]!=0)
					sb.append(""+alphabet[j]);
			}
			sb.append(")");
			names[i]=sb.toString();
		}
	}
	
	public void init_positions(boolean mouse){
		if(mouse){
			//mouse chromosomes are arranged 5 x 4 with the y chromosome added at bottom right
			karyowidth=cellwidth*6;
			karyoheight=cellheight*5;
			positions=new int[nchrom+1][2];
			int counter=0;
			for(int i=0;i<4;i++){
				for(int j=0;j<5;j++){
					positions[counter][0]=j*cellwidth+chromspacing;
					positions[counter][1]=i*cellheight+chromspacing;
					counter++;
				}
			}
			//now the y chromosome
			positions[counter][0]=5*cellwidth+chromspacing;
			positions[counter][1]=3*cellheight+chromspacing;
			counter++;
			//and the unassigned row
			positions[counter][0]=chromspacing;
			positions[counter][1]=4*cellheight+chromspacing;
		} else {
			//human chromosomes are arranged 7 x 4 with two gaps between 3 and 4 on the first row,
			//one gap between 15 and 16 on the third row,
			//and a half gap between 20 and 21 and between 22 and X on the fourth row
			karyowidth=cellwidth*7;
			karyoheight=cellheight*5;
			positions=new int[nchrom+1][2];
			int counter=0;
			//do the first row with two spaces
			for(int i=0;i<3;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=chromspacing;
				counter++;
			}
			for(int i=5;i<7;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=chromspacing;
				counter++;
			}
			//the second is complete
			for(int i=0;i<7;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=cellheight+chromspacing;
				counter++;
			}
			//the third has one space
			for(int i=0;i<3;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=2*cellheight+chromspacing;
				counter++;
			}
			for(int i=4;i<7;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=2*cellheight+chromspacing;
				counter++;
			}
			//finally two half spaces
			for(int i=0;i<2;i++){
				positions[counter][0]=i*cellwidth+chromspacing;
				positions[counter][1]=3*cellheight+chromspacing;
				counter++;
			}
			for(int i=3;i<5;i++){
				positions[counter][0]=i*cellwidth+chromspacing-cellwidth/2;
				positions[counter][1]=3*cellheight+chromspacing;
				counter++;
			}
			for(int i=6;i<8;i++){
				positions[counter][0]=i*cellwidth+chromspacing-cellwidth;
				positions[counter][1]=3*cellheight+chromspacing;
				counter++;
			}
			//and then the unassigned
			positions[counter][0]=chromspacing;
			positions[counter][1]=4*cellheight+chromspacing;
		}
	}

	public void update_karyotype(){
		if(positions==null){
			init_positions(nchrom==21);
		}
		if(karyoimp==null){
			karyoimp=new ImagePlus("Karyotype",new ColorProcessor(karyowidth,karyoheight));
			karyoimp.show();
		}
		if(flip==null){
			flip=new boolean[nobjects];
		}
		jutils.clearColorImp(karyoimp);
		Polygon[] objects2=fb.get_object_outlines(objects);//breaks multiobject
		for(int i=0;i<nobjects;i++){
			if(iscluster[i]){
				objects2[i]=fb.get_cluster_outline(unclust,objects,i+1);
			}
		}
		karyoobj=new float[karyowidth*karyoheight];
		karyounclust=new float[karyowidth*karyoheight];
		karyoimp.setColor(Color.white);
		for(int i=0;i<nchrom;i++){
			int xposition=positions[i][0];
			int yposition=positions[i][1];
			karyoimp.getProcessor().drawString(names[i],xposition,yposition-2*chromspacing+cellheight);
			karyoimp.getProcessor().drawRect(xposition-chromspacing,yposition-chromspacing,cellwidth,cellheight);
			for(int k=0;k<nobjects;k++){
				if(assignments[k]==i){
					Polygon rotroi=copy_rotate_object(oldcolorkissimp,objects2[k],karyoimp,xposition,yposition,(k+1),karyoobj,flip[k]);//breaks multiobject
					Rectangle r=rotroi.getBounds();//breaks multiobject
					xposition+=r.width+chromspacing;
				}
			}
		}
		int xposition=positions[nchrom][0];
		int yposition=positions[nchrom][1];
		karyoimp.getProcessor().drawString("Unassigned",xposition,yposition-2*chromspacing+cellheight);
		karyoimp.getProcessor().drawRect(xposition-chromspacing,yposition-chromspacing,karyowidth,cellheight);
		for(int i=0;i<nobjects;i++){
			if(assignments[i]<0){
				Polygon rotroi=copy_rotate_object(oldcolorkissimp,objects2[i],karyoimp,xposition,yposition,(i+1),karyoobj,flip[i]);//breaks multiobject
				Rectangle r=rotroi.getBounds();//breaks multiobject
				xposition+=r.width+chromspacing;
			}
		}
		karyofb=new findblobs3(karyowidth,karyoheight);
		karyounclust=getUnclust(karyoobj,karyofb);
		karyoimp.updateAndDraw();
	}

	public Polygon copy_rotate_object(ImagePlus source,Polygon outline,ImagePlus destination,int x,int y,int objid,float[] objimg,boolean flip1){
		//breaks multiobject
		// start by cropping the image and the profile
		Rectangle r=outline.getBounds();
		int bigger=r.width;
		if(r.height>bigger)
			bigger=r.height;
		int xstart=r.x-(2*bigger-r.width)/2;
		int width=2*bigger;
		int ystart=r.y-(2*bigger-r.height)/2;
		int height=2*bigger;
		int[] newxpts=new int[outline.npoints];
		int[] newypts=new int[outline.npoints];
		for(int i=0;i<outline.npoints;i++){
			newxpts[i]=outline.xpoints[i]-xstart;
			newypts[i]=outline.ypoints[i]-ystart;
		}
		int[] newimage=new int[width*height];
		int[] pixels=(int[])source.getProcessor().getPixels();
		int oldwidth=source.getWidth();
		for(int i=ystart;i<(ystart+height);i++){
			for(int j=xstart;j<(xstart+width);j++){
				if(outline.contains(j,i)){
					newimage[j-xstart+(i-ystart)*width]=pixels[j+i*oldwidth];
				}
			}
		}
		// new ImagePlus("test",new
		// ColorProcessor(width,height,newimage)).show();
		Polygon newoutline=new Polygon(newxpts,newypts,outline.npoints);
		// jutils.copyRoi(newoutline,newimage,width,height,destination,x,y,null);
		// return newoutline;
		// now get the ellipse parameters
		float[] ellipse=measure_object.get_ellipse_parameters(newoutline);
		float rotangle=0.5f*(float)Math.PI-ellipse[2];
		if(flip1) rotangle+=(float)Math.PI;
		// IJ.log(""+ellipse[2]);
		int[] rotimage2=interpolation.rotate_color_image(newimage,width,height,-rotangle,(int)ellipse[0],(int)ellipse[1]);
		Polygon rotroi=interpolation.rotate_polygon(newoutline,rotangle,(int)ellipse[0],(int)ellipse[1]);
		jutils.copyColorRoi(rotroi,rotimage2,width,height,destination,x,y,null);
		int dstwidth=destination.getWidth();
		int dstheight=destination.getHeight();
		Rectangle rr=rotroi.getBounds();
		for(int i=0;i<rr.height;i++){
			for(int j=0;j<rr.width;j++){
				if((x+j)>=0&&(x+j)<dstwidth&&(y+i)>=0&&(y+i)<dstheight){
					if(rotroi.contains(j+rr.x,i+rr.y)){
						objimg[j+x+(i+y)*dstwidth]=objid;
					}
				}
			}
		}
		return rotroi;
	}

	private void separate_objects(){
		// this is a bit painful from the programming standpoint but necessary
		// we need to separate two chromosomes
		// we need to completely redo object stats but don't change any
		// contributions or assignments except for the separated objects
		// will have to remap everything from old object index to new
		// object must be selected and have a linear (or polylinear) roi to
		// complete action
		// start by separating the object mask
		if(objsel){
			float[] newobj=objects.clone();
			Roi roi=colorkissimp.getRoi();
			if(roi==null){
				roi=threshimp.getRoi();
				if(roi==null)
					return;
			}
			Polygon poly;
			if(roi instanceof Line){
				int[] xpts={((Line)roi).x1,((Line)roi).x2};
				int[] ypts={((Line)roi).y1,((Line)roi).y2};
				// IJ.log(""+xpts[0]+" , "+ypts[0]+" , "+xpts[1]+" , "+ypts[1]);
				poly=new Polygon(xpts,ypts,2);
			}else{
				poly=roi.getPolygon();
			}
			fb.separateobjects(newobj,poly,false,currobj);
			int newnobjects=fb.nobjects;
			// create a map that relates each new object to the old one
			int[] newobjmap=new int[newnobjects];
			for(int i=0;i<newobj.length;i++){
				if(newobj[i]>0.0f&&newobjmap[(int)newobj[i]-1]==0){
					newobjmap[(int)newobj[i]-1]=(int)objects[i];
				}
			}
			// now regenerate all object stats

		}
	}

}
