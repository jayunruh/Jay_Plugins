/***************************
 * this is a collection of stitching tools giving similar results to the Grid/Collection stitching tools in Fiji
 * in principle this collection gives you more control over curation of coordinates
 * eventually I hope to implement a STICS/PIV method (optic flow in overlapped regions) of stitching and merge it into this workflow
 * the general sequence is this:
 * for coordinate optimization, select the appropriate data (e.g. max projection, selected slices, selected channels, etc)
 * initialize the coordinates with one of the getTileCoords utilities
 * create a list of pair indices with one of the getPairs utilities
 * then create an instance of this class with those coordinates
 * calculate the phase correlations with phaseCorrPad
 * find the peaks and their stats with getCorrPeaks (will need to rearrange stats from phasecorr,x,y,crosscorr to x,y,phasecorr,crosscorr)
 * pick the best peak
 * check and curate the coordinates with checkCoords
 * find the image to fix (needs to have good correlation)
 * preAlign the data set
 * get the average overlap with getAvgOverlap
 * (optional) stitch with stitch_frame
 * once coordinates are verified, stitch all slices and channels (and frames)
 */
package jalgs;

import java.util.ArrayList;
import java.util.List;

import jalgs.jfft.crosscorr2D;
import jalgs.jfft.fftutils;
import jalgs.jfft.padding;
import jalgs.jfft.po4realfft2D;
import jalgs.jseg.jsmooth;


public class stitching{

	public int width,height,typeindex,newwidth,newheight;
	public float xmin,xmax,ymin,ymax,foverlap,hoverlap,voverlap;
	crosscorr2D ccclass;
	public int fftdim;
	public float[] xvals,yvals;
	public boolean feather;
	public gui_interface gui;
	
	public stitching(int width,int height,float[] xvals,float[] yvals,int typeindex,gui_interface gui){
		this.width=width;
		this.height=height;
		this.typeindex=typeindex;
		xmin=jstatistics.getstatistic("Min",xvals,null);
		xmax=jstatistics.getstatistic("Max",xvals,null);
		ymin=jstatistics.getstatistic("Min",yvals,null);
		ymax=jstatistics.getstatistic("Max",yvals,null);
		xmax+=width;
		ymax+=height;
		newwidth=1+(int)(xmax-xmin);
		newheight=1+(int)(ymax-ymin);
		this.xvals=xvals;
		this.yvals=yvals;
		this.gui=gui;
	}
	
	public stitching(int width,int height,float[] xvals1,float[] yvals1,float psize,int typeindex,gui_interface gui){
		//here we use calibrated x and y values
		this.width=width;
		this.height=height;
		this.typeindex=typeindex;
		xvals=new float[xvals1.length];
		yvals=new float[yvals1.length];
		//convert the calibrated values to pixel values
		for(int i=0;i<xvals1.length;i++){
			xvals[i]=xvals1[i]/psize;
			yvals[i]=yvals1[i]/psize;
		}
		xmin=jstatistics.getstatistic("Min",xvals,null);
		xmax=jstatistics.getstatistic("Max",xvals,null);
		ymin=jstatistics.getstatistic("Min",yvals,null);
		ymax=jstatistics.getstatistic("Max",yvals,null);
		xmax+=width;
		ymax+=height;
		newwidth=1+(int)(xmax-xmin);
		newheight=1+(int)(ymax-ymin);
		this.gui=gui;
	}
	
	/************
	 * sets the global coordinate list for this class (used by stitch_frame), will convert coords to pixels
	 * @param xvals
	 * @param yvals
	 * @param psize: if coordinates are in pixels, set to 1
	 */
	public void setCoords(float[] xvals,float[] yvals,float psize){
		this.xvals=new float[xvals.length];
		this.yvals=new float[yvals.length];
		for(int i=0;i<xvals.length;i++){
			this.xvals[i]=xvals[i]/psize;
			this.yvals[i]=yvals[i]/psize;
		}
		xmin=jstatistics.getstatistic("Min",this.xvals,null);
		xmax=jstatistics.getstatistic("Max",this.xvals,null);
		ymin=jstatistics.getstatistic("Min",this.yvals,null);
		ymax=jstatistics.getstatistic("Max",this.yvals,null);
		xmax+=width;
		ymax+=height;
		newwidth=1+(int)(xmax-xmin);
		newheight=1+(int)(ymax-ymin);
	}
	
	/***********
	 * sets the global coordinate list for this class (used by stitch_frame), will convert coords to pixels
	 * @param coords: a list of coordinates
	 * @param psize: if coordinates are in pixels, set to 1
	 */
	public void setCoordsTranspose(float[][] coords,float psize){
		xvals=new float[coords.length];
		yvals=new float[coords.length];
		for(int i=0;i<xvals.length;i++){
			xvals[i]=coords[i][0]/psize;
			yvals[i]=coords[i][1]/psize;
		}
		xmin=jstatistics.getstatistic("Min",xvals,null);
		xmax=jstatistics.getstatistic("Max",xvals,null);
		ymin=jstatistics.getstatistic("Min",yvals,null);
		ymax=jstatistics.getstatistic("Max",yvals,null);
		xmax+=width;
		ymax+=height;
		newwidth=1+(int)(xmax-xmin);
		newheight=1+(int)(ymax-ymin);
	}
	
	/********************
	 * gets coordinates for a snake by rows starting at bottom left scanning pattern
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param width: image width
	 * @param height: image height
	 * @param overlap: fractional overlap between images
	 * @return a list of coordinates
	 */
	public static float[][] getTileCoords(int ximgs,int yimgs,int width,int height,float overlap){
		//this version implements snake by rows starting at bottom left
		return getTileCoords(ximgs,yimgs,width,height,overlap,0,yimgs-1,1,-1);
	}
	
	/********************
	 * gets coordinates for a general snake by rows scanning pattern
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param width: image width
	 * @param height: image height
	 * @param overlap: fractional overlap between images
	 * @param startxpos: starting x coordinate
	 * @param startypos: starting y coordinate
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @return a list of coordinates
	 */
	public static float[][] getTileCoords(int ximgs,int yimgs,int width,int height,float overlap,int startxpos,int startypos,int startxdir,int startydir){
		//assume order is snake by rows
		float[][] coords=new float[2][ximgs*yimgs];
		float imgwidth=width;
		float imgheight=height;
		float xspacing=imgwidth*(1.0f-overlap);
		float yspacing=imgheight*(1.0f-overlap);
		int ypos=startypos; //initial y position
		int xpos=startxpos; //initial x position
		int counter=0;
		int xdir=startxdir; //initial scan direction (pos = right)
		int ydir=startydir; //initial y scan direction (pos = right);
		for(int i=0;i<yimgs;i++){
			xpos=(xdir>0)?0:(ximgs-1);
			for(int j=0;j<ximgs;j++){
				coords[0][counter]=xspacing*xpos;
				coords[1][counter]=yspacing*ypos;
				xpos+=xdir; counter++;
			}
			xdir*=-1;
			ypos+=ydir;
		}
		return coords;
	}
	
	/*****************
	 * gets coordinates for a general snake by rows scanning pattern (scaled version)
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param width: image width
	 * @param height: image height
	 * @param overlap: fractional overlap between images
	 * @param startxpos: starting x coordinate
	 * @param startypos: starting y coordinate
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @param psize: pixel size
	 * @return a list of coordinates
	 */
	public static float[][] getTileCoords(int ximgs,int yimgs,int width,int height,float overlap,int startxpos,int startypos,int startxdir,int startydir,float psize){
		float[][] coords=getTileCoords(ximgs,yimgs,width,height,overlap,startxpos,startypos,startxdir,startydir);
		if(psize==1.0f) return coords;
		for(int i=0;i<coords.length;i++){
			for(int j=0;j<coords[i].length;j++){
				coords[i][j]*=psize;
			}
		}
		return coords;
	}
	
	/*******************
	 * gets a list of pairs for snake by rows scanning configuration
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @param diagonal: include diagonal pairs?
	 * @return a list of pair indices
	 */
	public static int[][] getPairs(int ximgs,int yimgs,int startxdir,int startydir,boolean diagonal){
		//corners have two pairs
		//first and last rows and columns have 3 pairs
		//all others have 4
		//start by generating the relative snake positions
		int xdir=startxdir;
		int ydir=startydir;
		int xpos=0;
		int ypos=0;
		int counter=0;
		int[][] positions=new int[ximgs*yimgs][];
		for(int i=0;i<yimgs;i++){
			xpos=(xdir>0)?0:(ximgs-1);
			for(int j=0;j<ximgs;j++){
				positions[counter]=new int[]{xpos,ypos};
				xpos+=xdir; counter++;
			}
			xdir*=-1;
			ypos+=ydir;
		}
		List<int[]> pairs=new ArrayList<int[]>();
		//if we do this right, we should never get redundant pairs
		for(int i=0;i<positions.length;i++){
			for(int j=(i+1);j<positions.length;j++){
				if(Math.abs(positions[j][0]-positions[i][0])==1 && Math.abs(positions[j][1]-positions[i][1])==0){
					//horizontal pairs
					pairs.add(new int[]{i,j});
				}
				if(Math.abs(positions[j][0]-positions[i][0])==0 && Math.abs(positions[j][1]-positions[i][1])==1){
					//vertical pairs
					pairs.add(new int[]{i,j});
				}
				if(diagonal){
    				if(Math.abs(positions[j][0]-positions[i][0])==1 && Math.abs(positions[j][1]-positions[i][1])==1){
    					//diagonal
    					pairs.add(new int[]{i,j});
    				}
				}
			}
		}
		int[][] pairs2=new int[pairs.size()][];
		for(int i=0;i<pairs.size();i++) pairs2[i]=pairs.get(i);
		return pairs2;
	}
	
	/*******************
	 * gets a list of pairs for the raster scanning configuration
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @param diagonal: include diagonal pairs?
	 * @return a list of pair indices
	 */
	public static int[][] getPairsRaster(int ximgs,int yimgs,int startxdir,int startydir,boolean diagonal){
		//corners have two pairs
		//first and last rows and columns have 3 pairs
		//all others have 4
		//start by generating the relative snake positions
		int xdir=startxdir;
		int ydir=startydir;
		int xpos=0;
		int ypos=0;
		int counter=0;
		int[][] positions=new int[ximgs*yimgs][];
		for(int i=0;i<yimgs;i++){
			xpos=(xdir>0)?0:(ximgs-1);
			for(int j=0;j<ximgs;j++){
				positions[counter]=new int[]{xpos,ypos};
				xpos+=xdir; counter++;
			}
			//xdir*=-1;
			ypos+=ydir;
		}
		List<int[]> pairs=new ArrayList<int[]>();
		//if we do this right, we should never get redundant pairs
		for(int i=0;i<positions.length;i++){
			for(int j=(i+1);j<positions.length;j++){
				if(Math.abs(positions[j][0]-positions[i][0])==1 && Math.abs(positions[j][1]-positions[i][1])==0){
					//horizontal pairs
					pairs.add(new int[]{i,j});
				}
				if(Math.abs(positions[j][0]-positions[i][0])==0 && Math.abs(positions[j][1]-positions[i][1])==1){
					//vertical pairs
					pairs.add(new int[]{i,j});
				}
				if(diagonal){
    				if(Math.abs(positions[j][0]-positions[i][0])==1 && Math.abs(positions[j][1]-positions[i][1])==1){
    					//diagonal
    					pairs.add(new int[]{i,j});
    				}
				}
			}
		}
		int[][] pairs2=new int[pairs.size()][];
		for(int i=0;i<pairs.size();i++) pairs2[i]=pairs.get(i);
		return pairs2;
	}
	
	/******************
	 * gets a general list of image pairs for arbitrary coordinates
	 * @param coords: a list of image coordinates
	 * @param width: image width
	 * @param height: image height
	 * @param diagonal: include diagonal pairs?
	 * @return a list of pair indices
	 */
	public static float[][] getPairs(float[][] coords,int width,int height,boolean diagonal){
		//here we retroactively get the pairs and their distances from a set of coordinates
		List<float[]> pairs=new ArrayList<float[]>();
		//if we do this right, we should never get redundant pairs
		for(int i=0;i<coords.length;i++){
			for(int j=(i+1);j<coords.length;j++){
				float xshift=coords[j][0]-coords[i][0];
				float yshift=coords[j][1]-coords[i][1];
				float absxshift=Math.abs(xshift);
				float absyshift=Math.abs(yshift);
				float relabsxshift=absxshift/(float)width;
				float relabsyshift=absyshift/(float)height;
				int rndxshift=0; if(relabsxshift>0.5f && relabsxshift<1.0f) rndxshift=1; if(relabsxshift>=1.0f) rndxshift=2;
				int rndyshift=0; if(relabsyshift>0.5f && relabsyshift<1.0f) rndyshift=1; if(relabsyshift>=1.0f) rndyshift=2;
				if(rndxshift==1 || rndyshift==1){
					//we have a neighbor, now where is it?
    				if(rndxshift==1 && rndyshift==0){
    					//horizontal pairs
    					pairs.add(new float[]{i,j,xshift,yshift});
    				}
    				if(rndxshift==0 && rndyshift==1){
    					//vertical pairs
    					pairs.add(new float[]{i,j,xshift,yshift});
    				}
    				if(diagonal){
        				if(rndxshift==1 && rndyshift==1){
        					//diagonal
        					pairs.add(new float[]{i,j,xshift,yshift});
        				}
    				}
				}
			}
		}
		float[][] pairs2=new float[pairs.size()][];
		for(int i=0;i<pairs.size();i++) pairs2[i]=pairs.get(i);
		return pairs2;
	}
	
	/***********************
	 * gets the average vertical and horizontal overlap for a set of coordinates
	 * @param coords: the coordinates
	 * @param width: image width
	 * @param height: image height
	 * @return the horizontal and vertical overlap
	 */
	public static float[] getAvgOverlap(float[][] coords,float width,float height){
		float hover=0.0f; int nhpairs=0;
		float vover=0.0f; int nvpairs=0;
		float halfwidth=0.5f*(float)width;
		float halfheight=0.5f*(float)height;
		for(int i=0;i<coords[0].length;i++){
			for(int j=(i+1);j<coords[0].length;j++){
				float xdist=Math.abs(coords[0][j]-coords[0][i]);
				float ydist=Math.abs(coords[1][j]-coords[1][i]);
				if(xdist>halfwidth && xdist<(float)width && ydist<halfheight){
					hover+=xdist; nhpairs++;
				}
				if(xdist<halfwidth && ydist>halfheight && ydist<(float)height){
					vover+=ydist; nvpairs++;
				}
			}
		}
		if(nhpairs>0) hover/=(float)nhpairs;
		else hover=(float)height;
		if(nvpairs>0) vover/=(float)nvpairs;
		else vover=(float)width;
		return new float[]{(float)width-hover,(float)height-vover};
	}
	
	/********************
	 * gets coordinates for a general raster by rows scanning pattern
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param width: image width
	 * @param height: image height
	 * @param overlap: fractional overlap between images
	 * @param startxpos: starting x coordinate
	 * @param startypos: starting y coordinate
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @return a list of coordinates
	 */
	public static float[][] getTileCoordsRaster(int ximgs,int yimgs,int width,int height,float overlap,int startxpos,int startypos,int startxdir,int startydir){
		//assume order is raster by rows
		float[][] coords=new float[2][ximgs*yimgs];
		float imgwidth=width;
		float imgheight=height;
		float xspacing=imgwidth*(1.0f-overlap);
		float yspacing=imgheight*(1.0f-overlap);
		int ypos=startypos; //initial y position
		int xpos=startxpos; //initial x position
		int counter=0;
		int xdir=startxdir; //initial scan direction (pos = right)
		int ydir=startydir; //initial y scan direction (pos = right);
		for(int i=0;i<yimgs;i++){
			xpos=(xdir>0)?0:(ximgs-1);
			for(int j=0;j<ximgs;j++){
				coords[0][counter]=xspacing*xpos;
				coords[1][counter]=yspacing*ypos;
				xpos+=xdir; counter++;
			}
			//xdir*=-1;
			ypos+=ydir;
		}
		return coords;
	}
	
	/**************************
	 * gets coordinates for a general raster by rows scanning pattern (scaled this time)
	 * @param ximgs: x dimension of tiling pattern
	 * @param yimgs: y dimension of tiling pattern
	 * @param width: image width
	 * @param height: image height
	 * @param overlap: fractional overlap between images
	 * @param startxpos: starting x coordinate
	 * @param startypos: starting y coordinate
	 * @param startxdir: x tiling direction (-1 or +1)
	 * @param startydir: y tiling direction (-1 or +1)
	 * @param psize: the pixel size
	 * @return a list of coordinates
	 */
	public static float[][] getTileCoordsRaster(int ximgs,int yimgs,int width,int height,float overlap,int startxpos,int startypos,int startxdir,int startydir,float psize){
		float[][] coords=getTileCoordsRaster(ximgs,yimgs,width,height,overlap,startxpos,startypos,startxdir,startydir);
		if(psize==1.0f) return coords;
		for(int i=0;i<coords.length;i++){
			for(int j=0;j<coords[i].length;j++){
				coords[i][j]*=psize;
			}
		}
		return coords;
	}
	
	/*****************
	 * returns a list of rectangles corresponding the global coordinate list
	 * @return the list of rectangles (x,y,width,height)
	 */
	public int[][] getRois(){
		int[][] rois=new int[xvals.length][];
		for(int i=0;i<xvals.length;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			rois[i]=new int[]{xstart,ystart,width,height};
		}
		return rois;
	}
	
	/***************
	 * stitches an image set according to the global coordinate list (no feathering)
	 * @param input: the image set
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(Object[] input){
		//assume that xvals and yvals are in pixel units
		return stitch_frame(input,false,0.0f);
	}
	
	/*********************
	 * stitches an image set according to the global coordinate list, feathered with the specified overlap
	 * @param input: the image set
	 * @param feathered
	 * @param foverlap: fractional overlap
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(Object[] input,boolean feathered,float foverlap){
		//assume that xvals and yvals are in pixel units
		Object stitched=null;
		this.feather=feathered;
		if(typeindex==0) stitched=new byte[newwidth*newheight];
		if(typeindex==1) stitched=new short[newwidth*newheight];
		if(typeindex==2) stitched=new float[newwidth*newheight];
		this.foverlap=foverlap;
		hoverlap=foverlap*width;
		voverlap=foverlap*height;
		for(int i=0;i<input.length;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			pasteSubImage(stitched,input[i],xstart,ystart);
			gui.showProgress(i,input.length);
		}
		return stitched;
	}
	
	/*********************
	 * stitches an image set according to the global coordinate list, feathered with the specified overlap
	 * @param input: the image set
	 * @param feathered
	 * @param hoverlap: fractional horizontal overlap
	 * @param voverlap: fractional vertical overlap
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(Object[] input,boolean feathered,float hoverlap,float voverlap){
		//assume that xvals and yvals are in pixel units
		Object stitched=null;
		this.feather=feathered;
		if(typeindex==0) stitched=new byte[newwidth*newheight];
		if(typeindex==1) stitched=new short[newwidth*newheight];
		if(typeindex==2) stitched=new float[newwidth*newheight];
		this.foverlap=hoverlap/width;
		this.hoverlap=hoverlap;
		this.voverlap=voverlap;
		for(int i=0;i<input.length;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			pasteSubImage(stitched,input[i],xstart,ystart);
			gui.showProgress(i,input.length);
		}
		return stitched;
	}
	
	/******************
	 * stitches an image set according to the global coordinate list, feathered with the average overlap
	 * @param input: the image set
	 * @param feathered
	 * @param pairs: list of pair indices
	 * @param pairStats: the pair statistics (x,y,phasecorr,crosscorr)
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(Object[] input,boolean feathered,int[][] pairs,float[][] pairStats){
		//haven't tested this one
		//assume that xvals and yvals are in pixel units
		Object stitched=null;
		this.feather=feathered;
		if(typeindex==0) stitched=new byte[newwidth*newheight];
		if(typeindex==1) stitched=new short[newwidth*newheight];
		if(typeindex==2) stitched=new float[newwidth*newheight];
		this.hoverlap=0.0f; //this is just for feathering, set it for individuals
		this.voverlap=0.0f;
		for(int i=0;i<input.length;i++){
			if(feathered){ //find the connections, then use their overlap values
				int[][] connections=findConnectedImages(pairs,i,input.length);
				int nhcons=0;
				int nvcons=0;
				for(int j=0;j<connections.length;j++){
					int[] orient=pairOrientation(pairStats[connections[j][1]],width,height);
					if(orient[0]==0){
						this.hoverlap+=(float)width-Math.abs(pairStats[connections[j][1]][1]);
						nhcons++;
					}
					if(orient[0]==1){
						this.voverlap+=(float)height-Math.abs(pairStats[connections[j][1]][2]);
						nvcons++;
					}
				}
				this.hoverlap/=(float)nhcons;
				this.voverlap/=(float)nvcons;
			}
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			pasteSubImage(stitched,input[i],xstart,ystart);
			gui.showProgress(i,input.length);
		}
		return stitched;
	}
	
	/*********************
	 * stitches an image set according to the global coordinate list, feathered with the specified overlap
	 * @param input: a frame interface dealing the image set
	 * @param nframes
	 * @param feathered
	 * @param foverlap: fractional overlap
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(FrameInterface input,int nframes,boolean feathered,float foverlap){
		//assume that xvals and yvals are in pixel units
		Object stitched=null;
		this.feather=feathered;
		if(typeindex==0) stitched=new byte[newwidth*newheight];
		if(typeindex==1) stitched=new short[newwidth*newheight];
		if(typeindex==2) stitched=new float[newwidth*newheight];
		this.foverlap=foverlap;
		hoverlap=foverlap*width;
		voverlap=foverlap*height;
		for(int i=0;i<nframes;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			pasteSubImage(stitched,input.getNextFrame(),xstart,ystart);
		}
		return stitched;
	}
	
	/**********************
	 * stitches an image set according to the global coordinate list, feathered with the specified overlap
	 * @param input: a frame interface dealing the image set
	 * @param nframes
	 * @param feathered
	 * @param hoverlap: fractional horizontal overlap
	 * @param voverlap: fractional vertical overlap
	 * @return a stitched image of newwidth by newheight size
	 */
	public Object stitch_frame(FrameInterface input,int nframes,boolean feathered,float hoverlap,float voverlap){
		//assume that xvals and yvals are in pixel units
		Object stitched=null;
		this.feather=feathered;
		if(typeindex==0) stitched=new byte[newwidth*newheight];
		if(typeindex==1) stitched=new short[newwidth*newheight];
		if(typeindex==2) stitched=new float[newwidth*newheight];
		this.foverlap=hoverlap/(float)width;
		this.hoverlap=hoverlap;
		this.voverlap=voverlap;
		for(int i=0;i<nframes;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			pasteSubImage(stitched,input.getNextFrame(),xstart,ystart);
		}
		return stitched;
	}
	
	/*****************
	 * pastes a sub image into a bigger image with the global settings
	 * @param image: the final image
	 * @param subarray: the subimage
	 * @param xoff: horizontal offset
	 * @param yoff: vertical offset
	 */
	public void pasteSubImage(Object image,Object subarray,int xoff,int yoff){
		if(typeindex==0) pasteSubImage((byte[])image,(byte[])subarray,xoff,yoff);
		if(typeindex==1) pasteSubImage((short[])image,(short[])subarray,xoff,yoff);
		if(typeindex==2) pasteSubImage((float[])image,(float[])subarray,xoff,yoff);
	}
	
	/***************
	 * pastes a sub image into a bigger image with the global settings
	 * @param image: the final image
	 * @param subarray: the subimage
	 * @param xoff: horizontal offset
	 * @param yoff: vertical offset
	 */
	public void pasteSubImage(float[] image,float[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(subarray[k+j*width],image[temp],getWeight(k,j));
						} else {
							image[temp]=subarray[k+j*width];
						}
					}
				}
			}
		}
	}
	
	/***************
	 * pastes a sub image into a bigger image with the global settings
	 * @param image: the final image
	 * @param subarray: the subimage
	 * @param xoff: horizontal offset
	 * @param yoff: vertical offset
	 */
	public void pasteSubImage(short[] image,short[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(subarray[k+j*width],image[temp],getWeight(k,j));
						} else {
							image[temp]=subarray[k+j*width];
						}
					}
				}
			}
		}
	}
	
	/***************
	 * pastes a sub image into a bigger image with the global settings
	 * @param image: the final image
	 * @param subarray: the subimage
	 * @param xoff: horizontal offset
	 * @param yoff: vertical offset
	 */
	public void pasteSubImage(byte[] image,byte[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(subarray[k+j*width],image[temp],getWeight(k,j));
						} else {
							image[temp]=subarray[k+j*width];
						}
					}
				}
			}
		}
	}
	
	public byte combine_values(byte val1,byte val2,float weight){
		float comb=(val1&0xff)*weight+(val2&0xff)*(1.0f-weight);
		int comb2=(int)comb;  if(comb2>255) comb2=255; if(comb2<0) comb2=0;
		return (byte)comb2;
	}
	
	public short combine_values(short val1,short val2,float weight){
		float comb=(val1&0xffff)*weight+(val2&0xffff)*(1.0f-weight);
		int comb2=(int)comb;  if(comb2>65535) comb2=65535; if(comb2<0) comb2=0;
		return (short)comb2;
	}
	
	public float combine_values(float val1,float val2,float weight){
		//the weight is the weight of the first value
		return val1*weight+val2*(1.0f-weight);
	}
	
	/********************
	 * gets the image combination weights based on overlap values and settings
	 * @param x: the x coordinate of the subimage
	 * @param y: the y coordinate of the subimage
	 * @return
	 */
	public float getWeight(int x,int y){
		if(!feather) return 1.0f;
		float hweight=1.0f;
		if((float)x<hoverlap){
			hweight=(float)x/hoverlap;
		} else if((float)x>((float)width-1.0f-hoverlap)) {
			hweight=(float)(width-1-x)/hoverlap;
		}
		float vweight=1.0f;
		if((float)y<voverlap){
			vweight=(float)y/voverlap;
		} else if((float)y>((float)height-1.0f-voverlap)) {
			vweight=(float)(height-1-y)/voverlap;
		}
		//return Math.min(hweight,vweight);
		if(hweight==1.0f && vweight==1.0f) return 1.0f;
		//if(hweight==1.0f) return vweight;
		//if(vweight==1.0f) return hweight;
		return Math.min(hweight,vweight);
	}
	
	/***************
	 * a pairwise phase correlation
	 * @param pix1
	 * @param pix2
	 * @return the correlation (fftdim x fftdim size)
	 */
	public float[] phaseCorrPad(Object pix1,Object pix2){
		if(ccclass==null){
			int largestdim=width; if(height>largestdim) largestdim=height;
			int[] bestindex=fftutils.get_best_index((int)((float)largestdim*1.2f),true);
			fftdim=bestindex[1];
			ccclass=new crosscorr2D(fftdim,fftdim,bestindex[0],bestindex[0]);
		}
		float[] pix11=algutils.convert_arr_float2(pix1);
		float[] pix21=algutils.convert_arr_float2(pix2);
		//start by padding with linear decay
		float[] pad1=padding.pad_xy_mirrored(pix11,width,height,fftdim,fftdim,true);
		float[] pad2=padding.pad_xy_mirrored(pix21,width,height,fftdim,fftdim,true);
		//now cross correlate
		float[][] cc=ccclass.phaseCorr(pad1,pad2,false);
		//zero the origin: it is not important
		cc[0][0]=0.0f;
		cc[0]=jsmooth.smooth2D(cc[0],fftdim,fftdim);
		return cc[0];
	}
	
	/****************
	 * does a multithreaded calculation of pair phase correlations
	 * @param pix1: the image array
	 * @param pairs: the pair indices
	 * @param nthreads: number of thread workers to use
	 * @return the array of correlations (fftdim x fftdim size)
	 */
	public float[][] phaseCorrPad(Object[] pix1,int[][] pairs,int nthreads){
		//this is the multithreaded version
		if(ccclass==null){
			int largestdim=width; if(height>largestdim) largestdim=height;
			int[] bestindex=fftutils.get_best_index((int)((float)largestdim*1.2f),true);
			fftdim=bestindex[1];
			ccclass=new crosscorr2D(fftdim,fftdim,bestindex[0],bestindex[0]);
		}
		float[][] real=new float[pix1.length][];
		float[][] im=new float[pix1.length][fftdim*fftdim];
		for(int i=0;i<pix1.length;i++){
			real[i]=padding.pad_xy_mirrored(algutils.convert_arr_float2(pix1[i]),width,height,fftdim,fftdim,true);
		}
		//now do the forward fft
		po4realfft2D fft=ccclass.fft;
		fftutils.dorealfft2D(real,im,false,nthreads,fft);
		float[][] ccr=new float[pairs.length][];
		float[][] cci=new float[pairs.length][];
		for(int i=0;i<pairs.length;i++){
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
			float[][] temp=ccclass.phaseCorr(real[pair2],im[pair2],real[pair1],im[pair1]);
			ccr[i]=temp[0];
			cci[i]=temp[1];
		}
		fftutils.dorealfft2D(ccr,cci,true,nthreads,fft);
		cci=null;
		//zero the origin: it is not important
		for(int i=0;i<pairs.length;i++){
			ccr[i][0]=0.0f;
			ccr[i]=jsmooth.smooth2D(ccr[i],fftdim,fftdim);
		}
		return ccr;
	}
	
	/***********************
	 * find the cross correlation peaks by max not mask
	 * @param cc: the cross correlation (fftdim size)
	 * @param npeaks: number of peaks to find
	 * @param pix1: the first image of the correlation (for crosscorr calc)
	 * @param pix2: the second image of the correlation (for crosscorr calc)
	 * @param guessshift: the guessed horizontal and vertical shifts
	 * @return an array of peaks, each with phasecorr-max,x,y,phasecorr-neighborsavg,crosscorr
	 */
	public float[][] getCorrPeaks(float[] cc,int npeaks,Object pix1,Object pix2,float[] guessshift){
		//need to find npeaks maxima
		float[][] maxima=new float[npeaks][5]; //has maxval,xval,yval,neighborsavg, cross-corr
		int minindex=0;
		float minval=-1000.0f;
		for(int i=0;i<fftdim;i++){
			for(int j=0;j<fftdim;j++){
				float value=cc[j+i*fftdim];
				if(value>minval){ //our value made the top list
    				float[] neighbors=(float[])algutils.getNeighborsWrapped(cc,j,i,fftdim,fftdim);
    				float max=getMax(neighbors);
    				float avg=getAvg(neighbors,value);
    				if(max<value){ //we have a local maximum, replace the minval
    					maxima[minindex]=new float[]{value,(float)j,(float)i,avg,0.0f};
    					minindex=0; minval=maxima[0][0];
    					for(int k=1;k<npeaks;k++){ //refind the minval and its index
    						if(maxima[k][0]<minval){minval=maxima[k][0]; minindex=k;}
    					}
    				}
				}
			}
		}
		//now need to distinguish "noise" peaks from actual peaks
		//actual peaks will smoothly fall off so difference between neighborhood avg and peak should not be large
		//sort the peaks by height minus neighbor difference (actually just the smoothed value)
		//if we have points beyond half the image size, shift them to negative
		//also need to make predicted negative ones negative
		float[] metrics=new float[npeaks];
		for(int i=0;i<npeaks;i++){
			//metrics[i]=maxima[i][0];
			if(maxima[i][1]>(float)width) maxima[i][1]-=(float)fftdim;
			if(maxima[i][2]>(float)height) maxima[i][2]-=(float)fftdim;
			if(guessshift[0]<-1.0f && maxima[i][1]>1.0f) maxima[i][1]-=(float)fftdim;
			if(guessshift[1]<-1.0f && maxima[i][2]>1.0f) maxima[i][2]-=(float)fftdim;
			maxima[i][4]=get2DCorr(pix1,pix2,(int)maxima[i][1],(int)maxima[i][2]);
			metrics[i]=maxima[i][0];
		}
		int[] order=jsort.get_javasort_order(metrics); //this is from low to high
		float[][] maxima2=new float[npeaks][];
		for(int i=0;i<npeaks;i++) maxima2[i]=maxima[order[npeaks-i-1]];
		return maxima2;
	}
	
	/****************
	 * calculates the correlation between image1 and image2 (not by fft) for a single shift
	 * @param image1
	 * @param image2
	 * @param xshift
	 * @param yshift
	 * @return the crosscorr value
	 */
	public float get2DCorr(Object image1,Object image2,int xshift,int yshift){
		//this calculates the pearson correlation for only the overlapped region for a set of tolerance x tolerance shifts
		//xshift and yshift should always be in range for the image, if not, send a warning
		if(xshift>(width-1) || yshift>(height-1) || xshift<(1-width) || yshift<(1-height)){
			gui.showMessage("warning: correlation shift outside of image");
			return 0.0f;
		}
		int dtype=algutils.get_array_type(image1);
		float corr=0.0f;
		float var1=0.0f;
		float var2=0.0f;
		float avg1=0.0f;
		float avg2=0.0f;
		int xstart1=xshift;
        int xstart2=0;
        int xpts=width-xshift;
        if(xstart1<0){
        	xstart1=0; xstart2=-xshift; xpts=width+xshift;
        }
        int ystart1=yshift;
        int ystart2=0;
        int ypts=height-yshift;
        if(ystart1<0){
        	ystart1=0; ystart2=-yshift; ypts=height+yshift;
        }
        for(int i=0;i<ypts;i++){
        	for(int j=0;j<xpts;j++){
        		float val1=algutils.getPixelVal(image1,j+xstart1,i+ystart1,width,height,dtype);
        		float val2=algutils.getPixelVal(image2,j+xstart2,i+ystart2,width,height,dtype);
        		avg1+=val1;
				avg2+=val2;
				corr+=val1*val2;
				var1+=val1*val1;
				var2+=val2*val2;
        	}
        }
        float factor=(float)(xpts*ypts);
        avg1/=factor;
		avg2/=factor;
		var1/=factor;
		var2/=factor;
		corr/=factor;
		var1-=avg1*avg1;
		var2-=avg2*avg2;
		if(var1>0.0f && var2>0.0f) corr=(corr-avg1*avg2)/(float)Math.sqrt(var1*var2);
		else corr=0.0f;
		return corr;
	}
	
	public float getMax(float[] arr){
		float max=arr[0];
		for(int i=1;i<arr.length;i++) if(arr[i]>max) max=arr[i];
		return max;
	}
	
	public float getAvg(float[] arr,float extra){
		float avg=extra;
		for(int i=0;i<arr.length;i++) avg+=arr[i];
		return avg/(float)(arr.length+1);
	}
	
	/************************
	 * this is a non-global optimization of the alignment based on correlation shifts and statistics
	 * @param pairStats: the correlation statistics for pairs (x,y,phasecorr,crosscorr)
	 * @param fixedindex: the index of the fixed image (often the upper left but make sure it has significant correlation)
	 * @param pairs: the pair indices
	 * @param nimages
	 * @param corrthresh: correlation threshold
	 * @return the coordinate list for each image
	 */
	public static float[][] preAlign(float[][] pairStats,int fixedindex,int[][] pairs,int nimages, float corrthresh){
		//start with the fixed pair and propagate backwards and forwards until all are aligned
		int nunaligned=nimages-1;
		float[][] outcoords=new float[nimages][];
		outcoords[fixedindex]=new float[]{0.0f,0.0f};
		int prevunaligned=nunaligned;
		while(nunaligned>0){
			for(int i=0;i<outcoords.length;i++){
				if(outcoords[i]!=null){
					//we found an aligned one, see if it has neighbors with valid phase correlation
					//if it does, add those to the aligned pool and the coordinate system
					int[][] targets=findConnectedImages(pairs,i,nimages);
					if(targets==null) continue;
					for(int j=0;j<targets.length;j++){
						if(outcoords[targets[j][0]]==null){
    						//ComparePair pair=pairs2.get(targets[j][1]);
    						float[] pair=pairStats[targets[j][1]];
    						int targetmember=pairs[targets[j][1]][1];
    						float[] targetshift=new float[]{pair[0],pair[1]};
    						if(targetmember==i){//the query is the second member of the pair
    							targetmember=pairs[targets[j][1]][0];
    							targetshift[0]*=-1.0f;
    							targetshift[1]*=-1.0f;
    						}
    						if(pair[3]>corrthresh){
    							float[] querycoords=outcoords[i];
    							float[] newcoords={querycoords[0]+targetshift[0],querycoords[1]+targetshift[1]};
    							outcoords[targetmember]=newcoords;
    							nunaligned--;
    						}
						}
					}
				}
			}
			if(nunaligned>0 && nunaligned==prevunaligned){
				//we have isolated images, allow them to be aligned by pairs with invalid phase correlation
				//use the first found partner which is already aligned
				//just to one isolated image per round to let the others potentially catch up if possible
				for(int i=0;i<outcoords.length;i++){
					if(outcoords[i]==null){
						int[][] targets=findConnectedImages(pairs,i,nimages);
						for(int j=0;j<targets.length;j++){
							if(outcoords[targets[j][0]]!=null){
								float[] pair=pairStats[targets[j][1]];
								int targetmember=pairs[targets[j][1]][1]; //this time the target is the aligned image
	    						float[] targetshift=new float[]{pair[0],pair[1]};
	    						targetshift[0]*=-1.0f; //this time we are aligning the query to the target
	    						targetshift[1]*=-1.0f;
	    						if(targetmember==i){//the query is the second member of the pair
	    							targetmember=pairs[targets[j][1]][0];
	    							targetshift[0]*=-1.0f;
	    							targetshift[1]*=-1.0f;
	    						}
	    						float[] targetcoords=outcoords[targetmember];
	    						float[] newcoords={targetcoords[0]+targetshift[0],targetcoords[1]+targetshift[1]};
	    						outcoords[i]=newcoords;
	    						nunaligned--;
	    						break;
							}
						}
					}
					if(nunaligned!=prevunaligned) break;
				}
			}
			prevunaligned=nunaligned;
		}
		return outcoords;
	}
	
	/*public void globOptIterJay(Vector<ComparePair> pairs2,float[][] coords,int fixedindex,int[][] pairs,int nimages,float ftoler){
    	//here we perform one iteration of my poor-man's global optimization
    	//each set of coordinates is shifted in the direction of the weighted average of the surrounding shifts
    	//don't allow shifts greater than ftoler*dimension_size
    	//ignore the fixed image for now, and shift everything by its values at the end
    	//start with the image with highest overall weighted distance
    	float maxcorr=0.0f;
    	for(int i=0;i<pairs.length;i++){
    		float corr=pairs2.get(i).getCrossCorrelation();
    		if(corr>maxcorr) corr=maxcorr;
    	}
    	float[][] stats=new float[nimages][3];
    	for(int i=0;i<nimages;i++){
    		int[][] connected=findConnectedImages(pairs,i,nimages);
    		for(int j=0;j<connected.length;j++){
    			ComparePair pair=pairs2.get(connected[j][1]);
    			float corr=pair.getCrossCorrelation();
    			float currxshift=coords[connected[j][0]][0]-coords[connected[]]
    		}
    	}
	}*/
	
	/************************
	 * find all of the images connected to the query image
	 * @param pairs: image pair listing
	 * @param query: index of query image
	 * @param nimages
	 * @return
	 */
	public static int[][] findConnectedImages(int[][] pairs,int query,int nimages){
		int[][] targets=new int[nimages][2];  //first value is other image index, second is pair index
		int ntargets=0;
		for(int i=0;i<pairs.length;i++){
			if(pairs[i][0]==query){
				targets[ntargets][0]=pairs[i][1]; targets[ntargets][1]=i; ntargets++;
			}
			if(pairs[i][1]==query){
				targets[ntargets][0]=pairs[i][0]; targets[ntargets][1]=i; ntargets++;
			}
		}
		if(ntargets==0) return null;
		int[][] targets2=new int[ntargets][2];
		for(int i=0;i<ntargets;i++) targets2[i]=targets[i];
		return targets2;
	}
	
	/*******************
	 * creates of list of relative pair orientations and signs
	 * @param pairStats: the pair statistics (phasecorr, x, y, crosscorr)
	 * @param width: image width
	 * @param height: image height
	 * @return a list of pair types (1:h, 2:v, 3:d) and directions
	 */
	public static int[] pairOrientation(float[] pairStats,int width,int height){
		float halfwidth=0.5f*(float)width; float halfheight=0.5f*(float)height;
		int pairtype=0; //0=horizontal,1=vertical,2=diagonal
		int pairsign=1;
		if(Math.abs(pairStats[1])<halfwidth && Math.abs(pairStats[2])>halfheight){
			pairtype=1;
			if(pairStats[2]<-halfheight) pairsign=-1;
		}
		if(Math.abs(pairStats[1])>halfwidth && Math.abs(pairStats[2])>halfheight){
			pairtype=2;
		}
		if(pairtype!=1 && pairStats[1]<-halfwidth) pairsign=-1;
		return new int[]{pairtype,pairsign};
	}
	
	/****************
	 * check for internal consistency among coordinates, filling in median values when the correlation isn't good enough
	 * @param pairStats: the list of pair statistics (x,y,phasecorr,crosscorr)
	 * @param pairs: pair indices
	 * @param ftoler: fraction of the dimension size to tolerate
	 * @param corrthresh: correlation threshold
	 * @return two objects: 1-the median shift values (hh,hv,vh,vv,dh,dv); 2-the pair types for each pair (1:h, 2:v, 3:d)
	 */
	public Object[] checkCoords(float[][] pairStats,int[][] pairs,float ftoler,float corrthresh){
		//this subroutine uses the median displacements for vertical and horizontal pairs
		//if values are more than ftoler*dimension different than the median, they are converted to the median
		//they should also be within ftoler of the guess values
		int halfwidth=width/2;
		int halfheight=height/2;
		float[][] vertvals=new float[2][pairs.length];
		int numvert=0;
		float[][] horvals=new float[2][pairs.length];
		int numhor=0;
		float[][] diagvals=new float[2][pairs.length];
		int numdiag=0;
		int[][] pairtype=new int[pairs.length][2];
		for(int i=0;i<pairs.length;i++){
			//ComparePair pair=pairs2.get(i);
			float[] pair=pairStats[i];
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
			float xdiff=xvals[pair2]-xvals[pair1];
			float ydiff=yvals[pair2]-yvals[pair1];
			float axdiff=Math.abs(xdiff);
			float aydiff=Math.abs(ydiff);
    			if(axdiff>halfwidth && aydiff<halfheight){ //horizontal pairs
    				pairtype[i][0]=0; pairtype[i][1]=(xdiff>0)?1:-1; //forward or backwards pair
    				float mult=(float)pairtype[i][1];
    				if(pair[3]>=corrthresh){
        				horvals[0][numhor]=mult*pair[0]; horvals[1][numhor]=mult*pair[1]; numhor++; 
    				} else{pair[2]=0.0f; pair[3]=0.0f;}
    			}
    			if(axdiff<halfwidth && aydiff>halfheight){ //vertical pairs
    				pairtype[i][0]=1; pairtype[i][1]=(ydiff>0)?1:-1; //up or down pair
    				float mult=(float)pairtype[i][1];
    				if(pair[3]>=corrthresh){
        				vertvals[0][numvert]=mult*pair[0]; vertvals[1][numvert]=mult*pair[1]; numvert++;
    				} else{pair[2]=0.0f; pair[3]=0.0f;}
    			}
    			if(axdiff>halfwidth && aydiff>halfheight){ //diagonal pairs
    				pairtype[i][0]=2; pairtype[i][1]=(xdiff>0)?1:-1; //down-right or up-left pairs
    				float mult=(float)pairtype[i][1];
    				if(pair[3]>=corrthresh){
        				diagvals[0][numdiag]=mult*pair[0]; diagvals[1][numdiag]=mult*pair[1]; numdiag++; 
    				} else{pair[2]=0.0f; pair[3]=0.0f;}
    			}
		}
		float htoler=ftoler*(float)width; //pixel tolerance for difference from median and guess
		float vtoler=ftoler*(float)height;
		float dtoler=(float)Math.sqrt(htoler*htoler+vtoler*vtoler);
		boolean hout=false,vout=false,dout=false;
		float vguessxdiff=0.0f, vguessydiff=0.0f, hguessxdiff=0.0f, hguessydiff=0.0f, dguessxdiff=0.0f, dguessydiff=0.0f;
		//start by getting the expected guess shift values
		for(int i=0;i<pairs.length;i++){
			if(pairtype[i][0]==0){
				float xdiff=xvals[pairs[i][1]]-xvals[pairs[i][0]];
				float ydiff=yvals[pairs[i][1]]-yvals[pairs[i][0]];
				hguessxdiff=(float)pairtype[i][1]*xdiff; hguessydiff=(float)pairtype[i][1]*ydiff;
				break;
			}
		}
		for(int i=0;i<pairs.length;i++){
			if(pairtype[i][0]==1){
				float xdiff=xvals[pairs[i][1]]-xvals[pairs[i][0]];
				float ydiff=yvals[pairs[i][1]]-yvals[pairs[i][0]];
				vguessxdiff=(float)pairtype[i][1]*xdiff; vguessydiff=(float)pairtype[i][1]*ydiff;
				break;
			}
		}
		for(int i=0;i<pairs.length;i++){
			if(pairtype[i][0]==2){
				float xdiff=xvals[pairs[i][1]]-xvals[pairs[i][0]];
				float ydiff=yvals[pairs[i][1]]-yvals[pairs[i][0]];
				dguessxdiff=(float)pairtype[i][1]*xdiff; dguessydiff=(float)pairtype[i][1]*ydiff;
				break;
			}
		}
		//get the median statistics for our data
		float hmedxdiff=0.0f, hmedydiff=0.0f, hdist=0.0f;
		if(numhor>1){
    		hmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(horvals[0],0,numhor),null);
    		hmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(horvals[1],0,numhor),null);
    		hdist=(float)Math.sqrt((hmedxdiff-hguessxdiff)*(hmedxdiff-hguessxdiff)+(hmedydiff-hguessydiff)*(hmedydiff-hguessydiff));

		}
		if(hdist>htoler || numhor<2){
			hmedxdiff=hguessxdiff; hmedydiff=hguessydiff;
			hout=true;
			gui.showMessage("horizontal median beyond tolerance, setting to guess");
		}
		gui.showMessage("horizontal median shift");
		gui.showMessage(""+hmedxdiff+" , "+hmedydiff);
		float vmedxdiff=0.0f, vmedydiff=0.0f, vdist=0.0f;;
		if(numvert>1){
			vmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(vertvals[0],0,numvert),null);
			vmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(vertvals[1],0,numvert),null);
			vdist=(float)Math.sqrt((vmedxdiff-vguessxdiff)*(vmedxdiff-vguessxdiff)+(vmedydiff-vguessydiff)*(vmedydiff-vguessydiff));
			//gui.showMessage(""+vguessydiff+" , "+vdist);

		}
		if(vdist>vtoler || numvert<2){
			vmedxdiff=vguessxdiff; vmedydiff=vguessydiff;
			vout=true;
			gui.showMessage("vertical median beyond tolerance, setting to guess");
		}
		gui.showMessage("vertical median shift");
		gui.showMessage(""+vmedxdiff+" , "+vmedydiff);
		float dmedxdiff=0.0f; float dmedydiff=0.0f, ddist=0.0f;;
		if(numdiag>1){ //we might not have diagonal pairs
			dmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(diagvals[0],0,numdiag),null);
			dmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(diagvals[1],0,numdiag),null);
			ddist=(float)Math.sqrt((dmedxdiff-dguessxdiff)*(dmedxdiff-dguessxdiff)+(dmedydiff-dguessydiff)*(dmedydiff-dguessydiff));
		}
		if(ddist>dtoler || numdiag<2){
			dmedxdiff=dguessxdiff; dmedydiff=dguessydiff;
			dout=true;
			gui.showMessage("diagonal median beyond tolerance, setting to guess");
		}
		gui.showMessage("diagonal median shift");
		gui.showMessage(""+dmedxdiff+" , "+dmedydiff);

		int hcounter=0; int vcounter=0; int dcounter=0;
		for(int i=0;i<pairs.length;i++){
			//ComparePair pair=pairs2.get(i);
			float[] pair=pairStats[i];
			if(pair[3]>=corrthresh){
    			if(pairtype[i][0]==0){ //horizontal
    				if(hout || Math.abs(horvals[0][hcounter]-hmedxdiff)>htoler){
    					//the pair is out of bounds, set its offset to the median and zero the correlation
    					float mult=(float)pairtype[i][1];
    					pair[0]=mult*hmedxdiff; pair[1]=mult*hmedydiff;
    					pair[2]=0.0f; pair[3]=0.0f;
    				}
    				hcounter++; //counts the valid shifts
    			}
    			if(pairtype[i][0]==1){ //vertical
    				if(vout || Math.abs(vertvals[1][vcounter]-vmedydiff)>vtoler){
    					float mult=(float)pairtype[i][1];
    					pair[0]=mult*vmedxdiff; pair[1]=mult*vmedydiff;
    					pair[2]=0.0f; pair[3]=0.0f;
    				}
    				vcounter++;
    			}
    			if(pairtype[i][0]==2){ //diagonal
    				if(dout || Math.abs(diagvals[0][dcounter]-dmedxdiff)>htoler){
    					float mult=(float)pairtype[i][1];
    					pair[0]=mult*dmedxdiff; pair[1]=mult*dmedydiff;
    					pair[2]=0.0f; pair[3]=0.0f;
    				}
    				dcounter++;
    			}
			} else {
				if(pairtype[i][0]==0){pair[0]=(float)pairtype[i][1]*hmedxdiff; pair[1]=(float)pairtype[i][1]*hmedydiff;}
				if(pairtype[i][0]==1){pair[0]=(float)pairtype[i][1]*vmedxdiff; pair[1]=(float)pairtype[i][1]*vmedydiff;}
				if(pairtype[i][0]==2){pair[0]=(float)pairtype[i][1]*dmedxdiff; pair[1]=(float)pairtype[i][1]*dmedydiff;}
			}
		}
		float[] medvals={hmedxdiff,hmedydiff,vmedxdiff,vmedydiff,dmedxdiff,dmedydiff};
		return new Object[]{medvals,pairtype};
	}
	
	/*********************
	 * returns a projection (max for now) of phase correlations for each indicated image pair along its overlapped seam
	 * @param pix1
	 * @param pairs
	 * @param guessshifts
	 * @param fftdim1: subimage size to use for the fft--assume power of 2
	 * @return
	 */
	public float[][] phaseCorrSTICS(Object[] pix1,int[][] pairs,float[][] guessshifts,int fftdim1,int nthreads,String projstat){
		int[] bestindex=fftutils.get_best_index(fftdim1,false,19);
		crosscorr2D ccclass1=new crosscorr2D(fftdim1,fftdim1,bestindex[0],bestindex[0]);
		float[][] ccproj=new float[pairs.length][fftdim1*fftdim1];
		for(int i=0;i<pairs.length;i++){
			int[] orientation=pairOrientation(new float[]{0.0f,guessshifts[i][0],guessshifts[i][1]},width,height);
			float[][] seamsubpairs=getSeamSubImagePairs(pix1[pairs[i][0]],pix1[pairs[i][1]],orientation,guessshifts[i],fftdim1);
			float[][] im=new float[seamsubpairs.length][fftdim1*fftdim1];
			po4realfft2D fft=ccclass1.fft;
			fftutils.dorealfft2D(seamsubpairs,im,false,nthreads,fft);
			float[][] ccr=new float[seamsubpairs.length/2][];
			float[][] cci=new float[seamsubpairs.length/2][];
			for(int j=0;j<seamsubpairs.length/2;j++){
				float[][] temp=ccclass1.phaseCorr(seamsubpairs[2*j],im[2*j],seamsubpairs[2*j+1],im[2*j+1]);
				ccr[j]=temp[0];
				cci[j]=temp[1];
			}
			fftutils.dorealfft2D(ccr,cci,true,nthreads,fft);
			cci=null;
			//zero the origin: it is not important
			for(int j=0;j<seamsubpairs.length/2;j++){
				ccr[j][0]=0.0f;
				ccr[j]=jsmooth.smooth2D(ccr[j],fftdim1,fftdim1);
			}
			ccproj[i]=algutils.get_stack_proj_stat(projstat,ccr,fftdim1,fftdim1,ccr.length,null);
			/*for(int j=0;j<seamsubpairs.length;j++){
				float[] corr=ccclass1.phaseCorr(algutils.convert_arr_float2(seamsubpairs[j][0]),algutils.convert_arr_float2(seamsubpairs[j][1]),false)[0];
				corr[0]=0.0f;
				jsmooth.smooth2D(corr,fftdim1,fftdim1);
				for(int k=0;k<fftdim1*fftdim1;k++){
					if(ccproj[i][k]<corr[k]) ccproj[i][k]=corr[k];
				}
			}*/
		}
		return ccproj;
	}
	
	/*********************
	 * this gets subimages a long the common seam of two images, converts the subimages to floats
	 * @param image1
	 * @param image2
	 * @param orientation: from pairOrientation
	 * @param guessshifts
	 * @param subsize: the subsize of the seam images
	 * @return an array of image pairs
	 */
	public float[][] getSeamSubImagePairs(Object image1,Object image2,int[] orientation,float[] guessshift,int subsize){
		//should probably evenly space the images but for now just shift the last one so it doesn't overextend
		if(orientation[0]==0){ //horizontal pairs
			Object image11=image1; Object image12=image2;
			float[] shifts=guessshift.clone();
			if(orientation[1]<0){ //if the second image is to the left, swap them
				image11=image2; image12=image1;
				shifts[0]=-guessshift[0]; shifts[1]=-guessshift[1];
			}
			int xstart1=(int)shifts[0]; int xstart2=0;
			int ystart1=(int)shifts[1]; if(shifts[1]<0.0f) ystart1=0;
			int yend1=height-1; if(shifts[1]<0.0f) yend1+=(int)shifts[1];
			float fnoverlaps=(float)(yend1-ystart1)/(float)subsize;
			int noverlaps=(int)fnoverlaps;
			if(fnoverlaps>(float)noverlaps) noverlaps++;
			float[][] overlaps=new float[noverlaps*2][];
			for(int i=0;i<noverlaps;i++){
				int y1=ystart1+i*subsize; if((y1+subsize)>yend1) y1=yend1-subsize;
				int y2=y1+(int)shifts[1];
				overlaps[2*i]=algutils.convert_arr_float2(padding.get_region2_padmirrored(image11,xstart1,y1,subsize,subsize,width,height));
				overlaps[2*i+1]=algutils.convert_arr_float2(padding.get_region2_padmirrored(image12,xstart2,y2,subsize,subsize,width,height));
			}
			return overlaps;
		} else if(orientation[0]==1){ //vertical pairs
			Object image11=image1; Object image12=image2;
			float[] shifts=guessshift.clone();
			if(orientation[1]<0){ //if the second image is above, swap them
				image11=image2; image12=image1;
				shifts[0]=-guessshift[0]; shifts[1]=-guessshift[1];
			}
			int xstart1=(int)shifts[0]; if(shifts[0]<0.0f) xstart1=0;
			int xend1=width-1; if(shifts[0]<0.0f) xend1+=(int)shifts[0];
			int ystart1=(int)shifts[1]; int ystart2=0;
			float fnoverlaps=(float)(xend1-xstart1)/(float)subsize;
			int noverlaps=(int)fnoverlaps;
			if(fnoverlaps>(float)noverlaps) noverlaps++;
			float[][] overlaps=new float[noverlaps*2][];
			for(int i=0;i<noverlaps;i++){
				int x1=xstart1+i*subsize; if((x1+subsize)>xend1) x1=xend1-subsize;
				int x2=x1+(int)shifts[0];
				overlaps[2*i]=algutils.convert_arr_float2(padding.get_region2_padmirrored(image11,x1,ystart1,subsize,subsize,width,height));
				overlaps[2*i+1]=algutils.convert_arr_float2(padding.get_region2_padmirrored(image12,x2,ystart2,subsize,subsize,width,height));
			}
			return overlaps;
		} else {
			//for now, don't consider diagonal pairs
			return null;
		}
	}
	
	/******************
	 * here we get the peaks from a STICS projection--need to back calculate to the original image coordinates
	 * @param cc
	 * @param npeaks
	 * @param pix1
	 * @param pix2
	 * @param guessshift
	 * @param fftdim1
	 * @return
	 */
	public float[][] getCorrSTICSPeaks(float[] cc,int npeaks,Object pix1,Object pix2,float[] guessshift,int fftdim1){
		//need to find npeaks maxima
		float[][] maxima=new float[npeaks][5]; //has maxval,xval,yval,neighborsavg, cross-corr
		int minindex=0;
		float minval=-1000.0f;
		for(int i=0;i<fftdim1;i++){
			for(int j=0;j<fftdim1;j++){
				float value=cc[j+i*fftdim1];
				if(value>minval){ //our value made the top list
    				float[] neighbors=(float[])algutils.getNeighborsWrapped(cc,j,i,fftdim1,fftdim1);
    				float max=getMax(neighbors);
    				float avg=getAvg(neighbors,value);
    				if(max<value){ //we have a local maximum, replace the minval
    					maxima[minindex]=new float[]{value,(float)j,(float)i,avg,0.0f};
    					minindex=0; minval=maxima[0][0];
    					for(int k=1;k<npeaks;k++){ //refind the minval and its index
    						if(maxima[k][0]<minval){minval=maxima[k][0]; minindex=k;}
    					}
    				}
				}
			}
		}
		//our stics images should only be shifted by the difference between the guessed shifts and the real shifts
		//if the orientation direction is negative, the images were swapped
		int[] orientation=pairOrientation(new float[]{0.0f,guessshift[0],guessshift[1]},width,height);
		Object image11=pix1; Object image12=pix2;
		float[] shifts=guessshift.clone();
		if(orientation[1]<0){ //if the second image is to the left, swap them
			image11=pix2; image12=pix1;
			shifts[0]=-guessshift[0]; shifts[1]=-guessshift[1];
		}
		float[] metrics=new float[npeaks];
		for(int i=0;i<npeaks;i++){
			//start by getting the negative shifts (resolve wrapping)
			if(maxima[i][1]>fftdim1/2) maxima[i][1]=maxima[i][1]-fftdim1;
			if(maxima[i][2]>fftdim1/2) maxima[i][2]=maxima[i][2]-fftdim1;
			//now add to the guessed shifts
			maxima[i][1]=maxima[i][1]+(int)shifts[0];
			maxima[i][2]=maxima[i][2]+(int)shifts[1];
			maxima[i][4]=get2DCorr(image11,image12,(int)maxima[i][1],(int)maxima[i][2]);
			metrics[i]=maxima[i][0];
		}
		int[] order=jsort.get_javasort_order(metrics); //this is from low to high
		float[][] maxima2=new float[npeaks][];
		for(int i=0;i<npeaks;i++) maxima2[i]=maxima[order[npeaks-i-1]];
		//now go back and make the shifts negative if the oriention was negative
		if(orientation[1]<0){
    		for(int i=0;i<npeaks;i++){
    			maxima2[i][1]=-maxima2[i][1];
    			maxima2[i][2]=-maxima2[i][2];
    		}
		}
		return maxima2;
	}
	
}
