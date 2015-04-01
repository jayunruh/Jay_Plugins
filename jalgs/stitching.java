package jalgs;

import jalgs.jfft.crosscorr2D;
import jalgs.jfft.fftutils;

public class stitching{
	public int width,height,typeindex,newwidth,newheight;
	public float xmin,xmax,ymin,ymax,foverlap,hoverlap,voverlap;
	public float[] xvals,yvals;
	public boolean feather;
	
	public stitching(int width,int height,float[] xvals,float[] yvals,int typeindex){
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
	}
	
	public stitching(int width,int height,float[] xvals1,float[] yvals1,float psize,int typeindex){
		//here we use calibrated x and y values
		this.width=width;
		this.height=height;
		this.typeindex=typeindex;
		xvals=new float[xvals1.length];
		yvals=new float[yvals1.length];
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
	}
	
	public static float[][] getTileCoords(int ximgs,int yimgs,int width,int height,float overlap){
		//this version implements snake by rows starting at bottom left
		return getTileCoords(ximgs,yimgs,width,height,overlap,0,yimgs-1,1,-1);
	}
	
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
	
	public int[][] getRois(){
		int[][] rois=new int[xvals.length][];
		for(int i=0;i<xvals.length;i++){
			int xstart=(int)(xvals[i]-xmin);
			int ystart=(int)(yvals[i]-ymin);
			rois[i]=new int[]{xstart,ystart,width,height};
		}
		return rois;
	}
	
	public Object stitch_frame(Object[] input){
		//assume that xvals and yvals are in pixel units
		return stitch_frame(input,false,0.0f);
	}
	
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
		}
		return stitched;
	}
	
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
	
	public void pasteSubImage(Object image,Object subarray,int xoff,int yoff){
		if(typeindex==0) pasteSubImage((byte[])image,(byte[])subarray,xoff,yoff);
		if(typeindex==1) pasteSubImage((short[])image,(short[])subarray,xoff,yoff);
		if(typeindex==2) pasteSubImage((float[])image,(float[])subarray,xoff,yoff);
	}
	
	public void pasteSubImage(float[] image,float[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(image[temp],subarray[k+j*width],getWeight(k,j));
						} else {
							image[temp]=subarray[k+j*width];
						}
					}
				}
			}
		}
	}
	
	public void pasteSubImage(short[] image,short[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(image[temp],subarray[k+j*width],getWeight(k,j));
						} else {
							image[temp]=subarray[k+j*width];
						}
					}
				}
			}
		}
	}
	
	public void pasteSubImage(byte[] image,byte[] subarray,int xoff,int yoff){
		for(int j=0;j<height;j++){
			int ypos=yoff+j;
			if(ypos>=0 && ypos<newheight){
				for(int k=0;k<width;k++){
					int xpos=xoff+k;
					if(xpos>=0 && xpos<newwidth){
						int temp=xpos+ypos*newwidth;
						if(image[temp]!=0.0f){
							image[temp]=combine_values(image[temp],subarray[k+j*width],getWeight(k,j));
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
		return val1*weight+val2*(1.0f-weight);
	}
	
	public float getWeight(int x,int y){
		if(!feather) return 1.0f;
		float hweight=1.0f;
		if(x<hoverlap){
			hweight=(hoverlap-x)/hoverlap;
		} else if(x>(width-1-hoverlap)) {
			hweight=(width-1-x)/hoverlap;
		}
		float vweight=1.0f;
		if(y<voverlap){
			vweight=(voverlap-y)/voverlap;
		} else if(y>(height-1-voverlap)) {
			vweight=(height-1-y)/voverlap;
		}
		return Math.min(hweight,vweight);
	}
	
	public float[] get2DOffsets(Object image1,Object image2,int xoff,int yoff,int tolerance){
		//here we calculate a tolerance x tolerance sized correlation image between image1 and image2 and find its maximum
		int dtype=algutils.get_array_type(image1);
		int halftoler=(int)(0.5f*tolerance);
		float[] corr=get2DCorr(image1,image2,xoff,yoff,tolerance);
		float maxcorr=corr[0];
		int maxx=0;
		int maxy=0;
		int counter=0;
		for(int i=0;i<tolerance;i++){
			for(int j=0;j<tolerance;j++){
				if(corr[counter]>maxcorr){
					maxcorr=corr[counter];
					maxx=j;
					maxy=i;
				}
				counter++;
			}
		}
		return new float[]{maxx-halftoler,maxy-halftoler};
	}
	
	public float[] get2DCorr(Object image1,Object image2,int xoff,int yoff,int tolerance){
		int dtype=algutils.get_array_type(image1);
		float[] corr=new float[tolerance*tolerance];
		int halftoler=(int)(0.5f*tolerance);
		for(int i=0;i<tolerance;i++){
			int yshift=yoff-halftoler+i;
			int overy=height-yshift;
			for(int j=0;j<tolerance;j++){
				int xshift=xoff-halftoler+j;
				int overx=width-xshift;
				int temp=j+i*tolerance;
				for(int k=yshift;k<height;k++){
					for(int l=xshift;l<width;l++){
						corr[temp]+=algutils.getPixelVal(image1,l,k,width,height,dtype)*algutils.getPixelVal(image2,l-xshift,k-yshift,width,height,dtype);
					}
				}
				corr[temp]/=overy*overx;
			}
		}
		return corr;
	}
	
	public float[] get2DCorrFFT(Object image1,Object image2,int xoff,int yoff,int tolerance){
		int minsize=Math.min(width,height);
		int[] bestindex=fftutils.get_best_index(minsize,false,19);
		//int po2width=fftutils.trim_length(width,false,4);
		//int po2height=fftutils.trim_length(height,false,4);
		//int po2size=po2width; if(po2height<po2size) po2size=po2height;
		int po2size=bestindex[1];
		float[][] cropped=bestcrop(image1,image2,xoff,yoff,po2size);
		crosscorr2D cc=new crosscorr2D(po2size,po2size,bestindex[0],bestindex[0]);
		float[] corr=cc.docrosscorr2D(cropped[0],cropped[1],false,false,false,false);
		//need to shift the correlation to the new xoff,yoff
		//manipulate_quads.shiftxycenter(corr,po2size/2-newxoff,po2size/2-newyoff,po2size,po2size);
		return algutils.get_region(corr,po2size/2,po2size/2,tolerance,tolerance,po2size,po2size);
	}
	
	public float[][] bestcrop(Object image1,Object image2,int xoff,int yoff,int po2size){
		//need to return the new xoff and yoff
		if(yoff==0){
			if(xoff>0){
				//take the right of image1 and the left of image2
				int startx=width-po2size-1;
				int starty=height/2-po2size/2;
				float[] region1=algutils.get_region2(image1,startx,starty,po2size,po2size,width,height);
				float[] region2=algutils.get_region2(image2,0,starty,po2size,po2size,width,height);
				return new float[][]{region1,region2};
			} else {
				//take the left of image1 and the right of image2
				int startx=width-po2size-1;
				int starty=height/2-po2size/2;
				float[] region1=algutils.get_region2(image1,0,starty,po2size,po2size,width,height);
				float[] region2=algutils.get_region2(image2,startx,starty,po2size,po2size,width,height);
				return new float[][]{region1,region2};
			}
		}
		if(xoff==0){
			if(yoff>0){
				//take the bottom of image1 and the top of image2
				int startx=width/2-po2size/2;
				int starty=height-po2size-1;
				float[] region1=algutils.get_region2(image1,startx,starty,po2size,po2size,width,height);
				float[] region2=algutils.get_region2(image2,startx,0,po2size,po2size,width,height);
				return new float[][]{region1,region2};
			} else {
				//take the top of image1 and the bottom of image2
				int startx=width/2-po2size/2;
				int starty=height-po2size-1;
				float[] region1=algutils.get_region2(image1,startx,0,po2size,po2size,width,height);
				float[] region2=algutils.get_region2(image2,startx,starty,po2size,po2size,width,height);
				return new float[][]{region1,region2};
			}
		}
		//here the offsets are diagonal
		if(xoff>0 && yoff>0){
			//take the bottom right of image 1 and top left of image2
			int startx=width-po2size-1;
			int starty=height-po2size-1;
			float[] region1=algutils.get_region2(image1,startx,starty,po2size,po2size,width,height);
			float[] region2=algutils.get_region2(image2,0,0,po2size,po2size,width,height);
			return new float[][]{region1,region2};
		}
		if(xoff<0 && yoff>0){
			//take the bottom left of image 1 and top right of image2
			int startx=width-po2size-1;
			int starty=height-po2size-1;
			float[] region1=algutils.get_region2(image1,0,starty,po2size,po2size,width,height);
			float[] region2=algutils.get_region2(image2,startx,0,po2size,po2size,width,height);
			return new float[][]{region1,region2};
		}
		if(xoff<0 && yoff>0){
			//take the bottom left of image 1 and top right of image2
			int startx=width-po2size-1;
			int starty=height-po2size-1;
			float[] region1=algutils.get_region2(image1,0,starty,po2size,po2size,width,height);
			float[] region2=algutils.get_region2(image2,startx,0,po2size,po2size,width,height);
			return new float[][]{region1,region2};
		}
		if(xoff>0 && yoff<0){
			//take the top right of image 1 and the bottom left of image2
			int startx=width-po2size-1;
			int starty=height-po2size-1;
			float[] region1=algutils.get_region2(image1,startx,0,po2size,po2size,width,height);
			float[] region2=algutils.get_region2(image2,0,starty,po2size,po2size,width,height);
			return new float[][]{region1,region2};
		}
		if(xoff<0 && yoff<0){
			//take the top right of image 1 and the bottom left of image2
			int startx=width-po2size-1;
			int starty=height-po2size-1;
			float[] region1=algutils.get_region2(image1,0,0,po2size,po2size,width,height);
			float[] region2=algutils.get_region2(image2,startx,starty,po2size,po2size,width,height);
			return new float[][]{region1,region2};
		}
		return null;
	}
	
}
