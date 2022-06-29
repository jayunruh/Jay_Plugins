import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import jalgs.*;
import jguis.*;

public class jay_stitching_plugin_jru_v1 implements PlugIn, gui_interface {

	public void run(String arg) {
		//ImagePlus imp = IJ.getImage();
		//int width=imp.getWidth(); int height=imp.getHeight();
		//start by getting the image and the position plot (if we have it)
		ImagePlus[] imps=jutils.selectImages(true,2,new String[]{"Mosaic_Stack","Position_Plot (optional)"});
		if(imps==null) return;
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		ImageStack stack=imps[0].getStack();
		int slices=stack.getSize();
		//get the global options for the stitching
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Find Shifts?",true);
		gd.addNumericField("Median_Shift_Tolerance (%)",10,5,15,null);
		gd.addNumericField("Correlation_Threshold",0.3,5,15,null);
		gd.addNumericField("Number of threads",4,0);
		gd.addCheckbox("Output Correlation",true);
		String[] selectoptions={"Highest_Crosscorr","Highest_Phasecorr","Closest_to_guess"};
		gd.addChoice("Peak Selection",selectoptions,selectoptions[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean findshifts=gd.getNextBoolean();
		float medshifttoler=(float)gd.getNextNumber()/100.0f;
		float corrthresh=(float)gd.getNextNumber();
		int nthreads=(int)gd.getNextNumber();
		boolean outcorr=gd.getNextBoolean();
		int peaksel=gd.getNextChoiceIndex();
		
		//now initialize the coordinates if we didn't import them
		PlotWindow4 coordplot=null;
		//coords is an 2 x ntiles array containing the coordinates in pixel units
		float[][] coords=null;
		//pairs is a ntiles x 2 array containing the pair identities
		int[][] pairs=null;
		if(imps[1]==null){
			//get all of the parameters to initialize the coordinates
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addNumericField("Overlap Percent",20.0,5,15,null);
			gd2.addNumericField("X Images",3,0);
			gd2.addNumericField("Y Images",3,0);
			String[] scantype={"Snake_by_rows","Raster_by_rows"};
			gd2.addChoice("Scan_Type",scantype,scantype[0]);
			String[] startoptions={"Upper_Left","Upper_Right","Lower_Left","Lower_Right"};
			gd2.addChoice("Scan_Start",startoptions,startoptions[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			float overlap=(float)gd2.getNextNumber()/100.0f;
			int ximgs=(int)gd2.getNextNumber();
			int yimgs=(int)gd2.getNextNumber();
			int scantypeoption=gd2.getNextChoiceIndex();
			int scandiroption=gd2.getNextChoiceIndex();
			int xstart=0; int ystart=yimgs-1; int xscandir=1; int yscandir=-1;
			if(scandiroption==0){ystart=0; yscandir=1;}
			if(scandiroption==1){xstart=ximgs-1; xscandir=-1; ystart=0; yscandir=1;}
			if(scandiroption==3){xstart=ximgs-1; xscandir=-1;}
			if(scantypeoption==0){
				coords=stitching.getTileCoords(ximgs,yimgs,width,height,overlap,xstart,ystart,xscandir,yscandir,1.0f);
				pairs=stitching.getPairs(ximgs,yimgs,xscandir,yscandir,false);
			} else{
				coords=stitching.getTileCoordsRaster(ximgs,yimgs,width,height,overlap,xstart,ystart,xscandir,yscandir,1.0f);
				pairs=stitching.getPairsRaster(ximgs,yimgs,xscandir,yscandir,false);
			}
			coordplot=new PlotWindow4("Stitching_Coordinates","x","y",coords[0],coords[1]);
			coordplot.draw();
		} else {
			float psize=(float)jutils.get_psize(imps[0]);
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addCheckbox("Plot_Has_Units",false);
			gd2.addNumericField("Pixel_Size",psize,5,15,null);
			gd2.showDialog(); if(gd2.wasCanceled()) return;
			boolean hasunits=gd2.getNextBoolean();
			psize=(float)gd2.getNextNumber();
			coordplot=jutils.getPW4SelCopy(imps[1].getWindow());
			//coordplot.setTitle("Coordinates");
			coordplot.getImagePlus().setTitle("Stitching_Coordinates");
			coords=new float[2][];
			coords[0]=coordplot.getXValues()[0];
			coords[1]=coordplot.getYValues()[0];
			float[][] tempcoords=new float[coords[0].length][2];
			for(int i=0;i<coords[0].length;i++){
				tempcoords[i][0]=coords[0][i];
				tempcoords[i][1]=coords[1][i];
			}
			if(hasunits){
				float minx=coords[0][0]/psize; float miny=coords[1][0]/psize;
				for(int i=0;i<coords[0].length;i++){
					coords[0][i]/=psize; coords[1][i]/=psize;
					if(coords[0][i]<minx) minx=coords[0][i];
					if(coords[1][i]<miny) miny=coords[1][i];
				}
				for(int i=0;i<coords[0].length;i++){
					coords[0][i]-=minx; coords[1][i]-=miny;
					tempcoords[i][0]=coords[0][i]; tempcoords[i][1]=coords[1][i];
				}
				coordplot.updatePlot();
				coordplot.autoscale();
			}
			float[][] temp=stitching.getPairs(tempcoords,width,height,false);
			pairs=new int[temp.length][2];
			for(int i=0;i<temp.length;i++){pairs[i][0]=(int)temp[i][0]; pairs[i][1]=(int)temp[i][1];}
			IJ.log(""+pairs.length+" pairs found");
		}

		stitching s=new stitching(width,height,coords[0],coords[1],algutils.get_array_type(stack.getPixels(1)),this);

		//Object[] pix=jutils.stack2array(stack);
		int chans=imps[0].getNChannels();
		int slices2=imps[0].getNSlices();
		int frames=imps[0].getNFrames();
		int currchan=imps[0].getC();
		int currslice=imps[0].getZ();
		Object[] pix=jutils.get3DTSeries(stack,currslice-1,currchan-1,frames,slices2,chans);
		if(pix.length<2) pix=jutils.get3DZSeries(stack,currchan-1,0,frames,slices2,chans);
		Object stitched=null;
		if(!findshifts){
			float[] overlaps=getAvgOverlap(new float[][]{coords[0],coords[1]},width,height);
			float hoverlap=overlaps[0]; float voverlap=overlaps[1];
			IJ.log("overlaps: "+hoverlap+" , "+voverlap);
			stitched=s.stitch_frame(pix,true,hoverlap,voverlap);
		} else {
			IJ.showStatus("calculating correlation");
			ImageStack ccstack=null;
			float[][] stats=new float[pairs.length][8];
			float[][] pairStats=new float[pairs.length][4];
			float[][] cc=s.phaseCorrPad(pix,pairs,nthreads); //this is multithreaded
			IJ.showStatus("finding peaks");
			IJ.log("peak statistics: corr,phasecorr,x,y");
			for(int i=0;i<pairs.length;i++){
				int pair1=pairs[i][0]; int pair2=pairs[i][1];
				IJ.log(""+pair1+ " <- "+pair2+":");
				float predxshift=coords[0][pair2]-coords[0][pair1];
    				float predyshift=coords[1][pair2]-coords[1][pair1];
				stats[i][0]=(float)pair1; stats[i][1]=(float)pair2;
				stats[i][2]=predxshift; stats[i][3]=predyshift;
				if(ccstack==null) ccstack=new ImageStack(s.fftdim,s.fftdim);
				ccstack.addSlice(""+pair1+","+pair2,cc[i]);
				float[][] peaks=s.getCorrPeaks(cc[i],5,pix[pair1],pix[pair2],new float[]{predxshift,predyshift}); //these are phasecorr,x,y,crosscorr
				//peaks are ranked in order of decreasing correlation
				//float[] bestpeak=peaks[0];
				float[][] bestpeak=findBestPeak(peaks,new float[]{predxshift,predyshift},peaksel);
				IJ.log("chose peak "+((int)bestpeak[1][0]+1));
				for(int j=0;j<peaks.length;j++){
					IJ.log(""+peaks[j][4]+" , "+peaks[j][0]+" , "+peaks[j][1]+" , "+peaks[j][2]);
				}
				pairStats[i][0]=bestpeak[0][1]; pairStats[i][1]=bestpeak[0][2]; pairStats[i][2]=bestpeak[0][0]; pairStats[i][3]=bestpeak[0][4];  //this is x,y,phasecorr,crosscorr
				stats[i][4]=pairStats[i][0]; stats[i][5]=pairStats[i][1];
				stats[i][6]=pairStats[i][2]; stats[i][7]=pairStats[i][3];

			}
			IJ.log("found peaks:");
			for(int i=0;i<pairs.length;i++){
				IJ.log(""+pairs[i][0]+ " <- "+pairs[i][1]+":,"+stats[i][2]+" , "+stats[i][3]+" , "+pairStats[i][0]+ " , "+pairStats[i][1]+" , " + " correlation (PC,R)=" + pairStats[i][2]+" , "+pairStats[i][3]);
			}
			if(outcorr) new ImagePlus("Norm Correlation Images",ccstack).show();
			//check the coordinates
			//all vertical and horizontal absolute offsets should be within ftoler (10%) of their median absolute values
			//if not, set them to the median and zero their correlation values (they are now slave to the other offsets)
			IJ.showStatus("checking shifts");
			Object[] checkStats=s.checkCoords(pairStats,pairs,medshifttoler,corrthresh);
			IJ.log("corrected shifts");
			for(int i=0;i<pairs.length;i++){
				IJ.log(""+pairs[i][0]+ " <- "+pairs[i][1]+":,"+stats[i][2]+" , "+stats[i][3]+" , "+pairStats[i][0]+ " , "+pairStats[i][1]+" , " + " correlation (PC,R)=" + pairStats[i][2]+" , "+pairStats[i][3]);
			}
			float[] medshifts=(float[])checkStats[0];
			//find the first pair with significant correlation
			int fixedindex=0;
			for(int i=0;i<pairs.length;i++){
				if(pairStats[i][3]>corrthresh){
					fixedindex=pairs[i][0];
					IJ.log("fixing frame "+fixedindex);
					break;
				}
			}
			IJ.showStatus("creating alignment");
			//now place them on the grid with the fixed image as 0,0
			float[][] newcoords=stitching.preAlign(pairStats,fixedindex,pairs,coords[0].length,corrthresh);
			float[] newxvals=new float[newcoords.length];
			float[] newyvals=new float[newcoords.length];
			float xmin=newcoords[0][0]; float ymin=newcoords[0][1];
			for(int i=0;i<newcoords.length;i++){
				newxvals[i]=newcoords[i][0];
				newyvals[i]=newcoords[i][1];
				if(newxvals[i]<xmin) xmin=newxvals[i];
				if(newyvals[i]<ymin) ymin=newyvals[i];
			}
			for(int i=0;i<newcoords.length;i++){
				newxvals[i]-=xmin; newyvals[i]-=ymin;
			}
			s.setCoords(newxvals,newyvals,1.0f);
			coordplot.addPoints(newxvals,newyvals,true);
			float[] overlaps=getAvgOverlap(new float[][]{newxvals,newyvals},width,height);
			float hoverlap=overlaps[0]; float voverlap=overlaps[1];
			IJ.log("overlaps: "+hoverlap+" , "+voverlap);
			IJ.showStatus("fusing images");
			stitched=s.stitch_frame(pix,true,hoverlap,voverlap);
		}
		ImageStack outstack=new ImageStack(s.newwidth,s.newheight);
		outstack.addSlice("",stitched);
		new ImagePlus("Stitched Image",outstack).show();
		IJ.showStatus("done stitching");
	}

	public static float[] getAvgOverlap(float[][] coords,int width,int height){
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
		hover/=(float)nhpairs;
		vover/=(float)nvpairs;
		return new float[]{(float)width-hover,(float)height-vover};
	}

	public float[][] findBestPeak(float[][] peaks,float[] shift,int method){
		if(method==2){
			//find closest distance
			float mindist=getdist(new float[]{peaks[0][1],peaks[0][2]},shift);
			int minindex=0;
			for(int i=1;i<peaks.length;i++){
				float dist=getdist(new float[]{peaks[i][1],peaks[i][2]},shift);
				if(dist<mindist){
					mindist=dist; minindex=i;
				}
			}
			return new float[][]{peaks[minindex],new float[]{(float)minindex}};
		} else if(method==0){
			//find best crosscorr
			float maxcorr=peaks[0][3];
			int minindex=0;
			for(int i=1;i<peaks.length;i++){
				if(peaks[i][3]>maxcorr){
					maxcorr=peaks[i][3];
					minindex=i;
				}
			}
			return new float[][]{peaks[minindex],new float[]{(float)minindex}};
		} else {
			//find best phasecorr
			float maxcorr=peaks[0][0];
			int minindex=0;
			for(int i=1;i<peaks.length;i++){
				if(peaks[i][0]>maxcorr){
					maxcorr=peaks[i][0];
					minindex=i;
				}
			}
			return new float[][]{peaks[minindex],new float[]{(float)minindex}};
		}
	}

	public float getdist(float[] coords1,float[] coords2){
		return (float)Math.sqrt((coords1[0]-coords2[0])*(coords1[0]-coords2[0])+(coords1[1]-coords2[1])*(coords1[1]-coords2[1]));
	}

	public void showMessage(String message){
		IJ.log(message);
	}

	public void showProgress(int currpos,int finalpos){
		IJ.showProgress(currpos,finalpos);
	}

}
