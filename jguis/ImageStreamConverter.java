package jguis;

import java.awt.Polygon;

import ij.ImagePlus;
import jalgs.algutils;
import jalgs.interpolation;
import jalgs.jstatistics;
import jalgs.jseg.findblobs3;
import jalgs.jseg.jsmooth;
import jalgs.jseg.jsobel;
import jalgs.jseg.measure_object;

public class ImageStreamConverter{

	public static void main(String[] args){
		//this is a stand-alone program to convert ImageStream images to FCS files
		//need to decide which parameters get recorded
		//the first arg is the file path, the second is the fcs path, the third is the mask/trans channel number (consecutive)
		LOCI_series_reader lsr=new LOCI_series_reader(args[0],false);
		int nrows=lsr.nseries/2;
		int maskch=Integer.parseInt(args[2])-1; int transch=maskch;
		ImagePlus imp1=(ImagePlus)lsr.getNextFrame();
		ImagePlus imp2=(ImagePlus)lsr.getNextFrame();
		int nch=imp1.getStackSize();
		int[] outchans=new int[nch];
		for(int i=0;i<nch;i++) outchans[i]=i; //by default get all of the channel parameters
		//int[] outchans={0,6}; //this just gets the meos channels
		String[] collabels=getColLabels(outchans);
		//int ncols=collabels.length;
		float[][] analysis=new float[nrows][];
		System.out.println("image "+0+" of "+nrows+" converting");
		for(int i=0;i<nrows;i++){
			if(i!=0) {
				imp1=(ImagePlus)lsr.getNextFrame();
				imp2=(ImagePlus)lsr.getNextFrame();
			}
			if((i%100)==0) System.out.print("\rimage "+(i+1)+" of "+nrows+" converting");
			Object[] temp=analyzeImage(imp1,imp2,maskch,transch,null,true,0.0f,50,true);
			float[] params=(float[])temp[0];
			params=(float[])algutils.combine_arrays(new float[]{i},params);
			analysis[i]=params;
		}
		lsr.dispose();
		boolean temp=(new export_flowcyte()).write_table(analysis,collabels,args[1]);
		if(!temp) System.out.println("Data Export Error");
		return;
	}
	
	public static String[] getColLabels(int[] outchans){
		String[] statlabels= {"avg","std","min","max","punctaarea","punctaavg","punctacirc"};
		String[] labels= {"object","area","transavg","gradcv","gradavg","gradmax","major","aspectratio","circ"};
		String[] combined=new String[labels.length+outchans.length*statlabels.length];
		for(int i=0;i<labels.length;i++) combined[i]=labels[i];
		for(int j=0;j<outchans.length;j++){
    		for(int i=0;i<statlabels.length;i++) {
    			combined[labels.length+i+j*statlabels.length]="ch"+(outchans[j]+1)+statlabels[i];
    		}
		}
		return combined; 
	}
	
	public static Object[] analyzeImage(ImagePlus imp1,ImagePlus imp2,int maskch,int transch,int[] outchans,boolean backsub,float fillval,int maxdim,boolean noimages){
		//produces a background subtracted (optional) multi-channel image of maxdim x maxdim dimension with the last channel being the mask
		Object[] stack1=jutils.stack2array(imp1.getStack());
		int width=imp1.getWidth(); int height=imp1.getHeight();
		Object[] stack=stack1;
		if(outchans!=null) {
			stack=new Object[outchans.length];
			for(int i=0;i<outchans.length;i++) {
				stack[i]=stack1[outchans[i]];
			}
		}
		Object[] masks=jutils.stack2array(imp2.getStack());
		findblobs3 fb=new findblobs3(width,height);
		float[] object=fb.dofindblobs((byte[])masks[maskch]);
		float[] back=new float[object.length];
		for(int i=0;i<back.length;i++) if(object[i]==0.0f) back[i]=1.0f;
		Polygon outline=fb.get_object_outline(object,1);
		float[] grmsd=gradRMSD(stack1[transch],outline,width,height); //5 params with area first
		float[] ellipse=ellipseParams(outline); //3 parameters with circ last
		float[] params=(float[])algutils.combine_arrays(grmsd,ellipse);
		float[] backs=new float[stack.length];
		for(int i=0;i<stack.length;i++) backs[i]=fb.get_object_stats(back,1,stack[i],"Avg");
		//float maskarea=fb.get_object_stats(object,1,stack[0],"Count");
		int nstats=7;
		int offset=params.length;
		float[] temp=(float[])algutils.expand_array(params,offset+stack.length*nstats);
		float divbystdev=8.0f; float rollballrad=3.0f; float procthresh=0.5f;
		float[][] punctastats=getPunctaStats(stack,object,width,height,divbystdev,rollballrad,procthresh,back);
		for(int i=0;i<stack.length;i++){
			float avg=fb.get_object_stats(object,1,stack[i],"Avg");
			temp[i*nstats+offset]=avg-backs[i];
			float std=fb.get_object_stats(object,1,stack[i],"stdev");
			float min=fb.get_object_stats(object,1,stack[i],"min");
			float max=fb.get_object_stats(object,1,stack[i],"max");
			temp[i*nstats+offset+1]=std; temp[i*nstats+offset+2]=min-backs[i]; temp[i*nstats+offset+3]=max-backs[i];
			temp[i*nstats+offset+4]=punctastats[i][0]; temp[i*nstats+offset+5]=punctastats[i][1]; temp[i*nstats+offset+6]=punctastats[i][2];
		}
		if(!noimages){
			float[][] images=new float[outchans.length+1][maxdim*maxdim];
			float xshift=0.5f*(float)(maxdim-width);
			float yshift=0.5f*(float)(maxdim-height);
			for(int i=0;i<outchans.length;i++){
				for(int j=0;j<maxdim*maxdim;j++) images[i][j]=fillval;
				int selch=outchans[i];
				float[] tempch=algutils.convert_arr_float(stack1[selch]);
				if(backsub) for(int j=0;j<tempch.length;j++) tempch[j]-=backs[selch];
				interpolation.shift_copy_image(tempch,width,height,images[i],maxdim,maxdim,xshift,yshift);
			}
			for(int j=0;j<maxdim*maxdim;j++) images[outchans.length][j]=fillval;
			interpolation.shift_copy_image(object,width,height,images[outchans.length],maxdim,maxdim,xshift,yshift);
			return new Object[]{temp,images};
		} else {
			return new Object[]{temp};
		}
	}

	public static float[] ellipseParams(Polygon outline){
		//find the centroid
		//float[] centroid=measure_object.centroid(outline);
		if(outline==null) return new float[3];
		float[] ellipseparams=measure_object.get_ellipse_parameters(outline);
		//these are 0x,1y,2angle,3major,4minor
		//change minor to aspect ratio
		ellipseparams[4]/=ellipseparams[3];
		float circularity=measure_object.circularity(outline);
		float[] params=(float[])algutils.expand_array(ellipseparams,ellipseparams.length+1);
		params[params.length-1]=circularity;
		float[] params2= {params[3],params[4],params[5]}; //these are major, aspect, circularity
		return params2;
	}

	public static float[] gradRMSD(Object pix,Polygon outline1,int width,int height){
		float[] pix2=algutils.convert_arr_float2(pix);
		float[] sobel=(new jsobel(width,height)).do_sobel(pix2)[0];
		Polygon outline=new Polygon(new int[]{0,width-1,width-1,0},new int[]{0,0,height-1,height-1},4);
		if(outline1!=null) outline=outline1;
		float cnt=jstatistics.getstatistic("Count",pix,width,height,outline,null);
		float avg=jstatistics.getstatistic("Avg",pix,width,height,outline,null);
		float sobelcv=jstatistics.getstatistic("stdev",sobel,width,height,outline,null);
		sobelcv/=avg;
		float sobelmax=jstatistics.getstatistic("Max",pix,width,height,outline,null);
		float sobelavg=jstatistics.getstatistic("Avg",pix,width,height,outline,null);
		//note pixel size is around 330 nm per pixel (area is 0.11 um^2 per pixel)
		return new float[]{cnt*0.11f,avg,sobelcv,sobelavg,sobelmax};
	}
	
	public static float[][] getPunctaStats(Object[] pix,float[] object,int width,int height,float divbystdev,float rollballrad,float procthresh,float[] back){
		//for each pix image, get the puncta area, intensity, and circularity
		float[][] stats=new float[pix.length][3];
		for(int i=0;i<pix.length;i++) {
			float[] temp=algutils.convert_arr_float(pix[i]);
			for(int j=0;j<temp.length;j++) temp[j]-=back[i];
			float[] blurred=temp.clone();
			//blur the image
			jsmooth.blur2D(temp,divbystdev,width,height);
			//divide by blurred to eliminate the variation in puncta intensity
			//avoid dividing by 0
			//put the original image back in temp for reference
			for(int j=0;j<width*height;j++){
				float temp2=blurred[j];
				if(temp[j]<0.1f) blurred[j]/=temp[j];
				temp[j]=temp2;
			}
			//subtract a rolling ball
			blurred=jutils.sub_roll_ball_back(blurred,rollballrad,width,height);
			//now find the avg and area of regions above procthresh (and inside the object)
			float[] obj=new float[width*height];
			for(int j=0;j<width*height;j++) {
				if(object[j]>0.0f && blurred[j]>procthresh) {
					stats[i][0]+=1.0f;
					obj[j]=1.0f;
					stats[i][1]+=temp[j];
				}
			}
			if(stats[i][0]>0.0f) { //find puncta parameters only if there is a punctum
				stats[i][1]/=stats[i][0];
    			findblobs3 fb=new findblobs3(width,height);
    			obj=fb.dofindblobs(obj,0.5f);
    			float[][] apc=fb.get_area_perim_circ(obj);
    			float avgcirc=jstatistics.getstatistic("Avg",apc[2],null);
    			stats[i][2]=avgcirc;
			}
		}
		return stats;
	}

}
