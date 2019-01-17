package jguis;

import ij.CompositeImage;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.RGBStackConverter;
import ij.process.ColorProcessor;
import ij.process.LUT;
import ij.process.StackConverter;
import jalgs.algutils;
import jalgs.jstatistics;

import java.io.File;

public class CLI_file_converter{

	public static void main(String[] args){
		//here we want to max project the first frame and convert to RGB
		//the first arg is the incoming file name, the second is the saving file name
		LOCI_file_reader lfr=new LOCI_file_reader();
		ImagePlus imp=lfr.get_loci_imp(args[0]);
		if(imp==null){
			System.out.println("image did not load");
			return;
		}
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nchan=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		//now do the z projection over all the channels in the first frame if necessary
		float[][] proj=new float[nchan][];
		for(int i=0;i<nchan;i++){
			proj[i]=jutils.get3DProjZStat(stack,0,i,nframes,nslices,nchan,"Max");
		}
		//if there are more than 4 channels, do the projection over the channels as well and output a grayscale image
		/*LUT[] lut=((CompositeImage)imp).getLuts();
		for(int i=0;i<lut.length;i++){
			int index=jutils.get_LUT_color_index(lut[i]);
			System.out.println(jutils.colornames[index]);
		}*/
		int[] colorpix=new int[width*height];
		if(nchan>4){
			float[] proj2=algutils.get_stack_proj_stat("Avg",proj,width,height,nslices,null);
			float min=jstatistics.getstatistic("Min",proj2,null);
			float max=jstatistics.getstatistic("Max",proj2,null);
			for(int i=0;i<proj2.length;i++){
				int intval=(int)(255.0f*(proj2[i]-min)/(max-min));
				colorpix[i]=jutils.rgb2intval(intval,intval,intval);
			}
			/*ImagePlus colorimp=new ImagePlus("temp",new ColorProcessor(width,height,colorpix));
			FileSaver fs=new FileSaver(colorimp);
			fs.saveAsTiff(args[1]);
			imp.close();
			System.out.println(args[0]+"=>"+args[1]+" complete");*/
		} else { //try to use the saved luts
			/*byte[][] chanstack=new byte[nchan][width*height];
			for(int i=0;i<nchan;i++){
				float min=jstatistics.getstatistic("Min",proj[i],null);
				float max=jstatistics.getstatistic("Max",proj[i],null);
				for(int j=0;j<proj[i].length;j++){
					int intval=(int)(255.0f*(proj[i][j]-min)/(max-min));
					chanstack[i][j]=(byte)intval;
				}
			}*/
			ImageStack chanstack2=jutils.array2stack(proj,width,height);
			ImagePlus chanimp=jutils.create_hyperstack("chanimp",chanstack2,imp,1,1,nchan);
			((CompositeImage)chanimp).setDisplayMode(CompositeImage.COMPOSITE);
			/*FileSaver fs=new FileSaver(chanimp);
			fs.saveAsTiff(args[1]);
			imp.close();
			System.out.println(args[0]+"=>"+args[1]+" complete");*/
			//StackConverter sc=new StackConverter(chanimp);
			//sc.convertToRGB();
			RGBStackConverter.convertToRGB(chanimp);
			colorpix=(int[])chanimp.getProcessor().getPixels();
			chanimp.close();
		}
		ImagePlus colorimp=new ImagePlus("temp",new ColorProcessor(width,height,colorpix));
		FileSaver fs=new FileSaver(colorimp);
		fs.saveAsTiff(args[1]);
		imp.close();
		System.out.println(args[0]+"=>"+args[1]+" complete");
	}

}
