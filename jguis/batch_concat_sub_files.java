package jguis;

import java.awt.Color;
import java.io.File;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.IJ;
import ij.ImagePlus;
import ij.process.LUT;
import jalgs.algutils;
import jalgs.jdataio;

public class batch_concat_sub_files {
	//this class batch processes screening image concatenation with projection and subtraction

	public static void main(String[] args) {
		//args are indir,outdir,platename,outplatename,channels,dapichan, slices, threads
		//hard code colors
		//typical channels is 4, dapichan is 4, slices is 5 or 6, threads is 4
		//for this version don't save the full image stack (saveStack is false)
		int chan=Integer.parseInt(args[4]);
		int dapichan=Integer.parseInt(args[5]);
		int slices=Integer.parseInt(args[6]);
		int threads=Integer.parseInt(args[7]);
		exec(args[0],args[1],args[2],args[3],null,chan,dapichan,slices,threads,false);
	}
	
	public static void exec(String indir,String outdir,String pname,String outplate,Color[] colors,
			int channels,int dapichan,int slices,int nthreads,boolean saveStack) {
		String stackdir=outdir+outplate+File.separator+"analysis"+File.separator+"stacks"+File.separator;
		String dapidir=outdir+outplate+File.separator+"analysis"+File.separator+"dapi"+File.separator;
		String measdir=outdir+outplate+File.separator+"analysis"+File.separator+"measurement_channels"+File.separator;
		new File(stackdir).mkdirs();
		new File(dapidir).mkdirs();
		new File(measdir).mkdirs();
		int start_row=1;
		int end_row=8;
		int start_col=1;
		int end_col=12;
		ExecutorService executor=Executors.newFixedThreadPool(nthreads);
		for(int i=start_row;i<=end_row;i++) {
			for(int j=start_col;j<=end_col;j++) {
				Runnable worker=new wellRunner(indir,outdir,pname,outplate,
						i,j,-1,slices,channels,dapichan,colors,saveStack);
				//processWell(indir,outdir,pname,outplate,i,j,-1,slices,channels,dapichan,colors);
				executor.execute(worker);
			}
		}
		executor.shutdown();
		while(!executor.isTerminated()) {}
		System.out.println("process finished for "+pname);
	}
	
	public static String letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
	public static void processWell(String indir,String outdir,String pname,String outplate,
			int row,int col,int nframes,int nslices,int nchan,int dapichan,Color[] colors,
			boolean saveStack) {
		String rowstring="r0"+row;
		String colstring="c";
		if(col>9) colstring+=""+col;
		else colstring+="0"+col;
		String wellstring=rowstring+colstring;
		String new_wellstring="_well_"+letters.charAt(row-1)+colstring.substring(1);
		jdataio jdio=new jdataio();
		System.out.println(wellstring);
		System.out.println(new_wellstring);
		String tdir=indir+pname+File.separator+"Images"+File.separator;
		//check if our well exists
		if(!(new File(tdir+wellstring+"f01p01-ch1sk1fk1fl1.tiff")).exists()) return;
		String[] fnames=jdio.get_numeric_sorted_string_list(tdir,wellstring);
		//get the number of frames
		int tnframes=(int)(fnames.length/(nslices*nchan));
		if(nframes>0) tnframes=nframes;
		System.out.println("n frames: "+nframes);
		//check if we have the right number of frames
		if(fnames.length!=tnframes*nslices*nchan) return;
		Object[] slicepix=new Object[fnames.length];
		int width=-1;
		int height=-1;
		for(int i=0;i<fnames.length;i++) {
			ImagePlus imp=IJ.openImage(tdir+fnames[i]);
			if(width<0) {
				width=imp.getWidth();
				height=imp.getHeight();
			}
			//slicepix[i]=jutils.stack2array(imp.getImageStack());
			slicepix[i]=imp.getProcessor().getPixels();
			//if(i%100==0) System.out.print("frame "+i+" read\r");
		}
		System.out.println("finished "+wellstring+" read");
		Color[] tcolors=new Color[] {Color.RED, Color.GREEN,Color.MAGENTA, Color.BLUE};
		if(colors!=null) tcolors=colors;
		LUT[] luts=new LUT[nchan];
		for(int i=0;i<nchan;i++) {
			luts[i]=jutils.get_lut_for_color(tcolors[i]);
		}
		//make the hyperstack
		ImagePlus imp=jutils.create_hyperstack(outplate+new_wellstring+".tif", 
				jutils.array2stack(slicepix,width,height),tnframes,nslices,nchan,true,luts);
		//and save it
		String stackdir=outdir+outplate+File.separator+"analysis"+File.separator+"stacks"+File.separator;
		String newname=outplate+new_wellstring+".tif";
		//new Tiff_Writer(imp,fnames.length,null).saveAsTiffStack2(stackdir+newname);
		if(saveStack) IJ.saveAsTiff(imp, stackdir+newname);
		
		//now make the dapi max projection
		Object[] dapimax=maxProjDAPI(slicepix,tnframes,nslices,nchan,dapichan-1);
		//wrap it and save it
		String dapidir=outdir+outplate+File.separator+"analysis"+File.separator+"dapi"+File.separator;
		ImagePlus dapiimp=jutils.create_hyperstack("dapi", 
				jutils.array2stack(dapimax,width,height),tnframes,1,1,false,null);
		//new Tiff_Writer(dapiimp,dapimax.length,null).saveAsTiffStack2(dapidir+newname);
		IJ.saveAsTiff(dapiimp, dapidir+newname);
		
		//now make the sum projection
		float[][] sumproj=sumProjStack(slicepix,tnframes,nslices,nchan);
		//run the rolling ball sub
		float[][] sumproj2=jutils.sub_roll_ball_back(sumproj, 1000.0f, width, height);
		sumproj=null;
		//wrap it and save it
		String measdir=outdir+outplate+File.separator+"analysis"+File.separator+"measurement_channels"+File.separator;
		ImagePlus measimp=jutils.create_hyperstack("measurement", 
				jutils.array2stack(sumproj2,width,height),tnframes,1,nchan,true,null);
		//new Tiff_Writer(measimp,sumproj2.length,null).saveAsTiffStack2(measdir+newname);
		IJ.saveAsTiff(measimp, measdir+newname);
		System.out.println("finished well "+wellstring);
		return;
	}
	
	public static float[][] sumProjStack(Object[] stack,int nframes,int nslices,int nchan){
		float[][] sumproj=new float[nframes*nchan][];
		for(int i=0;i<nframes;i++) {
			for(int j=0;j<nchan;j++) {
				int pos=i*nslices*nchan+j;
				float[] tproj=algutils.convert_arr_float(stack[pos]);
				for(int k=1;k<nslices;k++) {
					pos=i*nslices*nchan+j+k*nchan;
					float[] temp=algutils.convert_arr_float(stack[pos]);
					for(int l=0;l<temp.length;l++) tproj[l]+=temp[l];
				}
				sumproj[j+i*nchan]=tproj;
			}
		}
		return sumproj;
	}
	
	public static Object[] maxProjDAPI(Object[] stack,int nframes,int nslices,int nchan,int dapichan){
		Object[] maxproj=new Object[nframes];
		for(int i=0;i<nframes;i++) {
			int pos=i*nslices*nchan+dapichan;
			Object tproj=stack[pos];
			if(tproj instanceof float[]) {
				float[] ftproj=((float[])tproj).clone();
				for(int k=1;k<nslices;k++) {
					pos=i*nslices*nchan+dapichan+k*nchan;
					float[] temp=(float[])stack[pos];
					for(int l=0;l<temp.length;l++) ftproj[l]=fmax(ftproj[l],temp[l]);
				}
				maxproj[i]=ftproj;
			} else if(tproj instanceof short[]) {
				short[] ftproj=((short[])tproj).clone();
				for(int k=1;k<nslices;k++) {
					pos=i*nslices*nchan+dapichan+k*nchan;
					short[] temp=(short[])stack[pos];
					for(int l=0;l<temp.length;l++) ftproj[l]=smax(ftproj[l],temp[l]);
				}
				maxproj[i]=ftproj;
			} else { //this is byte
				byte[] ftproj=((byte[])tproj).clone();
				for(int k=1;k<nslices;k++) {
					pos=i*nslices*nchan+dapichan+k*nchan;
					byte[] temp=(byte[])stack[pos];
					for(int l=0;l<temp.length;l++) ftproj[l]=bmax(ftproj[l],temp[l]);
				}
				maxproj[i]=ftproj;
			}
		}
		return maxproj;
	}
	
	public static float fmax(float v1,float v2) {
		if(v1>v2) return v1;
		else return v2;
	}
	
	public static short smax(short v1,short v2) {
		if(v1>v2) return v1;
		else return v2;
	}
	
	public static byte bmax(byte v1,byte v2) {
		if(v1>v2) return v1;
		else return v2;
	}

}

class wellRunner implements Runnable{
	private String indir,outdir,pname,outplate;
	private int row,col,nframes,nslices,nchan,dapichan;
	private Color[] colors;
	private boolean saveStack;
	
	public wellRunner(String indir,String outdir,String pname,String outplate,
			int row,int col,int nframes,int nslices,int nchan,int dapichan,Color[] colors,
			boolean saveStack) {
		this.indir=indir;
		this.outdir=outdir;
		this.pname=pname;
		this.outplate=outplate;
		this.row=row;
		this.col=col;
		this.nframes=nframes;
		this.nslices=nslices;
		this.nchan=nchan;
		this.colors=colors;
		this.dapichan=dapichan;
	}
	
	public void run() {
		batch_concat_sub_files.processWell(indir,outdir,pname,outplate,
				row,col,nframes,nslices,nchan,dapichan,colors,saveStack);
	}
}
