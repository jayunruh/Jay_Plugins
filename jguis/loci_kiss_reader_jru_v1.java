package jguis;

import java.io.File;
import java.io.IOException;

import ome.units.quantity.Length;
import ome.units.quantity.Time;
import ome.units.UNITS;
import ij.ImageListener;
import ij.ImagePlus;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.meta.MetadataStore;
import ij.gui.*;
import ij.*;

import java.io.*;

public class loci_kiss_reader_jru_v1 extends FormatReader implements ImageListener {
	public String path,fname;
	public KissPanel sp;

	public loci_kiss_reader_jru_v1(){
		super("KISS Analysis",new String[]{"kiss","sky"});
		suffixSufficient=true;
		suffixNecessary=true;
		domains=new String[]{FormatTools.UNKNOWN_DOMAIN};
	}

	public boolean isThisType(RandomAccessInputStream stream) throws IOException{
		return false;
	}

	public byte[] openBytes(int no,byte[] buf,int x,int y,int w,int h){
		//for(int i=0;i<buf.length;i++) buf[i]=(byte)255;
		return buf;
	}

	public void close(boolean fileOnly) throws IOException{
		super.close(fileOnly);
		sp=null;
		path=null;
	}

	protected void initFile(String id) throws FormatException, IOException{
		super.initFile(id);
		//path=Location.getMappedId(id);
		path=new Location(id).getAbsolutePath();
		fname=(new File(path)).getName();
		//IJ.log(path);
		//CoreMetadata cm=core.get(0);
		CoreMetadata cm=getCoreMetadataList().get(0);
		//CoreMetadata cm=getCoreMetadata()[0];
		cm.sizeX=20;
		cm.sizeY=20;
		cm.sizeZ=1;
		cm.sizeC=1;
		cm.sizeT=1;
		cm.imageCount=1;
		cm.dimensionOrder="XYCZT";
		cm.rgb=true;
		cm.interleaved=true;
		cm.littleEndian=true;
		cm.indexed=false;
		cm.falseColor=false;
		cm.metadataComplete=true;
		cm.pixelType=FormatTools.UINT32;
		cm.thumbnail=false;
		cm.bitsPerPixel=32;
		MetadataStore store = makeFilterMetadata();
		MetadataTools.populatePixels(store, this);
		store.setImageName(path,0);
		store.setPixelsPhysicalSizeX(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsPhysicalSizeY(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsPhysicalSizeZ(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsTimeIncrement(new Time(1.0,UNITS.SECOND),0);
		if(path.endsWith(".kiss") || path.endsWith(".sky")){
			sp=getKISSPanel(path);
			KissPanel.launch_frame(sp);
		} else {
			sp=null;
			throw new FormatException("unsupported or corrupted spectral karyotyping file");
		}
		if(sp!=null) ImagePlus.addImageListener(this);
		//IJ.runPlugIn("import_plot_jru_v1",path);
	}

	public void imageOpened(ImagePlus imp){
		//String dir=imp.getOriginalFileInfo().directory;
		//String tpath=dir+imp.getTitle();
		//IJ.log(tpath);
		if(imp.getTitle().equalsIgnoreCase(fname)){
			//KissPanel.launch_frame(sp);
			imp.close();
			ImagePlus.removeImageListener(this);
		}
	}

	public void imageClosed(ImagePlus imp){}

	public void imageUpdated(ImagePlus imp){}

	public KissPanel getKISSPanel(String path){
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(path));
			KissPanel sp=new KissPanel();
			int nch=5;
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addNumericField("Area Accuracy (percent)",30,0);
			for(int i=0;i<nch;i++){
				gd2.addNumericField("Ch_"+(i+1)+"_Contr_Thresh",0.35,5,15,null);
			}
			//gd2.addNumericField("Contribution Threshold",0.35,5,15,null);
			gd2.addCheckbox("Mouse?",false);
			Object[] codes=KissPanel.getCustomCodes();
			String[] codenames=new String[]{"none"};
			if(codes!=null){
				String[] temp=(String[])codes[0];
				codenames=new String[temp.length+1]; codenames[0]="none";
				for(int i=0;i<temp.length;i++) codenames[i+1]=temp[i];
			}
			gd2.addChoice("Custom_Code",codenames,codenames[0]);
			gd2.addNumericField("Box_Width",150,0);
			gd2.addNumericField("Box_Height",100,0);
			gd2.addCheckbox("Output_Unmixed?",false);
			gd2.showDialog(); if(gd2.wasCanceled()){return null;}
			sp.areathresh=(float)gd2.getNextNumber();
			sp.objthresh2=new float[nch];
			for(int i=0;i<nch;i++) sp.objthresh2[i]=(float)gd2.getNextNumber();
			//sp.objthresh=(float)gd2.getNextNumber();
			boolean mouse=gd2.getNextBoolean();
			int codeindex=gd2.getNextChoiceIndex();
			int bwidth=(int)gd2.getNextNumber();
			int bheight=(int)gd2.getNextNumber();
			boolean outunmixed=gd2.getNextBoolean();
			int[] colorindices={4,1,2,6,3};
			GenericDialog gd3=new GenericDialog("Color Options");
			for(int i=0;i<5;i++) gd3.addChoice("Ch"+(i+1)+" Color",KissPanel.colornames,KissPanel.colornames[colorindices[i]]);
			gd3.showDialog(); if(gd3.wasCanceled()) return null;
			for(int i=0;i<5;i++) colorindices[i]=gd3.getNextChoiceIndex();
			sp.colorindices=colorindices;
			sp.nch=5;
			sp.dapilast=false;
			sp.cellwidth=bwidth;
			sp.cellheight=bheight;
			int[][] custcode=null;
			if(codeindex>0) custcode=(int[][])codes[codeindex+1];
			sp.init(is,mouse,custcode);
			is.close();
			if(outunmixed){
				ImageStack unstack=jutils.array2stack(sp.unmixed,sp.threshimp.getWidth(),sp.threshimp.getHeight());
				ImagePlus unimp=new ImagePlus("Unmixed",unstack);
				unimp.setOpenAsHyperStack(true);
				unimp.setDimensions(sp.unmixed.length,1,1);
				new CompositeImage(unimp,CompositeImage.COLOR).show();
			}
			return sp;
		} catch(IOException e){
			return null;
		}
	}

}