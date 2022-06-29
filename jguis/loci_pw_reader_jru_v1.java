package jguis;

import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;

import java.io.File;
import java.io.IOException;

import ome.units.UNITS;
import ome.units.quantity.Length;
import ome.units.quantity.Time;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.meta.MetadataStore;

public class loci_pw_reader_jru_v1 extends FormatReader implements ImageListener{
	public String path,fname;
	public Object plot;
	public int plotindex;

	public loci_pw_reader_jru_v1(){
		super("Plot Window",new String[]{"pw","pw2"});
		suffixSufficient=true;
		suffixNecessary=true;
		domains=new String[]{FormatTools.UNKNOWN_DOMAIN};
	}

	public boolean isThisType(RandomAccessInputStream stream) throws IOException{
		return false;
	}

	public byte[] openBytes(int no,byte[] buf,int x,int y,int w,int h){
		//IJ.log("attempting open");
		//for(int i=0;i<buf.length;i++) buf[i]=(byte)255;
		return buf;
	}

	public void close(boolean fileOnly) throws IOException{
		super.close(fileOnly);
		plot=null;
		path=null;
	}
	
	protected void initFile(String id) throws FormatException, IOException{
		//IJ.log("initializing superclass");
		super.initFile(id);
		//IJ.log("initializing subclass");
		//path=Location.getMappedId(id);
		path=new Location(id).getAbsolutePath();
		fname=(new File(path)).getName();
		//path=id;
		//IJ.log("location id path:"+path);
		//CoreMetadata cm=core.get(0);
		//CoreMetadata cm=super.getCoreMetadataList().get(0);
		//CoreMetadata cm=getCoreMetadata()[0];
		CoreMetadata cm=getCoreMetadataList().get(0);
		cm.sizeX=Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN;
		cm.sizeY=Plot4.HEIGHT+Plot4.BOTTOM_MARGIN+Plot4.TOP_MARGIN;
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
		cm.pixelType=FormatTools.UINT8;
		cm.thumbnail=false;
		cm.bitsPerPixel=8;
		MetadataStore store = makeFilterMetadata();
		MetadataTools.populatePixels(store, this);
		store.setImageName(path,0);
		store.setPixelsPhysicalSizeX(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsPhysicalSizeY(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsPhysicalSizeZ(new Length(1.0,UNITS.MICROMETER),0);
		store.setPixelsTimeIncrement(new Time(1.0,UNITS.SECOND),0);
		if(path.endsWith(".pw") || Plot4.is_this(path)){
			plot=new Plot4(path); plotindex=0;
			(new PlotWindow4(fname+"temp",(Plot4)plot)).draw();
		} else if(Plot3D.is_this(path)){
			plot=new Plot3D(path); plotindex=1;
			(new PlotWindow3D(fname+"temp",(Plot3D)plot)).draw();
		} else if(Traj3D.is_this(path)){
			plot=new Traj3D(path); plotindex=2;
			(new PlotWindow3D(fname+"temp",(Traj3D)plot)).draw();
		} else if(PlotHist.is_this(path)){
			plot=new PlotHist(path); plotindex=3;
			(new PlotWindowHist(fname+"temp",(PlotHist)plot)).draw();
		} else if(Plot2DHist.is_this(path)){
			plot=new Plot2DHist(path); plotindex=4;
			(new PlotWindow2DHist(fname+"temp",(Plot2DHist)plot)).draw();
		} else {
			plot=null;
			throw new FormatException("unsupported or corrupted pw2 file");
		}
		if(plot!=null) ImagePlus.addImageListener(this);
		else IJ.log("plot failed to load");
		//if(plot!=null) IJ.log("populating image");
		//IJ.runPlugIn("import_plot_jru_v1",path);
	}
	
	public void imageOpened(ImagePlus imp){
		//IJ.log("rendering plot");
		//String dir=imp.getOriginalFileInfo().directory;
		//String tpath=dir+imp.getTitle();
		//IJ.log(tpath);
		//IJ.log(path);
		if(imp.getTitle().equalsIgnoreCase(fname)){
			/*if(plotindex==0){
				//new PlotWindow4(imp,(Plot4)plot).draw();
				new PlotWindow4(imp.getTitle(),(Plot4)plot).draw();
			} else if(plotindex>1 && plotindex<4){
				new PlotWindow3D(imp,(Plot3D)plot).draw();
			} else if(plotindex==4){
				new PlotWindowHist(imp,(PlotHist)plot).draw();
			} else {
				new PlotWindow2DHist(imp,(Plot2DHist)plot).draw();
			}
			imp.close();*/
			imp.close();
			ImagePlus.removeImageListener(this);
			ImagePlus temp=WindowManager.getImage(fname+"temp");
			temp.setTitle(fname);
		}
	}

	public void imageClosed(ImagePlus imp){}

	public void imageUpdated(ImagePlus imp){}

}
