package jguis;

import loci.common.services.ServiceFactory;
import loci.formats.ChannelSeparator;
import loci.formats.FormatTools;
import loci.formats.FormatWriter;
import loci.formats.IFormatWriter;
import loci.formats.ImageWriter;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.formats.out.APNGWriter;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import ij.ImagePlus;

public class LOCI_file_writer{

	public boolean save_imp_as_png(ImagePlus imp,String path){
		//IFormatWriter w=new APNGWriter();
		try{
    		IFormatWriter w = new ImageWriter().getWriter(path);
    		int[] pix=(int[])imp.getProcessor().getPixels();
    		/*byte[] pix2=new byte[pix.length*3];
    		for(int i=0;i<pix.length;i++){
    			pix2[i*3+2]=(byte)pix[i];
    			pix2[i*3+1]=(byte)(pix[i]>>8);
    			pix2[i*3]=(byte)(pix[i]>>16);
    			//pix2[i*4]=(byte)(pix[i]>>24);
    		}*/
    		byte[][] pixtemp=jutils.intval2rgb(pix);
    		byte[] pix2=new byte[pix.length*3];
    		System.arraycopy(pixtemp[0],0,pix2,0,pix.length);
    		System.arraycopy(pixtemp[1],0,pix2,pix.length,pix.length);
    		System.arraycopy(pixtemp[2],0,pix2,2*pix.length,pix.length);
    
    	    ServiceFactory factory = new ServiceFactory();
    	    OMEXMLService service = factory.getInstance(OMEXMLService.class);
    	    IMetadata meta = service.createOMEXMLMetadata();
    	    int pixelType = FormatTools.UINT8;
    	    MetadataTools.populateMetadata(meta, 0, null, false, "XYZCT",
    	      FormatTools.getPixelTypeString(pixelType), imp.getWidth(), imp.getHeight(), 1, 3, 1, 3);
    	    w.setMetadataRetrieve(meta);
    	    w.setId(path);
    		w.saveBytes(0,pix2);
    		w.close();
		} catch (Exception e){
			return false;
		}
		return true;
	}

}
