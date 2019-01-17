package jguis;

import ij.IJ;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.RoiDecoder;
import ij.io.RoiEncoder;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

public class multi_roi_writer{
	
	public static boolean writeRois(Roi[] rois,String path){
		if(rois.length<1) return false;
		try{
			ZipOutputStream zos=new ZipOutputStream(new FileOutputStream(path));
			DataOutputStream out=new DataOutputStream(new BufferedOutputStream(zos));
			RoiEncoder re=new RoiEncoder(out);
			for(int i=0;i<rois.length;i++){
				if(rois[i]!=null){
					String label=rois[i].getName();
					if(label==null) label=getLabel(rois[i]);
					if(!label.endsWith(".roi")) label+=".roi";
					zos.putNextEntry(new ZipEntry(label));
					re.write(rois[i]);
					out.flush();
				}
			}
			out.close();
			return true;
		} catch(IOException e){
			IJ.log(e.toString());
			return false;
		}
	}
	
	public static String getLabel(Roi roi) {
		Rectangle r=roi.getBounds();
		int xc=r.x+r.width/2;
		int yc=r.y+r.height/2;
		if(xc<0) xc=0;
		if(yc<0) yc=0;
		String label=padIntString(yc,4)+"-"+padIntString(xc,4);
		if(xc>=10000) label=padIntString(yc,5)+"-"+padIntString(xc,5);
		return label;
	}
	
	public static String padIntString(int val,int ndigits) {
		int val2=val;
		if(val2<0) val2=0;
		String sval=""+val2;
		int slength=sval.length();
		if(slength>=ndigits) {
			return sval.substring(0,ndigits);
		} else {
			for(int i=slength;i<ndigits;i++) {
				sval="0"+sval;
			}
			return sval;
		}
	}
	
	public static boolean writeRois(Polygon[] polys,String path) {
		Roi[] rois=new Roi[polys.length];
		for(int i=0;i<polys.length;i++) rois[i]=new PolygonRoi(polys[i],Roi.FREEROI);
		return writeRois(rois,path);
	}
	
	public static boolean writeRoi(Roi roi,String path){
		RoiEncoder re=new RoiEncoder(path);
		try{
			re.write(roi);
			return true;
		} catch(IOException e){
			IJ.log(e.toString());
			return false;
		}
	}
	
	public static Roi readRoi(String path){
		RoiDecoder rd=new RoiDecoder(path);
		try{
			return rd.getRoi();
		} catch(IOException e){
			IJ.log(e.toString());
			return null;
		}
	}
	
	public static Roi[] readRois(String path){
		ZipInputStream in=null;
		ByteArrayOutputStream out;
		try{
			in=new ZipInputStream(new FileInputStream(path));
			int nRois=0;
			while(in.getNextEntry()!=null) nRois++;
			in.close();
			Roi[] rois=new Roi[nRois];
			in=new ZipInputStream(new FileInputStream(path));
			byte[] buf = new byte[1024]; 
            int len;
            int counter=0;
			ZipEntry entry=in.getNextEntry();
			while(entry!=null){
				String name=entry.getName();
				out=new ByteArrayOutputStream();
				while((len=in.read(buf))>0) out.write(buf,0,len);
				out.close();
				byte[] data=out.toByteArray();
				rois[counter]=(new RoiDecoder(data,name)).getRoi();
			}
			in.close();
			return rois;
		} catch(IOException e){
			IJ.log(e.toString());
			return null;
		}
	}
	
	public static Roi[] makeRois(float[][] coords){
		Roi[] rois=new Roi[coords.length];
		for(int i=0;i<coords.length;i++){
			int x=(int)coords[i][0];
			int y=(int)coords[i][1];
			int z=1+(int)coords[i][2];
			Roi temp=new PointRoi(x,y);
			temp.setName(getRoiLabel(temp,z));
			rois[i]=temp;
		}
		return rois;
	}
	
	public static String getRoiLabel(Roi roi,int slice){
		//basically this is slice-yc-xc with at minimum 4 digits
		//note that slice is the slice index (channel+slice*channels+frame*slices*channels+1), not the actual z slice
		Rectangle r=roi.getBounds();
		int xc=r.x+r.width/2;
		int yc=r.y+r.height/2;
		String xval="0000"+xc;
		xval=xval.substring(xval.length()-4);
		String yval="0000"+yc;
		yval=yval.substring(yval.length()-4);
		String zval="0000"+slice;
		zval=zval.substring(zval.length()-4);
		return zval+"-"+yval+"-"+xval;
	}
	
}
