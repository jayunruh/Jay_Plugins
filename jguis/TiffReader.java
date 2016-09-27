package jguis;

import ij.IJ;
import ij.io.RandomAccessStream;
import jalgs.jdataio;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class TiffReader{
	
	//protected RandomAccessStream is;
	public InputStream is;
	public int newsubfiletype,width,height,bitspersample,photointerp,stripoffsets;
    public int orientation,samplesperpixel,xresolution,yresolution,resolutionunit,datetime;
    public int colormap,sampleformat,rowsperstrip,stripbytecount,off;
    public String imgdescription,path;
	
	public TiffReader(String path){
		try{
			this.path=path;
			is=new BufferedInputStream(new FileInputStream(path));
			//this.is=new RandomAccessStream(is1);
		} catch(IOException e){
			e.printStackTrace();
			return;
		}
	}
	
	public void closeStream(){
		try{
			is.close();
		} catch(IOException e){
			e.printStackTrace();
			return;
		}
	}
	
	/************************
	 * this reads the tiff tags from the file
	 * @return
	 */
	public List<List<String>> readTags(){
		jdataio jdio=new jdataio();
		//try{
			//is.seek(0);
			off=0;
			String head=jdio.readstring(is,4); off+=4;
			IJ.log("header = "+head);
			boolean lend=head.indexOf("II")>=0;
			List<List<String>> tags=new ArrayList<List<String>>();
			if(lend){
    			int ifdoff=jdio.readintelint(is); off+=4;
    			IJ.log("ifdoff = "+ifdoff);
    			jdio.skipstreambytes(is,ifdoff-off); off=ifdoff;
    			int nentries=jdio.readintelshort(is); off+=2;
    			IJ.log("nentries = "+nentries);
    			for(int i=0;i<nentries;i++){
    				List<String> entry=readIntelEntry();
    				if(entry==null) return null;
    				tags.add(entry);
    			}
    			int nextifd=jdio.readintelint(is); off+=4;
    			IJ.log("nextifd = "+nextifd);
			} else {
    			int ifdoff=jdio.readmotorolaint(is); off+=4;
    			IJ.log("ifdoff = "+ifdoff);
    			jdio.skipstreambytes(is,ifdoff-off); off=ifdoff;
    			int nentries=jdio.readmotorolashort(is); off+=2;
    			IJ.log("nentries = "+nentries);
    			for(int i=0;i<nentries;i++){
    				List<String> entry=readMotorolaEntry();
    				if(entry==null) return null;
    				tags.add(entry);
    			}
			}
			return tags;
		/*}catch(IOException e){
			e.printStackTrace();
			return null;
		}*/
		
	}
	
	/*********************
	 * reads an intel (little endian) tiff tag
	 * @param off
	 * @return
	 */
	public List<String> readIntelEntry(){
		//try{
			jdataio jdio=new jdataio();
			int tag=jdio.readintelshort(is); off+=2;
			int type=jdio.readintelshort(is); off+=2;
			int cnt=jdio.readintelint(is); off+=4;
			float val=0.0f;
			String temp="";
			if(cnt==1 && type==3){ //padded short
				val=jdio.readintelshort(is); off+=2;
				int zero=jdio.readintelshort(is); off+=2;
			} else if(type==4){ //int
				val=jdio.readintelint(is); off+=4;
			} else if(type==2){ //string (at offset)
				int stroff=jdio.readintelint(is); off+=4;
				temp=readOffsetString(stroff,cnt);
			} else if(type==5){ //rational data type (numerator and denominator ints at offset)
				int numoff=jdio.readintelint(is); off+=4;
				val=readOffsetIntelRational(numoff);
			} else {
				val=jdio.readintelint(is); off+=4;
			}
			List<String> entry=new ArrayList<String>();
			entry.add(""+tag); entry.add(""+type); entry.add(""+cnt);
			if(type==2) entry.add(temp);
			else entry.add(""+val);
			Object tagval=new Integer((int)val);
			if(type==2) tagval=temp;
			entry.add(getTagDesc(tag,tagval));
			return entry;
		/*} catch(IOException e){
			e.printStackTrace();
			return null;
		}*/
	}
	
	/******************
	 * reads a motorola byte order tiff tag
	 * @param off
	 * @return
	 */
	public List<String> readMotorolaEntry(){
		//try{
			jdataio jdio=new jdataio();
			int tag=jdio.readmotorolashort(is); off+=2;
			int type=jdio.readmotorolashort(is); off+=2;
			int cnt=jdio.readmotorolaint(is); off+=4;
			float val=0.0f;
			String temp="";
			if(cnt==1 && type==3){ //padded short
				val=jdio.readmotorolashort(is); off+=2;
				int zero=jdio.readmotorolashort(is); off+=2;
			} else if(type==4){ //int
				val=jdio.readmotorolaint(is); off+=4;
			} else if(type==2){ //string
				int stroff=jdio.readmotorolaint(is); off+=4;
				temp=readOffsetString(stroff,cnt);
			} else if (type==5) { //rational data type (numerator and denominator ints at offset)
				int numoff=jdio.readintelint(is); off+=4;
				val=readOffsetMotorolaRational(off,numoff);
			} else {
				val=jdio.readmotorolaint(is); off+=4;
			}
			List<String> entry=new ArrayList<String>();
			entry.add(""+tag); entry.add(""+type); entry.add(""+cnt);
			if(type==2) entry.add(temp);
			else entry.add(""+val);
			Object tagval=new Integer((int)val);
			if(type==2) tagval=temp;
			entry.add(getTagDesc(tag,tagval));
			return entry;
		/*} catch(IOException e){
			e.printStackTrace();
			return null;
		}*/
	}
	
	/****************
	 * here we read a string at some offset from the current position and then return
	 * @param curroff
	 * @param stroff
	 * @param strlen
	 * @return
	 */
	public String readOffsetString(int stroff,int strlen){
		try{
			jdataio jdio=new jdataio();
			jdio.skipstreambytes(is,stroff-off);
			String temp=jdio.readstring(is,strlen);
			//is.seek(curroff); //hopefully this takes us back to our original position
			is=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(is,off);
			return temp;
		} catch(IOException e){
			e.printStackTrace();
			return null;
		}
	}
	
	/****************
	 * here we read a rational (numerator and denominator ints) at some offset from the current position and then return
	 * @param curroff
	 * @param numoff
	 * @return
	 */
	public float readOffsetIntelRational(int numoff){
		try{
			jdataio jdio=new jdataio();
			jdio.skipstreambytes(is,numoff-off);
			int num=jdio.readintelint(is);
			int den=jdio.readintelint(is);
			//is.seek(curroff); //hopefully this takes us back to our original position
			is=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(is,off);
			return (float)num/(float)den;
		} catch(IOException e){
			e.printStackTrace();
			return Float.NaN;
		}
	}
	
	/****************
	 * here we read a rational (numerator and denominator ints) at some offset from the current position and then return
	 * @param curroff
	 * @param numoff
	 * @return
	 */
	public float readOffsetMotorolaRational(int curroff,int numoff){
		try{
			jdataio jdio=new jdataio();
			jdio.skipstreambytes(is,numoff-curroff);
			int num=jdio.readmotorolaint(is);
			int den=jdio.readmotorolaint(is);
			//is.seek(curroff); //hopefully this takes us back to our original position
			is=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(is,curroff);
			return (float)num/(float)den;
		} catch(IOException e){
			e.printStackTrace();
			return Float.NaN;
		}
	}
	
	public String getTagDesc(int tag,Object val){
		if(tag==254) {newsubfiletype=(Integer)val; return "newsubfiletype";}
        if(tag==256) {width=(Integer)val;return "width";}
        if(tag==257) {height=(Integer)val;return "height";}
        if(tag==258) {bitspersample=(Integer)val;return "bitspersample";}
        if(tag==262) {photointerp=(Integer)val;return "photointerp";}
        if(tag==270) {imgdescription=(String)val;return "imgdescription";}
        if(tag==273) {stripoffsets=(Integer)val;return "stripoffsets";}
        if(tag==274) {orientation=(Integer)val;return "orientation";}
        if(tag==277) {samplesperpixel=(Integer)val;return "samplesperpixel";}
        if(tag==283) {xresolution=(Integer)val;return "xresolution";}
        if(tag==284) {yresolution=(Integer)val;return "yresolution";}
        if(tag==296) {resolutionunit=(Integer)val;return "resolutionunit";}
        if(tag==306) {datetime=(Integer)val;return "datetime";}
        if(tag==320) {colormap=(Integer)val;return "colormap";}
        if(tag==339) {sampleformat=(Integer)val;return "sampleformat";}
        if(tag==278) {rowsperstrip=(Integer)val;return "rowsperstrip";}
        if(tag==279) {stripbytecount=(Integer)val;return "stripbytecount";}
        return "unknown";
	}

}
