/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/* DeltaVision File Opener by Fabrice Cordelieres  & Sebastien HUART
 Institut Curie
 UMR146
 ORSAY
 FRANCE
 */

public class DVFile{

	private static final byte IW_BYTE=0;
	private static final byte IW_SHORT=1;
	private static final byte IW_FLOAT=2;
	private static final byte IW_COMPLEX_SHORT=3;
	private static final byte IW_COMPLEX=4;
	private static final byte IW_USHORT=5;

	/*
	 * useless private static final byte ZTW_SEQ=0; private static final byte
	 * WZT_SEQ=1; private static final byte ZWT_SEQ=2;
	 */
	private static final String[] dvPixelTypes={"BYTE","SHORT","FLOAT","COMPLEX_SHORT","COMPLEX","SHORT","USHORT","LONG"};
	private static final String[] ImageSequenceTag={"ZTW","WZT","ZWT"};
	private static final String[] dvExtendedHeaderFloatsDesc={"PHOTOSENSOR_MEASUREMENT","ELAPSED_TIME","STAGE_X","STAGE_Y","STAGE_Z","MIN_INTENSITY","MAX_INTENSITY","MEAN_INTENSITY","EXPOSURE_TIME",
			"NEUTRAL_DENSITY","EXCITATION_WL","EMISSION_WL","INTENSITY_SCALING","ENERGY_CONVERSION"};

	private static short bigendiantag=(short)0xc0a0; // for file endianness
	// test
	private static final byte[] dvtypesizes={1,2,4,4,8,2,2,4};

	private String status;
	private String error;
	private byte[] dvHeader;
	private byte[] dvExtendedHeader; // raw extended Header
	private int dvExtendedHeaderSize;
	private int dvExtendedHeaderNumFloats;
	private int dvExtendedHeaderNumInts;
	private float[][] dvExtendedHeaderFloats;
	private int[][] dvExtendedHeaderInts;
	private int dvImageDataOffset;
	private int dvNumOfImages;
	private int NWL; // maximum number of wavelengths
	private int[] wavelength; // wavelengths in nanometers
	private float[] dvmin; // minimum intensities (1/wavelength)
	private float[] dvmax; // maximum intensities (idem)
	private int imagesequence; // img sequence tag
	private int dvImageWidth;
	private int dvImageHeight;
	private float dvPixelWidth;
	private float dvPixelHeight;
	private float dvPixelDepth;
	private int dvPixelType;
	private boolean fileisbigendian;
	private String dvMetaData;
	private float dvZorigin;
	private float dvYorigin;
	private float dvXorigin;
	private short dvLensID;
	private short dvNumTimes;

	/* 4 static methods for byte swapping and integer conversions */

	// swap signed byte array and convert to int

	public static final int sBtosI(boolean BIGENDIAN,byte b1,byte b2,byte b3,byte b4){
		if(!BIGENDIAN)
			return ((b4&0xff)<<24)|((b3&0xff)<<16)|((b2&0xff)<<8)|(b1&0xff);
		else
			return ((b1&0xff)<<24)|((b2&0xff)<<16)|((b3&0xff)<<8)|(b4&0xff);
	}

	public static final int sBtosI(boolean BIGENDIAN,byte[] barray,int offset){
		return sBtosI(BIGENDIAN,barray[offset],barray[offset+1],barray[offset+2],barray[offset+3]);
	}

	// same stuff for short ints

	public static final short sBtosS(boolean BIGENDIAN,byte b1,byte b2){
		if(!BIGENDIAN)
			return (short)(((b2&0xff)<<8)|(b1&0xff));
		else
			return (short)(((b1&0xff)<<8)|(b2&0xff));
	}

	public static final short sBtosS(boolean BIGENDIAN,byte[] barray,int offset){
		return sBtosS(BIGENDIAN,barray[offset],barray[offset+1]);
	}

	public String getStatus(){
		return status;
	}

	void setStatus(String txt){
		status=txt;
	}

	public String getError(){
		return error;
	}

	void setError(String message,boolean append){
		if(append){
			error+="\n"+message;
		}else
			error=message;
	}

	void setError(String message){
		error=message;
	}

	public String getMetaDataString(){
		return dvMetaData;
	}

	public int getPixelType(){
		return dvPixelType;
	}

	public String getPixelTypeString(){
		return dvPixelTypes[dvPixelType];
	}

	public float getPixelWidth(){
		return dvPixelWidth;
	}

	public float getPixelHeight(){
		return dvPixelHeight;
	}

	public float getPixelDepth(){
		return dvPixelDepth;
	}

	public int getImageWidth(){
		return dvImageWidth;
	}

	public int getImageHeight(){
		return dvImageHeight;
	}

	public int getNumOfImages(){
		return dvNumOfImages;
	}

	public int getNChannels(){
		return NWL;
	}

	public int getNFrames(){
		return dvNumTimes&0xffff;
	}

	public int getNSlices(){
		return (dvNumOfImages/NWL)/getNFrames();
	}

	public boolean intelByteOrder(){
		return(!fileisbigendian);
	}

	public int getImageDataOffset(){
		return dvImageDataOffset;
	}

	/*
	 * this method populates an array of float with some useful infos read in
	 * the extended header: stage xyz position, min, max, exposure time.... for
	 * each acquisition
	 */
	void parseExtendedHeader(){
		dvExtendedHeaderNumInts=sBtosS(fileisbigendian,dvHeader,128);
		dvExtendedHeaderNumFloats=sBtosS(fileisbigendian,dvHeader,130);
		int sectionSize=4*(dvExtendedHeaderNumFloats+dvExtendedHeaderNumInts);
		int sections=0;
		if(sectionSize>0) sections=dvExtendedHeaderSize/sectionSize;
		if(sections<dvNumOfImages){
			setError("Bad number of sections in Extended Header, will not parse...");
			return;
		}
		sections=dvNumOfImages;// last sections (>dvNumOfImages are zeroes)
		// only 14 floats/section are useful, ints and other floats = 0
		int parsed=14;
		dvExtendedHeaderFloats=new float[sections][parsed];

		for(int i=0;i<sections;i++){
			for(int k=0;k<parsed;k++)
				dvExtendedHeaderFloats[i][k]=Float.intBitsToFloat(sBtosI(fileisbigendian,dvExtendedHeader,i*sectionSize+(k+dvExtendedHeaderNumInts)*4));
		}
	}

	/* extended header for humans */
	String getExtendedHeaderInfos(){
		String xtdinfos="";
		StringBuffer sb=new StringBuffer(1024);

		if(dvExtendedHeaderFloats!=null){
			int sections=dvExtendedHeaderFloats.length;
    		setStatus("getting ExtendedHeaderInfos:\nsections:"+sections);
    		for(int i=0;i<sections;i++){
    			for(int j=0;j<dvExtendedHeaderFloats[0].length;j++)
    				sb.append("\nSection["+i+"]_"+dvExtendedHeaderFloatsDesc[j]+":"+dvExtendedHeaderFloats[i][j]);
    		}
		}
		xtdinfos=sb.toString();
		return xtdinfos;
	}

	public DVFile(String directory,String name) throws IOException{

		// Read Header
		InputStream is;
		dvHeader=new byte[1024];
		is=new FileInputStream(directory+name);
		is.read(dvHeader,0,1024);
		dvExtendedHeaderSize=sBtosI(fileisbigendian,dvHeader,92);
		if(dvExtendedHeaderSize>0){
			dvExtendedHeader=new byte[dvExtendedHeaderSize];
			is.read(dvExtendedHeader,0,dvExtendedHeaderSize);

		}
		is.close();
		// check endianness first...
		// Java's natural byte order is "Big Endian"
		// Intel based machines are "Little Endian"
		// MIPS(SGI),SPARC(Sun),PowerPC(Apple,IBM) are "Big Endian"
		// Deltavision format can be either Big or Little endian depending on
		// the acquisition system native order:
		short filestamp=(short)(((dvHeader[96]&0xff)<<8)|(dvHeader[97]&0xff));
		fileisbigendian=(filestamp==bigendiantag)?true:false;
		// Reading size infos inside of the header

		dvImageWidth=sBtosI(fileisbigendian,dvHeader,0);
		dvImageHeight=sBtosI(fileisbigendian,dvHeader,4);
		dvNumOfImages=sBtosI(fileisbigendian,dvHeader,8);
		dvPixelWidth=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,40));
		dvPixelHeight=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,44));
		dvPixelDepth=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,48));
		dvImageDataOffset=1024+(dvExtendedHeaderSize);
		dvPixelType=sBtosI(fileisbigendian,dvHeader,12);
		if((dvPixelType<0)||(dvPixelType>7)) // should never happen
		{
			throw new IOException("Unsupported pixel Type"+dvPixelType);
		}
		parseExtendedHeader();

		dvLensID=sBtosS(fileisbigendian,dvHeader,162);
		dvNumTimes=sBtosS(fileisbigendian,dvHeader,180);
		dvZorigin=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,208));
		dvXorigin=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,212));
		dvYorigin=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,216));
		// number of wavelengths
		NWL=sBtosS(fileisbigendian,dvHeader,196);
		dvMetaData="NumOfWavelengthes:"+NWL+"\n";

		wavelength=new int[NWL];
		for(int i=0;i<NWL;i++) // read wavelengths values
		{
			wavelength[i]=sBtosS(fileisbigendian,dvHeader[198+2*i],dvHeader[198+2*i+1]);
			dvMetaData+="W"+i+":"+wavelength[i]+"\n";
		}
		dvmin=new float[NWL];
		dvmax=new float[NWL];
		// min and max intensity for each wavelength: (not stored consecutively
		// in dv file!)
		dvmin[0]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,76));
		dvmax[0]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,80));
		if(NWL>1){
			int tmp=NWL;
			tmp=(tmp<5)?tmp:4;
			for(int i=1;i<tmp;i++){
				dvmin[i]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,128+8*i));
				dvmax[i]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,132+8*i));
			}
			if(NWL==5){
				dvmin[4]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,172));
				dvmax[4]=java.lang.Float.intBitsToFloat(sBtosI(fileisbigendian,dvHeader,176));
			}
		}
		for(int i=0;i<NWL;i++){
			dvMetaData+="W"+i+"minIntensity:"+dvmin[i]+"\nW"+i+"maxIntensity:"+dvmax[i]+"\n";
		}
		imagesequence=sBtosS(fileisbigendian,dvHeader,182);
		dvMetaData+="ImageSequence:"+ImageSequenceTag[imagesequence]+"\nXorigin:"+dvXorigin+"\nYorigin:"+dvYorigin+"\nZorigin:"+dvZorigin+"\nLensID:"+dvLensID+"\nNumTimes:"+dvNumTimes+"\nPixelType:"
				+getPixelTypeString()+"\nExtendedHeaderInfos:"+getExtendedHeaderInfos();

	}
}
