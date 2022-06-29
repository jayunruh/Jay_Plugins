package jguis;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import ij.IJ;
import jalgs.algutils;
import jalgs.jdataio;

public class NpyFileWriter{

	public static boolean writeArrayAsNpy(Object arr,String path) {
		int dtype=algutils.get_array_type(arr);
		if(dtype<0) return false;
		int[] arrshape=algutils.getArrayShape(arr,dtype);
		if(arrshape==null) return false;
		String fname=(new File(path)).getName();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			//first write the magic string
			//first character is 0x93
			jdataio jdio=new jdataio();
			os.write(new byte[] {(byte)0x93});
			String magic="NUMPY";
			os.write(magic.getBytes());
			//now write the major and minor versions of the format (1 and 0)
			if(!jdio.writebytearray(os,new byte[] {(byte)1,(byte)0})) {
				os.close();
				IJ.log("error writing file");
				return false;
			}
			//now find the header length and write it as a short
			String header=getHeader(arrshape,dtype);
			int headlength=header.length();
			//offset to first data point should be multiple of 64 (magic plus version plus this is 10)
			int padlength=(headlength+10)%64;
			if(padlength!=0) padlength=64-padlength;
			jdio.writeintelshort(os,(short)(headlength+padlength));
			jdio.writebytearray(os,header.getBytes());
			jdio.writebytearray(os,makeSpaces(padlength));
			//finally write the array
			writeArray(os,arr,arrshape,dtype);
			os.close();
		} catch(IOException e){
			IJ.log(e.getMessage());
			return false;
		}
		return true;
	}
	
	/*****************
	 * in this version we take a classic image stack and write as 3D array
	 * @param stack
	 * @param path
	 * @return
	 */
	public static boolean writeArrayAsNpy(Object[] stack,int width,int height,String path) {
		int dtype=algutils.get_array_type(stack[0]);
		if(dtype<0) return false;
		Object shaped=null;
		if(dtype==0) {
			byte[][][] temp=new byte[stack.length][][];
			for(int i=0;i<stack.length;i++) temp[i]=(byte[][])algutils.reshape((byte[])stack[i],new int[] {height,width});
			shaped=temp;
		}
		if(dtype==1) {
			short[][][] temp=new short[stack.length][][];
			for(int i=0;i<stack.length;i++) temp[i]=(short[][])algutils.reshape((short[])stack[i],new int[] {height,width});
			shaped=temp;
		}
		if(dtype==2) {
			float[][][] temp=new float[stack.length][][];
			for(int i=0;i<stack.length;i++) temp[i]=(float[][])algutils.reshape((float[])stack[i],new int[] {height,width});
			shaped=temp;
		}
		if(dtype==3) {
			double[][][] temp=new double[stack.length][][];
			for(int i=0;i<stack.length;i++) temp[i]=(double[][])algutils.reshape((double[])stack[i],new int[] {height,width});
			shaped=temp;
		}
		if(dtype==4) {
			int[][][] temp=new int[stack.length][][];
			for(int i=0;i<stack.length;i++) temp[i]=(int[][])algutils.reshape((int[])stack[i],new int[] {height,width});
			shaped=temp;
		}
		return writeArrayAsNpy(shaped,path);
	}
	
	public static boolean writeArray(OutputStream os,Object arr,int[] arrshape,int dtype) {
		if(dtype==0) {
			if(arrshape.length==1) return (new jdataio()).writebytearray(os,(byte[])arr);
			if(arrshape.length==2) return (new jdataio()).writebytearray(os,(byte[][])arr);
			if(arrshape.length==3) return (new jdataio()).writebytearray(os,(byte[][][])arr);
		}
		else if(dtype==1) {
			if(arrshape.length==1) return (new jdataio()).writeintelshortarray(os,(short[])arr);
			if(arrshape.length==2) return (new jdataio()).writeintelshortarray(os,(short[][])arr);
			if(arrshape.length==3) return (new jdataio()).writeintelshortarray(os,(short[][][])arr);
		}
		else if(dtype==2) {
			if(arrshape.length==1) return (new jdataio()).writeintelfloatarray(os,(float[])arr);
			if(arrshape.length==2) return (new jdataio()).writeintelfloatarray(os,(float[][])arr);
			if(arrshape.length==3) return (new jdataio()).writeintelfloatarray(os,(float[][][])arr);
		}
		else if(dtype==3) {
			if(arrshape.length==1) return (new jdataio()).writeinteldoublearray(os,(double[])arr);
			if(arrshape.length==2) return (new jdataio()).writeinteldoublearray(os,(double[][])arr);
			if(arrshape.length==3) return (new jdataio()).writeinteldoublearray(os,(double[][][])arr);
		}
		else if(dtype==4) {
			if(arrshape.length==1) return (new jdataio()).writeintelintarray(os,(int[])arr);
			if(arrshape.length==2) return (new jdataio()).writeintelintarray(os,(int[][])arr);
			if(arrshape.length==3) return (new jdataio()).writeintelintarray(os,(int[][][])arr);
		}
		return false;
	}
	
	public static byte[] makeSpaces(int nspaces) {
		byte[] spaces=new byte[nspaces];
		for(int i=0;i<nspaces;i++) spaces[i]=(byte)0x20;
		return spaces;
	}
	
	public static String getHeader(int[] shape,int dtype) {
		//this returns a dictionary reporting the shape and data type of the contained array
		//example for byte: {'descr': '|u1', 'fortran_order': False, 'shape': (100,), }
		//data types are byte: '|u1', short: '<u2', float: '<f4', double: '<f8', int: '<i4'
		//assume we always want non-fortran order
		//shapes are simply comma separated dimensions
		//if one dimensional, end with comma (as in example above)
		String header="{'descr': ";
		String[] typestrings= {"'|u1'","'<u2'","'<f4'","'<f8'","'<i4'"};
		header+=typestrings[dtype];
		header+=", 'fortran_order': False, 'shape': ";
		String shapestring="("+shape[0];
		for(int i=1;i<shape.length;i++) {
			shapestring+=", "+shape[i];
		}
		if(shape.length<2) shapestring+=",)";
		else shapestring+=")";
		header+=shapestring+", }";
		return header;
	}

}
