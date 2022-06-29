package jguis;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import ij.IJ;
import jalgs.algutils;
import jalgs.jdataio;

public class NpyFileReader{
	
	public static Object[] readNpyArray(String path) {
		//this reads a numpy array file
		//the return object array contains the full array (as a single array) and the shape array (as an int array)
		try{
			jdataio jdio=new jdataio();
			InputStream instream=new BufferedInputStream(new FileInputStream(path));
			//first read the magic number
			String magic=jdio.readstring(instream,6);
			//and the major and minor format
			byte[] majorminor=new byte[2];
			jdio.readintelbytefile(instream,2,majorminor);
			//now the header length
			int headlength=jdio.readintelshort(instream)&0xffff;
			//now read the header
			byte[] header=new byte[headlength];
			jdio.readintelbytefile(instream,headlength,header);
			int[][] headvalues=parseHeader(header);
			//now read the array
			if(headvalues==null) {
				IJ.log("error parsing header");
				instream.close();
				return null;
			}
			int arrlength=headvalues[1][0];
			for(int i=1;i<headvalues[1].length;i++) {
				arrlength*=headvalues[1][i];
			}
			Object[] arr=null;
			if(headvalues[0][0]==0) {
				//this is a byte array
				byte[] combined=new byte[arrlength];
				if(!jdio.readintelbytefile(instream,arrlength,combined)) {
					IJ.error("error reading file");
					instream.close();
					return null;
				}
				arr=new Object[] {combined,headvalues[1]};
			}
			if(headvalues[0][0]==1) {
				//this is a short array
				short[] combined=new short[arrlength];
				if(!jdio.readintelshortfile(instream,arrlength,combined)) {
					IJ.error("error reading file");
					instream.close();
					return null;
				}
				arr=new Object[] {combined,headvalues[1]};
			}
			if(headvalues[0][0]==2) {
				//this is a float array
				float[] combined=new float[arrlength];
				if(!jdio.readintelfloatfile(instream,arrlength,combined)) {
					IJ.error("error reading file");
					instream.close();
					return null;
				}
				arr=new Object[] {combined,headvalues[1]};
			}
			if(headvalues[0][0]==3) {
				//this is a double array
				double[] combined=new double[arrlength];
				if(!jdio.readinteldoublefile(instream,arrlength,combined)) {
					IJ.error("error reading file");
					instream.close();
					return null;
				}
				arr=new Object[] {combined,headvalues[1]};
			}
			if(headvalues[0][0]==4) {
				//this is a int array
				int[] combined=new int[arrlength];
				if(!jdio.readintelintfile(instream,arrlength,combined)) {
					IJ.error("error reading file");
					instream.close();
					return null;
				}
				arr=new Object[] {combined,headvalues[1]};
			}
			instream.close();
			return arr;
		} catch(IOException e){
			IJ.log(e.getMessage());
			return null;
		}
	}
	
	public static int[][] parseHeader(byte[] header) {
		//this returns an array reporting the data type and shape of the contained array
		//example for byte: {'descr': '|u1', 'fortran_order': False, 'shape': (100,), }
		//data types are byte: '|u1', short: '<u2', float: '<f4', double: '<f8', int: '<i4'
		//assume we always have non-fortran order
		//shapes are simply comma separated dimensions
		//if one dimensional, end with comma (as in example above)
		String[] typestrings= {"'|u1'","'<u2'","'<f4'","'<f8'","'<i4'"};
		String headstring=new String(header);
		String typestring=headstring.substring(10,15);
		IJ.log("type string = "+typestring);
		int dtype=-1;
		for(int i=0;i<typestrings.length;i++) {
			if(typestring.equals(typestrings[i])) dtype=i;
		}
		if(dtype==-1) return null;
		int shapeloc=headstring.indexOf("shape");
		//the shape id is surrounded by parentheses after the shape key word
		shapeloc=+6;
		int shapestart=headstring.indexOf("(",shapeloc);
		int shapeend=headstring.indexOf(")",shapestart+1);
		String shapestring=headstring.substring(shapestart+1,shapeend);
		IJ.log("shape string = "+shapestring);
		String[] shapes=table_tools.split(shapestring,",");
		int ndims=shapes.length;
		if(shapes[shapes.length-1]=="") ndims--;
		int[][] headvalues=new int[2][];
		headvalues[0]=new int[] {dtype};
		headvalues[1]=new int[ndims];
		for(int i=0;i<ndims;i++) {
			headvalues[1][i]=Integer.parseInt(shapes[i].trim());
		}
		return headvalues;
	}

}
