package jguis;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLNumericArray;

import ij.IJ;
import jalgs.algutils;
import jalgs.jdataio;

public class mat_file_reader{
	
	public static Map<String,MLArray> readContents(String path) {
		try {
			MatFileReader mfr=new MatFileReader(path);
			return mfr.getContent();
		}catch (IOException e) {
			IJ.log(new jdataio().getExceptionTrace(e));
			return null;
		}
	}
	
	public static String[] listObjects(Map<String,MLArray> contents) {
		Set<String> keys=contents.keySet();
		String[] keylist=new String[keys.size()];
		keys.toArray(keylist);
		return keylist;
	}
	
	public static boolean isComplex(Map<String,MLArray> contents,String arrname) {
		return contents.get(arrname).isComplex();
	}
	
	public static Object getRealNumericArray(Map<String,MLArray> contents,String arrname) {
		MLArray mla=contents.get(arrname);
		if(mla instanceof MLNumericArray) {
			MLNumericArray mla2=(MLNumericArray)mla;
			int size=mla2.getSize();
			Number first=mla2.getReal(0);
			int type=algutils.get_number_type(first);
			if(type<0) type=3; //convert everything to double if not defined
			Object arr=algutils.create_array(size,type);
			if(type==0) {
				for(int i=0;i<size;i++) ((byte[])arr)[i]=((MLNumericArray)mla).getReal(i).byteValue();
			}
			if(type==1) {
				for(int i=0;i<size;i++) ((short[])arr)[i]=((MLNumericArray)mla).getReal(i).shortValue();
			}
			if(type==2) {
				for(int i=0;i<size;i++) ((float[])arr)[i]=((MLNumericArray)mla).getReal(i).floatValue();
			}
			if(type==3) {
				for(int i=0;i<size;i++) ((double[])arr)[i]=((MLNumericArray)mla).getReal(i).doubleValue();
			}
			if(type==4) {
				for(int i=0;i<size;i++) ((int[])arr)[i]=((MLNumericArray)mla).getReal(i).intValue();
			}
			return arr;
		} else {
			return null;
		}
	}
	
	public static Object getImaginaryNumericArray(Map<String,MLArray> contents,String arrname) {
		MLArray mla=contents.get(arrname);
		if(mla instanceof MLNumericArray) {
			MLNumericArray mla2=(MLNumericArray)mla;
			if(!mla2.isComplex()) return null;
			int size=mla2.getSize();
			Number first=mla2.getImaginary(0);
			int type=algutils.get_number_type(first);
			if(type<0) type=3; //convert everything to double if not defined
			Object arr=algutils.create_array(size,type);
			if(type==0) {
				for(int i=0;i<size;i++) ((byte[])arr)[i]=((MLNumericArray)mla).getImaginary(i).byteValue();
			}
			if(type==1) {
				for(int i=0;i<size;i++) ((short[])arr)[i]=((MLNumericArray)mla).getImaginary(i).shortValue();
			}
			if(type==2) {
				for(int i=0;i<size;i++) ((float[])arr)[i]=((MLNumericArray)mla).getImaginary(i).floatValue();
			}
			if(type==3) {
				for(int i=0;i<size;i++) ((double[])arr)[i]=((MLNumericArray)mla).getImaginary(i).doubleValue();
			}
			if(type==4) {
				for(int i=0;i<size;i++) ((int[])arr)[i]=((MLNumericArray)mla).getImaginary(i).intValue();
			}
			return arr;
		} else {
			return null;
		}
	}
	
	public static int[] getDimensions(Map<String,MLArray> contents,String arrname) {
		//note that matlab dimension order indixes the last dimension slowest
		return contents.get(arrname).getDimensions();
	}
	
	/*****************
	 * here we expect a 3 dimensional array, the dimension that matches nslices will be the z dimension, complex arrays will be returned in interleaved order
	 * @param contents
	 * @param arrname
	 * @param nslices
	 * @return
	 */
	public static Object[] getImageArray(Map<String,MLArray> contents,String arrname,int nslices) {
		int[] dimensions=getDimensions(contents,arrname);
		MLArray mla=contents.get(arrname);
		Object rawarray=getRealNumericArray(contents,arrname);
		int type=algutils.get_array_type(rawarray);
		Object rawimarray=null;
		int complex=0;
		if(mla.isComplex()) {
			rawimarray=getImaginaryNumericArray(contents,arrname);
			complex=1;
		}
		if(nslices==1) {
			//here we have a single slice
			Object[] stack=new Object[1+complex];
			stack[0]=algutils.create_array(dimensions[0]*dimensions[1],type);
			for(int i=0;i<dimensions[1];i++) {
				Object col=algutils.get_subarray(rawarray,i*dimensions[0],dimensions[1]);
				algutils.set_image_col(col,stack[0],dimensions[1],dimensions[0],i);
			}
			if(rawimarray!=null) {
				stack[1]=algutils.create_array(dimensions[0]*dimensions[1],type);
				for(int i=0;i<dimensions[1];i++) {
					Object col=algutils.get_subarray(rawimarray,i*dimensions[0],dimensions[1]);
					algutils.set_image_col(col,stack[1],dimensions[1],dimensions[0],i);
				}
			}
			return stack;
		}
		if(dimensions[2]==nslices) {
			//this is z,x,y order
			Object[] stack=new Object[nslices+complex*nslices];
			int slicesize=dimensions[0]*dimensions[1];
			if(complex==0) {
    			for(int i=0;i<nslices;i++) {
    				stack[i]=algutils.create_array(slicesize,type);
    				for(int j=0;j<dimensions[1];j++) {//iterate over x
    					Object col=algutils.get_subarray(rawarray,i*slicesize+j*dimensions[0],dimensions[1]);
    					algutils.set_image_col(col,stack[i],dimensions[1],dimensions[0],j);
    					//stack[i]=algutils.get_subarray(rawarray,i*slicesize,slicesize);
    				}
    			}
			} else {
    			for(int i=0;i<nslices;i++) {
    				stack[2*i]=algutils.create_array(slicesize,type);
    				stack[2*i+1]=algutils.create_array(slicesize,type);
    				for(int j=0;j<dimensions[1];j++) {//iterate over x
    					Object col=algutils.get_subarray(rawarray,i*slicesize+j*dimensions[0],dimensions[1]);
    					algutils.set_image_col(col,stack[2*i],dimensions[1],dimensions[0],j);
    					Object imcol=algutils.get_subarray(rawimarray,i*slicesize+j*dimensions[0],dimensions[1]);
    					algutils.set_image_col(imcol,stack[2*i+1],dimensions[1],dimensions[0],j);
    					//stack[2*i]=algutils.get_subarray(rawarray,i*slicesize,slicesize);
    					//stack[2*i+1]=algutils.get_subarray(rawimarray,i*slicesize,slicesize);
    				}
    			}
			}
			return stack;
		}
		if(dimensions[0]==nslices) { 
			int slicesize=dimensions[1]*dimensions[2];
			int tempnslices=nslices+complex*nslices;
			Object[] stack=new Object[tempnslices];
			for(int i=0;i<tempnslices;i++) {
				stack[i]=algutils.create_array(slicesize,type);
			}
			if(complex==0) {
    			for(int i=0;i<slicesize;i++) {
    				Object column=algutils.get_subarray(rawarray,i*nslices,nslices);
    				algutils.set_stack_col(stack,dimensions[1],dimensions[2],i,nslices,column);
    			}
			} else {
    			for(int i=0;i<slicesize;i++) {
    				Object recolumn=algutils.get_subarray(rawarray,i*nslices,nslices);
    				Object imcolumn=algutils.get_subarray(rawimarray,i*nslices,nslices);
    				Object column=interleaveArrays(recolumn,imcolumn);
    				algutils.set_stack_col(stack,dimensions[1],dimensions[2],i,2*nslices,column);
    			}
			}
			return stack;
		}
		return null;
	}
	
	public static Object interleaveArrays(Object arr1,Object arr2) {
		int type=algutils.get_array_type(arr1);
		if(type<0) return null;
		int length=algutils.get_array_length(arr1);
		Object interleaved=algutils.create_array(length,type);
		if(type==0) {
			for(int i=0;i<length;i++) {
				((byte[])interleaved)[2*i]=((byte[])arr1)[i];
				((byte[])interleaved)[2*i+1]=((byte[])arr1)[i];
			}
		}
		if(type==1) {
			for(int i=0;i<length;i++) {
				((short[])interleaved)[2*i]=((short[])arr1)[i];
				((short[])interleaved)[2*i+1]=((short[])arr1)[i];
			}
		}
		if(type==2) {
			for(int i=0;i<length;i++) {
				((float[])interleaved)[2*i]=((float[])arr1)[i];
				((float[])interleaved)[2*i+1]=((float[])arr1)[i];
			}
		}
		if(type==3) {
			for(int i=0;i<length;i++) {
				((double[])interleaved)[2*i]=((double[])arr1)[i];
				((double[])interleaved)[2*i+1]=((double[])arr1)[i];
			}
		}
		if(type==4) {
			for(int i=0;i<length;i++) {
				((int[])interleaved)[2*i]=((int[])arr1)[i];
				((int[])interleaved)[2*i+1]=((int[])arr1)[i];
			}
		}
		return interleaved;
	}

}
