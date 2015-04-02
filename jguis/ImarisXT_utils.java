package jguis;
import jalgs.jdataio;
import ij.IJ;
import Imaris.Error;
import Imaris.IDataSetPrx;
import Imaris.ISpotsPrx;
import Imaris.tType;

import com.bitplane.xt.*;

public class ImarisXT_utils{
	
	public static float[][] getSpotsTraj(boolean pixunits){
		//ImarisLib vImarisLib=new ImarisLib();
		BPImarisLib lib=new BPImarisLib();
		Imaris.IApplicationPrx app=lib.GetApplication(0);
		if(app==null){IJ.log("Imaris is not open"); return null;}
		try{
			ISpotsPrx spots=app.GetFactory().ToSpots(app.GetSurpassSelection());
			float[][] coords=spots.GetPositionsXYZ(); //units are micron
			//IJ.log(table_tools.print_float_array(coords));

			IDataSetPrx dataset=app.GetDataSet();
			int chans=dataset.GetSizeC();
			int width=dataset.GetSizeX();
			int height=dataset.GetSizeY();
			int slices=dataset.GetSizeZ();
			float xoff=dataset.GetExtendMinX();
			float yoff=dataset.GetExtendMinY();
			float zoff=dataset.GetExtendMinZ();
			float psize=(dataset.GetExtendMaxX()-xoff)/(float)width;
			float pdepth=(dataset.GetExtendMaxZ()-zoff)/(float)slices;
			float[][] coords2=new float[3][coords.length];
			for(int i=0;i<coords.length;i++){
				coords2[0][i]=coords[i][0]-xoff;
				coords2[1][i]=coords[i][1]-yoff;
				coords2[2][i]=coords[i][2]-zoff;
			}
			if(pixunits){
				for(int i=0;i<coords.length;i++){
					coords2[0][i]/=psize;
					coords2[1][i]/=psize;
					coords2[2][i]/=psize;
				}
			}
			return coords2;
		}catch(Error e){
			IJ.log((new jdataio()).getExceptionTrace(e));
			//e.printStackTrace();
			return null;
		}
	}
	
	public static boolean setSpotsTraj(boolean pixunits,float[][] coords){
		//ImarisLib vImarisLib=new ImarisLib();
		BPImarisLib lib=new BPImarisLib();
		Imaris.IApplicationPrx app=lib.GetApplication(0);
		if(app==null){IJ.log("Imaris is not open"); return false;}
		try{
			IDataSetPrx dataset=app.GetDataSet();
			//int chans=dataset.GetSizeC();
			int width=dataset.GetSizeX();
			//int height=dataset.GetSizeY();
			//int slices=dataset.GetSizeZ();
			float xoff=dataset.GetExtendMinX();
			float yoff=dataset.GetExtendMinY();
			float zoff=dataset.GetExtendMinZ();
			float psize=(dataset.GetExtendMaxX()-xoff)/(float)width;
			//float pdepth=(dataset.GetExtendMaxZ()-zoff)/(float)slices;
			float[][] coords2=new float[coords[0].length][3];
			float[] radii=new float[coords[0].length];
			int[] indices=new int[coords[0].length];
			for(int i=0;i<coords[0].length;i++){
				if(pixunits){
    				coords2[i][0]=coords[0][i]*psize+xoff;
    				coords2[i][1]=coords[1][i]*psize+yoff;
    				coords2[i][2]=coords[2][i]*psize+zoff;
				} else {
    				coords2[i][0]=coords[0][i]+xoff;
    				coords2[i][1]=coords[1][i]+yoff;
    				coords2[i][2]=coords[2][i]+zoff;
				}
				radii[i]=2.0f;
				indices[i]=i;
			}
			ISpotsPrx spots=app.GetFactory().CreateSpots();
			spots.Set(coords2,indices,radii);
			return true;
		}catch(Error e){
			IJ.log((new jdataio()).getExceptionTrace(e));
			//e.printStackTrace();
			return false;
		}
	}
	
	public static Object[] getDataSet(){
		//ImarisLib vImarisLib=new ImarisLib();
		BPImarisLib lib=new BPImarisLib();
		Imaris.IApplicationPrx app=lib.GetApplication(0);
		if(app==null){IJ.log("Imaris is not open"); return null;}
		try{
			IDataSetPrx dataset=app.GetDataSet();
			int chans=dataset.GetSizeC();
			int frames=dataset.GetSizeT();
			int slices=dataset.GetSizeZ();
			int width=dataset.GetSizeX();
			int height=dataset.GetSizeY();
			float xoff=dataset.GetExtendMinX();
			float yoff=dataset.GetExtendMinY();
			float zoff=dataset.GetExtendMinZ();
			float psize=(dataset.GetExtendMaxX()-xoff)/(float)width;
			float pdepth=(dataset.GetExtendMaxZ()-zoff)/(float)slices;
			Object[] stack=new Object[frames*slices*chans];
			tType type=dataset.GetType();
			int typeindex=0; //bytes
			if(type.compareTo(tType.eTypeUInt16)==0) typeindex=1;
			if(type.compareTo(tType.eTypeFloat)==0) typeindex=2;
			int counter=0;
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<chans;k++){
						if(typeindex==0) stack[counter]=dataset.GetDataSubVolumeAs1DArrayBytes(0,0,j,k,i,width,height,1);
						if(typeindex==1) stack[counter]=dataset.GetDataSubVolumeAs1DArrayShorts(0,0,j,k,i,width,height,1);
						if(typeindex==2) stack[counter]=dataset.GetDataSubVolumeAs1DArrayFloats(0,0,j,k,i,width,height,1);
						//need to convert these from 2D to 1D
						/*if(typeindex==0) stack[counter]=convertto1D(dataset.GetDataSliceBytes(j,k,i));
						if(typeindex==1) stack[counter]=convertto1D(dataset.GetDataSliceShorts(j,k,i));
						if(typeindex==2) stack[counter]=convertto1D(dataset.GetDataSliceFloats(j,k,i));*/
						counter++;
					}
				}
			}
			return stack;
		}catch(Error e){
			IJ.log((new jdataio()).getExceptionTrace(e));
			//e.printStackTrace();
			return null;
		}
	}
	
	/*public static Object convertto1D(Object image){
		if(image instanceof byte[][]){
			byte[][] temp=(byte[][])image;
			byte[] temp2=new byte[temp.length*temp[0].length];
			int width=temp.length;
			int height=temp[0].length;
			int counter=0;
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					temp2[counter]=temp[j][i];
					counter++;
				}
			}
			return temp2;
		}
		if(image instanceof short[][]){
			short[][] temp=(short[][])image;
			short[] temp2=new short[temp.length*temp[0].length];
			int width=temp.length;
			int height=temp[0].length;
			int counter=0;
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					temp2[counter]=temp[j][i];
					counter++;
				}
			}
			return temp2;
		}
		if(image instanceof float[][]){
			float[][] temp=(float[][])image;
			float[] temp2=new float[temp.length*temp[0].length];
			int width=temp.length;
			int height=temp[0].length;
			int counter=0;
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					temp2[counter]=temp[j][i];
					counter++;
				}
			}
			return temp2;
		}
		return null;
	}*/
	
	public static float[] getDimensions(){
		//output is width, height, channels, slices, frames, psize, and pdepth
		BPImarisLib lib=new BPImarisLib();
		Imaris.IApplicationPrx app=lib.GetApplication(0);
		if(app==null){IJ.log("Imaris is not open"); return null;}
		try{
			IDataSetPrx dataset=app.GetDataSet();
			int chans=dataset.GetSizeC();
			int frames=dataset.GetSizeT();
			int slices=dataset.GetSizeZ();
			int width=dataset.GetSizeX();
			int height=dataset.GetSizeY();
			float xoff=dataset.GetExtendMinX();
			//float yoff=dataset.GetExtendMinY();
			float zoff=dataset.GetExtendMinZ();
			float psize=(dataset.GetExtendMaxX()-xoff)/(float)width;
			float pdepth=(dataset.GetExtendMaxZ()-zoff)/(float)slices;
			return new float[]{(float)width,(float)height,(float)chans,(float)slices,(float)frames,psize,pdepth};
		}catch(Error e){
			IJ.log((new jdataio()).getExceptionTrace(e));
			//e.printStackTrace();
			return null;
		}
		
	}
	
	public static String getSelectionName(){
		BPImarisLib lib=new BPImarisLib();
		Imaris.IApplicationPrx app=lib.GetApplication(0);
		if(app==null){IJ.log("Imaris is not open"); return null;}
		try{
			return app.GetSurpassSelection().GetName();
		}catch(Error e){
			IJ.log((new jdataio()).getExceptionTrace(e));
			//e.printStackTrace();
			return null;
		}
	}

}
