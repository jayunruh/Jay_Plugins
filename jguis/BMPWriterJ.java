package jguis;

import ij.IJ;
import ij.process.ColorProcessor;
import jalgs.jdataio;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class BMPWriterJ{
	
	public static void writeBMP32(int[] pix,int width,int height,String path){
		jdataio jdio=new jdataio();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			//the file header (14 bytes)
			byte[] type={(byte)'B',(byte)'M'};
			jdio.writebytearray(os, type); //type
			jdio.writeintelint(os,0); //size
			jdio.writeintelint(os,0); //reserve
			jdio.writeintelint(os,54); //offset
			//the info header (40 bytes)
			jdio.writeintelint(os,40); //headersize
			jdio.writeintelint(os,width); //width
			jdio.writeintelint(os,height); //height
			jdio.writeintelshort(os,(short)1); //planes
			jdio.writeintelshort(os,(short)32); //bits
			jdio.writeintelint(os,0); //compression
			jdio.writeintelint(os,0); //size
			jdio.writeintelint(os,0); //xres
			jdio.writeintelint(os,0); //yres
			jdio.writeintelint(os,0); //colorsused
			jdio.writeintelint(os,0); //colorsimportant
			int[] pixcopy=new int[width*height];
			int counter=0;
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					pixcopy[counter]=pix[j+(height-i-1)*width];
					counter++;
				}
			}
			jdio.writeintelintarray(os,pixcopy);
			os.close();
		} catch(IOException e){
			IJ.log(jdio.getExceptionTrace(e));
		}
	}
	
	public static void writeBMP24(ColorProcessor cp,String path){
		int[] pix=(int[])cp.getPixels();
		jdataio jdio=new jdataio();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			//the file header (14 bytes)
			byte[] type={(byte)'B',(byte)'M'};
			jdio.writebytearray(os, type); //type
			jdio.writeintelint(os,0); //size
			jdio.writeintelint(os,0); //reserve
			jdio.writeintelint(os,54); //offset
			//the info header (40 bytes)
			jdio.writeintelint(os,40); //headersize
			jdio.writeintelint(os,cp.getWidth()); //width
			jdio.writeintelint(os,cp.getHeight()); //height
			jdio.writeintelshort(os,(short)1); //planes
			jdio.writeintelshort(os,(short)24); //bits
			jdio.writeintelint(os,0); //compression
			jdio.writeintelint(os,0); //size
			jdio.writeintelint(os,0); //xres
			jdio.writeintelint(os,0); //yres
			jdio.writeintelint(os,0); //colorsused
			jdio.writeintelint(os,0); //colorsimportant
			int width=cp.getWidth(); int height=cp.getHeight();
			byte[] data=new byte[width*height*3];
			int counter=0;
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
    				int[] rgb=jutils.intval2rgb(pix[j+i*width]);
    				data[counter]=(byte)rgb[2]; counter++;
    				data[counter]=(byte)rgb[1]; counter++;
    				data[counter]=(byte)rgb[0]; counter++;
				}
			}
			jdio.writebytearray(os,data);
			os.close();
		} catch(IOException e){
			IJ.log(jdio.getExceptionTrace(e));
		}
	}

}
