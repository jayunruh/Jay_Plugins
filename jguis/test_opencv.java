package jguis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import jalgs.jseg.binary_processing;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfDouble;
import org.opencv.core.MatOfInt;
import org.opencv.core.MatOfRect;
import org.opencv.core.MatOfPoint;
import org.opencv.core.Point;
import org.opencv.core.Rect;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;
import org.opencv.objdetect.CascadeClassifier;
import org.opencv.objdetect.Objdetect;

public class test_opencv{
	
	public static void main(String[] args){
		System.out.println("native name = "+Core.NATIVE_LIBRARY_NAME);
		ImagePlus imp=IJ.openImage("http://imagej.nih.gov/ij/images/lena.jpg");
		byte[] pix=(byte[])imp.getProcessor().convertToByteProcessor().getPixels();
		Rectangle[] faces=find_faces(pix,imp.getWidth(),imp.getHeight(),1.1,5);
		System.out.println(""+faces[0].x+","+faces[0].y+","+faces[0].width+","+faces[0].height);
		ImagePlus imp2=IJ.openImage("http://imagej.nih.gov/ij/images/blobs.gif");
		ImageProcessor ip=imp2.getProcessor();
		ip.invertLut();
		byte[] pix2=(byte[])ip.convertToByteProcessor().getPixels();
		
		int width=imp2.getWidth(); int height=imp2.getHeight();
		
		int[] contour=track_cells_opencv(pix2,width,height,150,100);
		byte[] filled=create_filled_image(contour,width,height,true);
		int area=0;
		for(int i=0;i<filled.length;i++){
			if(filled[i]!=(byte)0) area++;
		}
		System.out.println("Cell Size = "+area);
	}
	
	public static Rectangle[] find_faces(byte[] image,int width,int height,double scaleFactor,int minNeighbors){
		String classifierpath="c:\\Program Files\\ImageJ\\haarcascade_frontalface_default.xml";
		return find_faces(image,width,height,scaleFactor,minNeighbors,classifierpath);
	}
	
	public static Rectangle[] find_faces(byte[] image,int width,int height,double scaleFactor,int minNeighbors,String classifierpath){
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		Mat gray=Mat.zeros(height,width,CvType.CV_8U);
		gray.put(0,0,image);
		CascadeClassifier facecascade=new CascadeClassifier(classifierpath);
		MatOfRect faces=new MatOfRect();
		//parameters are image, output, scalefactor, minneighbors, flags, minsize, maxsize
		facecascade.detectMultiScale(gray,faces,scaleFactor,minNeighbors,Objdetect.CASCADE_SCALE_IMAGE,new Size(10,10),new Size(width,height));
		if(!faces.empty()){
			Rect[] faces2=faces.toArray();
			Rectangle[] faces3=new Rectangle[faces2.length];
			for(int i=0;i<faces2.length;i++) faces3[i]=new Rectangle(faces2[i].x,faces2[i].y,faces2[i].width,faces2[i].height);
			return faces3;
		} else {return null;}
	}
	
	public static Object[] find_faces2(byte[] image,int width,int height,double scaleFactor,int minNeighbors,String classifierpath){
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		Mat gray=Mat.zeros(height,width,CvType.CV_8U);
		gray.put(0,0,image);
		CascadeClassifier facecascade=new CascadeClassifier(classifierpath);
		MatOfRect faces=new MatOfRect();
		//parameters are image, output, scalefactor, minneighbors, flags, minsize, maxsize
		//facecascade.detectMultiScale(gray,faces,scaleFactor,minNeighbors,Objdetect.CASCADE_SCALE_IMAGE,new Size(10,10),new Size(width,height));
		MatOfInt rejectLevels=new MatOfInt();
		MatOfDouble levelWeights=new MatOfDouble();
		facecascade.detectMultiScale3(gray,faces,rejectLevels,levelWeights,scaleFactor,minNeighbors,Objdetect.CASCADE_SCALE_IMAGE,new Size(10,10),new Size(width,height),true);
		if(!faces.empty()){
			Rect[] faces2=faces.toArray();
			Rectangle[] faces3=new Rectangle[faces2.length];
			for(int i=0;i<faces2.length;i++) faces3[i]=new Rectangle(faces2[i].x,faces2[i].y,faces2[i].width,faces2[i].height);
			return new Object[]{faces3,rejectLevels.toArray(),levelWeights.toArray()};
		} else {return null;}
	}
	
	public static int[] track_cells_opencv(byte[] image,int width,int height,float thresh1,float thresh2){
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		Mat gray=Mat.zeros(height,width,CvType.CV_8U);
		gray.put(0,0,image);
		Mat edges=new Mat();
		Imgproc.Canny(gray,edges,thresh1,thresh2,3,false);
		Imgproc.dilate(edges,edges,new Mat());
		List<MatOfPoint> contours=new ArrayList<MatOfPoint>();
		Imgproc.findContours(edges,contours,new Mat(),Imgproc.RETR_EXTERNAL,Imgproc.CHAIN_APPROX_SIMPLE);
		//find the biggest contour
		double maxarea=0f;
		int maxslice=0;
		if(contours!=null && contours.size()>0){
			for(int i=0;i<contours.size();i++){
				if(!contours.get(i).empty()){
					double tempmax=Imgproc.contourArea(contours.get(i));
					if(tempmax>maxarea){
						maxslice=i;
						maxarea=tempmax;
					}
				}
			}
		} else {return null;}
		Point[] maxcontour=contours.get(maxslice).toArray();
		int[] outdata=new int[2*maxcontour.length];
		for(int i=0;i<maxcontour.length;i++){
			outdata[2*i]=(int)maxcontour[i].x;
			outdata[2*i+1]=(int)maxcontour[i].y;
		}
		return outdata;
	}
	
	public static byte[] create_filled_image(int[] seqpoly,int width,int height,boolean filloutline){
		byte[] retimage=new byte[width*height];
		ByteProcessor bp=new ByteProcessor(width,height,retimage,null);
		int length=seqpoly.length/2;
		if(length>1){
			int[] xpts=new int[length];
			int[] ypts=new int[length];
			for(int i=0;i<length;i++){
				xpts[i]=seqpoly[2*i];
				ypts[i]=seqpoly[2*i+1];
			}
			PolygonRoi temp=new PolygonRoi(xpts,ypts,length,Roi.POLYGON);
			bp.setColor(Color.WHITE);
			bp.draw(temp);
		}
		if(filloutline){
			(new binary_processing(width,height)).fillholes(retimage);
		}
		return retimage;
	}

}
