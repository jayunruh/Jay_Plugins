/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;

import com.xuggle.xuggler.ICodec;
import com.xuggle.xuggler.IContainer;
import com.xuggle.xuggler.IPacket;
import com.xuggle.xuggler.IPixelFormat;
import com.xuggle.xuggler.IStream;
import com.xuggle.xuggler.IStreamCoder;
import com.xuggle.xuggler.IVideoPicture;
import com.xuggle.xuggler.IVideoResampler;
import com.xuggle.xuggler.video.ConverterFactory;
import com.xuggle.xuggler.video.IConverter;

public class Xuggler_file_reader{

	private IContainer container;
	private int videoStreamId;
	private IStreamCoder videoCoder;
	private IVideoResampler resampler;
	private IConverter converter;
	public int width,height,npackets;
	public double duration;
	private IPixelFormat.Type pformat;
	private IPixelFormat.Type rformat;
	private long[] packetstamps;
	private double ptimebase;
	public int nframes;
	public double frameInterval;
	public String name;
	public File fimg;
	public ImagePlus firstframe;

	public Xuggler_file_reader(File fimg){
		this.fimg=fimg;
		name=fimg.getName();
		setup_xuggler(fimg);
	}

	public Xuggler_file_reader(String directory,String filename){
		fimg=new File(directory,filename);
		name=fimg.getName();
		setup_xuggler(fimg);
	}

	public ImagePlus open(){
		return subopen(0,0,1,false,null,1);
	}

	public ImagePlus subopen(int startframe,int endframe,int frameinterval,boolean gray,Rectangle crop,int bin){
		ImageStack stack=null;
		if(crop!=null){
			if(bin>1){
				stack=new ImageStack(crop.width/bin,crop.height/bin);
			}else{
				stack=new ImageStack(crop.width,crop.height);
			}
		}else{
			if(bin>1){
				stack=new ImageStack(width/bin,height/bin);
			}else{
				stack=new ImageStack(width,height);
			}
		}
		IPacket packet=IPacket.make();
		int piccounter=0;
		int packetcounter=0;
		IJ.showStatus("decoding movie");
		fileloop: while(container.readNextPacket(packet)>=0){
			if(packet.getStreamIndex()==videoStreamId){
				IVideoPicture picture=IVideoPicture.make(pformat,width,height);
				int offset=0;
				while(offset<packet.getSize()){
					int bytesDecoded=videoCoder.decodeVideo(picture,packet,offset);
					if(bytesDecoded<0)
						throw new RuntimeException("could not decode video");
					offset+=bytesDecoded;
					if(picture.isComplete()){
						if(endframe>0&&piccounter>endframe){
							break fileloop;
						}
						if(piccounter>=startframe){
							if((piccounter-startframe)%frameinterval==0){
								IVideoPicture newPict=picture;
								if(resampler!=null){
									newPict=IVideoPicture.make(rformat,width,height);
									if(resampler.resample(newPict,picture)<0)
										throw new RuntimeException("could not resample video");
								}
								if(newPict.getPixelType()!=IPixelFormat.Type.BGR24){
									throw new RuntimeException("could not decode video");
								}
								BufferedImage javaImage=converter.toImage(newPict);
								ImageProcessor ip=new ColorProcessor(javaImage);
								if(crop!=null)
									ip=jutils.cropProcessor(ip,crop);
								if(gray)
									ip=ip.convertToByte(true);
								if(bin>1)
									ip=jutils.binProcessor(ip,bin,bin);
								stack.addSlice("",ip);
							}
						}
						piccounter++;
					}
				}
				packetcounter++;
				IJ.showProgress(packetcounter,npackets);
			}
		}
		reset_container();
		IJ.showStatus("");
		return new ImagePlus(name,stack);
	}

	public void setup_xuggler(File fimg){
		IJ.showStatus("Opening Movie");
		container=IContainer.make();
		if(container.open(fimg.getPath(),IContainer.Type.READ,null)<0)
			throw new IllegalArgumentException("could not open");
		int numStreams=container.getNumStreams();
		videoStreamId=-1;
		videoCoder=null;
		for(int i=0;i<numStreams;i++){
			IStream stream=container.getStream(i);
			IStreamCoder coder=stream.getStreamCoder();
			if(coder.getCodecType()==ICodec.Type.CODEC_TYPE_VIDEO){
				videoStreamId=i;
				videoCoder=coder;
				break;
			}
		}
		if(videoStreamId==-1){
			throw new RuntimeException("could not find video");
		}
		if(videoCoder.open()<0){
			throw new RuntimeException("got an exception opening the coder");
		}
		width=videoCoder.getWidth();
		height=videoCoder.getHeight();
		pformat=videoCoder.getPixelType();
		resampler=null;
		if(videoCoder.getPixelType()!=IPixelFormat.Type.BGR24){
			resampler=IVideoResampler.make(width,height,IPixelFormat.Type.BGR24,width,height,videoCoder.getPixelType());
			if(resampler==null)
				throw new RuntimeException("could not get resampler");
		}
		rformat=resampler.getOutputPixelFormat();
		converter=ConverterFactory.createConverter(ConverterFactory.XUGGLER_BGR_24,IPixelFormat.Type.BGR24,width,height);
		IJ.showStatus("Getting Packet Info");
		setPacketInfo();
		IJ.showStatus("Estimating Frame Interval");
		setFrameInfo();
		IJ.showStatus("Getting Preview Frame");
		getFirstFrame();
	}

	private void setFrameInfo(){
		if(duration<5.0){
			// in this case, we read the entire data set
			double[] temp=get_frame_interval(true,0.0,0.0);
			frameInterval=temp[0];
			nframes=(int)temp[1];
		}else{
			// here just try reading the first 5 seconds
			if(npackets<300){
				// this is unlikely, but if it is the case, we will just read
				// the first five seconds contiguously
				double[] temp=get_frame_interval(false,0.0,5.0);
				frameInterval=temp[0];
				nframes=(int)(duration/frameInterval);
			}else{
				// this is the typical scenario--movie longer than 10 seconds
				// with more than 600 packets
				// (usually we have around 1 packet per frame)
				// here we sample the video at regular intervals to estimate the
				// avg frame interval in 1 sec bursts
				// be careful to not read past the end
				double interval=duration/5.0;
				double sum=0.0;
				for(int i=0;i<5;i++){
					double time=interval*i;
					double[] temp=get_frame_interval(false,time,2.0);
					sum+=temp[0];
				}
				frameInterval=sum/5.0;
				nframes=(int)(duration/frameInterval);
			}
		}
		reset_container();
	}

	private double[] get_frame_interval(boolean alltime,double starttime,double length){
		IPacket packet=IPacket.make();
		if(!alltime){
			seek_packet(starttime,packet,-2);
		}
		boolean first=true;
		long thisstarttime=0L;
		long thisendtime=0L;
		long targettime=0L;
		int npictures=0;
		temploop: while(container.readNextPacket(packet)>=0){
			if(packet.getStreamIndex()==videoStreamId){
				IVideoPicture picture=IVideoPicture.make(pformat,width,height);
				int offset=0;
				while(offset<packet.getSize()){
					int bytesDecoded=videoCoder.decodeVideo(picture,packet,offset);
					if(bytesDecoded<0)
						throw new RuntimeException("could not decode video");
					offset+=bytesDecoded;
					if(picture.isComplete()){
						if(first){
							thisstarttime=picture.getTimeStamp();
							targettime=thisstarttime+(int)(length*1000000.0);
							first=false;
						}
						thisendtime=picture.getTimeStamp();
						npictures++;
						if(!alltime){
							if(thisendtime>=targettime){
								// the time is up
								break temploop;
							}
						}
					}
				}
			}
		}
		// first get the time interval in seconds
		double interval=(thisendtime-thisstarttime)/1000000.0;
		// now divide by the number of measured pictures
		interval/=npictures;
		double[] temp2={interval,npictures};
		return temp2;
	}

	private void getFirstFrame(){
		IPacket packet=IPacket.make();
		temploop: while(container.readNextPacket(packet)>=0){
			if(packet.getStreamIndex()==videoStreamId){
				IVideoPicture picture=IVideoPicture.make(pformat,width,height);
				int offset=0;
				while(offset<packet.getSize()){
					int bytesDecoded=videoCoder.decodeVideo(picture,packet,offset);
					if(bytesDecoded<0)
						throw new RuntimeException("could not decode video");
					offset+=bytesDecoded;
					if(picture.isComplete()){
						IVideoPicture newPict=picture;
						if(resampler!=null){
							newPict=IVideoPicture.make(rformat,width,height);
							if(resampler.resample(newPict,picture)<0)
								throw new RuntimeException("could not resample video");
						}
						if(newPict.getPixelType()!=IPixelFormat.Type.BGR24){
							throw new RuntimeException("could not decode video");
						}
						BufferedImage javaImage=converter.toImage(newPict);
						ImageProcessor ip=new ColorProcessor(javaImage);
						firstframe=new ImagePlus(name+" preview",ip);
						break temploop;
					}
				}
			}
		}
		reset_container();
	}

	private void setPacketInfo(){
		IPacket packet=IPacket.make();
		npackets=0;
		boolean first=true;
		long[] temp=new long[1000000];
		while(container.readNextPacket(packet)>=0){
			if(packet.getStreamIndex()==videoStreamId){
				temp[npackets]=packet.getTimeStamp();
				if(first){
					ptimebase=packet.getTimeBase().getDouble();
				}
				npackets++;
			}
		}
		packetstamps=new long[npackets];
		System.arraycopy(temp,0,packetstamps,0,npackets);
		duration=ptimebase*(packetstamps[npackets-1]-packetstamps[0]);
		reset_container();
	}

	private boolean seek_packet(double time,IPacket packet){
		// here we find the packet just before the proposed time
		return seek_packet(time,packet,-1);
	}

	private boolean seek_packet(double time,IPacket packet,int offset){
		// here we find the packet offset packets before the proposed time
		// seeking is always forwards from the beginning
		reset_container();
		if(time>duration){
			return seek_packet(npackets,packet);
		}else{
			for(int i=0;i<npackets;i++){
				double currtime=ptimebase*(packetstamps[i]-packetstamps[0]);
				if(currtime>=time){
					// System.out.println("selected packet = "+(i+offset));
					return seek_packet(i+offset,packet);
				}
			}
		}
		return false;
	}

	private boolean seek_packet(int n,IPacket packet){
		// note that there is no checking to see if the container has been reset
		// this is to allow for packet to packet seeking (could be dangerous)
		if(n<=0){
			return true;
		}
		int counter=0;
		temploop: while(container.readNextPacket(packet)>=0){
			if(packet.getStreamIndex()==videoStreamId){
				counter++;
				if(counter>=n){
					break temploop;
				}
			}
		}
		if(counter>=n)
			return true;
		else
			return false;
	}

	public boolean reset_container(){
		if(container.seekKeyFrame(videoStreamId,-1,0)>=0)
			return true;
		else
			return false;
	}

	public void close_xuggler(){
		if(videoCoder!=null){
			videoCoder.close();
			videoCoder=null;
		}
		if(container!=null){
			container.close();
			container=null;
		}
		videoStreamId=-1;
		IJ.showStatus("");
	}

}
