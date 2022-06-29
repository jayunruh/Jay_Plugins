/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class binary_processing{
	public int width,height;
	int maxStackSize=500; // will be increased as needed
	int[] xstack=new int[maxStackSize];
	int[] ystack=new int[maxStackSize];
	int stackSize;

	public binary_processing(int width1,int height1){
		width=width1;
		height=height1;
	}

	public void dilate(byte[] image){
		byte[] newimage=new byte[width*height];
		for(int i=1;i<height-1;i++){
			for(int j=1;j<width-1;j++){
				if(image[j+i*width]!=(byte)0){
					setNeighbors2(newimage,j,i,(byte)255);
				}
			}
		}
		System.arraycopy(newimage,0,image,0,width*height);
	}

	public void make_edge_connected(byte[] image){
		for(int i=1;i<height;i++){
			for(int j=1;j<width;j++){
				byte[] temparr=new byte[4];
				int temp=j-1+(i-1)*width;
				temparr[0]=image[temp];
				temp++;
				temparr[1]=image[temp];
				temp=j-1+i*width;
				temparr[2]=image[temp];
				temp++;
				temparr[3]=image[temp];
				if(temparr[0]==(byte)0&&temparr[1]==(byte)255&&temparr[2]==(byte)255&&temparr[3]==(byte)0){
					image[j+i*width]=(byte)255;
				}else{
					if(temparr[0]==(byte)255&&temparr[1]==(byte)0&&temparr[2]==(byte)0&&temparr[3]==(byte)255){
						image[j+(i-1)*width]=(byte)255;
					}
				}
			}
		}
	}

	public void erode(byte[] image){
		byte[] newimage=new byte[width*height];
		for(int i=0;i<width*height;i++){
			newimage[i]=(byte)255;
		}
		for(int i=1;i<height-1;i++){
			for(int j=1;j<width-1;j++){
				if(image[j+i*width]==(byte)0){
					setNeighbors2(newimage,j,i,(byte)0);
				}
			}
		}
		System.arraycopy(newimage,0,image,0,width*height);
	}

	public void clear_edges(byte[] image){
		int temp=(height-1)*width;
		for(int i=0;i<width;i++){
			image[i]=(byte)0;
			image[temp+i]=(byte)0;
		}
		for(int i=0;i<height;i++){
			image[i*width]=(byte)0;
			image[(i+1)*width-1]=(byte)0;
		}
	}

	static int[] table={0,0,0,1,0,0,1,3,0,0,3,1,1,0,1,3,0,0,0,0,0,0,0,0,2,0,2,0,3,0,3,3,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,3,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,0,0,0,0,3,0,0,0,3,0,2,0,0,0,3,1,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,2,3,1,3,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,0,1,0,0,0,1,0,0,0,0,0,0,0,0,3,3,0,1,0,0,0,0,2,2,0,0,2,0,0,0};

	public void skeletonize(byte[] image){
		// this algorithm and its functions is modified from the ImageJ version
		// in which in turn was based on
		// Zhang and Suen (CACM, March 1984, 236-239)
		// the algorithm came from the orignal BinaryProcessor ImageJ class
		int pass=0;
		int pixelsRemoved;
		clear_edges(image);
		do{
			pixelsRemoved=thin(pass++,table,image);
			pixelsRemoved=thin(pass++,table,image);
		}while(pixelsRemoved>0);
	}

	public int thin(int pass,int[] table,byte[] image){
		int p1,p2,p3,p4,p5,p6,p7,p8,p9;
		int inc=height/25;
		if(inc<1)
			inc=1;
		int bgColor=0;
		byte[] pixels2=image.clone();
		int v,index,code;
		int offset,rowOffset=width;
		int pixelsRemoved=0;
		for(int y=1;y<(height-1);y++){
			offset=1+y*width;
			for(int x=1;x<(width-1);x++){
				p5=pixels2[offset]&0xff;
				v=p5;
				if(v!=bgColor){
					p1=pixels2[offset-rowOffset-1]&0xff;
					p2=pixels2[offset-rowOffset]&0xff;
					p3=pixels2[offset-rowOffset+1]&0xff;
					p4=pixels2[offset-1]&0xff;
					p6=pixels2[offset+1]&0xff;
					p7=pixels2[offset+rowOffset-1]&0xff;
					p8=pixels2[offset+rowOffset]&0xff;
					p9=pixels2[offset+rowOffset+1]&0xff;
					index=0;
					if(p1!=bgColor)
						index|=1;
					if(p2!=bgColor)
						index|=2;
					if(p3!=bgColor)
						index|=4;
					if(p6!=bgColor)
						index|=8;
					if(p9!=bgColor)
						index|=16;
					if(p8!=bgColor)
						index|=32;
					if(p7!=bgColor)
						index|=64;
					if(p4!=bgColor)
						index|=128;
					code=table[index];
					if((pass&1)==1){ // odd pass
						if(code==2||code==3){
							v=bgColor;
							pixelsRemoved++;
						}
					}else{ // even pass
						if(code==1||code==3){
							v=bgColor;
							pixelsRemoved++;
						}
					}
				}
				image[offset++]=(byte)v;
			}
		}
		return pixelsRemoved;
	}

	public byte[] getNeighbors(byte[] image,int x,int y){
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		byte[] temp=new byte[8];
		int temp2=x-1+(y-1)*width;
		temp[0]=image[temp2];
		temp2++;
		temp[1]=image[temp2];
		temp2++;
		temp[2]=image[temp2];
		temp2+=(width-2);
		temp[3]=image[temp2];
		temp2+=2;
		temp[4]=image[temp2];
		temp2+=(width-2);
		temp[5]=image[temp2];
		temp2++;
		temp[6]=image[temp2];
		temp2++;
		temp[7]=image[temp2];
		return temp;
	}

	public void setNeighbors(byte[] image,int x,int y,byte value){
		if(x==0||x>=(width-1)){
			return;
		}
		if(y==0||y>=(height-1)){
			return;
		}
		int temp2=x-1+(y-1)*width;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2+=(width-2);
		image[temp2]=value;
		temp2+=2;
		image[temp2]=value;
		temp2+=(width-2);
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
	}

	public byte[] getNeighbors2(byte[] image,int x,int y){
		if(x==0||x>=(width-1)){
			return null;
		}
		if(y==0||y>=(height-1)){
			return null;
		}
		byte[] temp=new byte[9];
		int temp2=x-1+(y-1)*width;
		temp[0]=image[temp2];
		temp2++;
		temp[1]=image[temp2];
		temp2++;
		temp[2]=image[temp2];
		temp2+=(width-2);
		temp[3]=image[temp2];
		temp2++;
		temp[4]=image[temp2];
		temp2++;
		temp[5]=image[temp2];
		temp2+=(width-2);
		temp[6]=image[temp2];
		temp2++;
		temp[7]=image[temp2];
		temp2++;
		temp[8]=image[temp2];
		return temp;
	}

	public void setNeighbors2(byte[] image,int x,int y,byte value){
		if(x==0||x>=(width-1)){
			return;
		}
		if(y==0||y>=(height-1)){
			return;
		}
		int temp2=x-1+(y-1)*width;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2+=(width-2);
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2+=(width-2);
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
		temp2++;
		image[temp2]=value;
	}

	public boolean floodfill(byte[] data,byte value,int x,int y){
		stackSize=0;
		push(x,y);
		byte color=(byte)0;
		fillLine(data,value,x,x,y);
		data[x+y*width]=(byte)0;
		while(true){
			x=popx();
			if(x==-1)
				return true;
			y=popy();
			if(data[x+width*y]!=color)
				continue;
			int x1=x;
			int x2=x;
			while(x1>=0&&data[x1+width*y]==color)
				x1--; // find start of scan-line
			x1++;
			while(x2<width&&data[x2+width*y]==color)
				x2++; // find end of scan-line
			x2--;
			fillLine(data,value,x1,x2,y); // fill scan-line
			boolean inScanLine=false;
			for(int i=x1;i<=x2;i++){ // find scan-lines above this one
				if(!inScanLine&&y>0&&data[i+width*(y-1)]==color){
					push(i,y-1);
					inScanLine=true;
				}else if(inScanLine&&y>0&&data[i+width*(y-1)]!=color)
					inScanLine=false;
			}
			inScanLine=false;
			for(int i=x1;i<=x2;i++){ // find scan-lines below this one
				if(!inScanLine&&y<height-1&&data[i+width*(y+1)]==color){
					push(i,y+1);
					inScanLine=true;
				}else if(inScanLine&&y<height-1&&data[i+width*(y+1)]!=color)
					inScanLine=false;
			}
		}
	}

	public void fillholes(byte[] data){
		// here we fill all non-edge connected components
		// first fill the edge connected components with 127
		for(int y=0;y<height;y++){
			if(data[y*width]==(byte)0)
				floodfill(data,(byte)127,0,y);
			if(data[width-1+width*y]==(byte)255)
				floodfill(data,(byte)127,width-1,y);
		}
		for(int x=0;x<width;x++){
			if(data[x]==(byte)0)
				floodfill(data,(byte)127,x,0);
			if(data[x+width*(height-1)]==(byte)0)
				floodfill(data,(byte)127,x,height-1);
		}
		// if things were filled, set them back to zero, if not set them to 255
		for(int i=0;i<width*height;i++){
			if(data[i]==127){
				data[i]=(byte)0;
			}else{
				data[i]=(byte)255;
			}
		}
	}

	final void push(int x,int y){
		stackSize++;
		if(stackSize==maxStackSize){
			int[] newXStack=new int[maxStackSize*2];
			int[] newYStack=new int[maxStackSize*2];
			System.arraycopy(xstack,0,newXStack,0,maxStackSize);
			System.arraycopy(ystack,0,newYStack,0,maxStackSize);
			xstack=newXStack;
			ystack=newYStack;
			maxStackSize*=2;
		}
		xstack[stackSize-1]=x;
		ystack[stackSize-1]=y;
	}

	final int popx(){
		if(stackSize==0)
			return -1;
		else
			return xstack[stackSize-1];
	}

	final int popy(){
		int value=ystack[stackSize-1];
		stackSize--;
		return value;
	}

	final void fillLine(byte[] data,byte color,int x1,int x2,int y){
		if(x1>x2){
			int t=x1;
			x1=x2;
			x2=t;
		}
		for(int x=x1;x<=x2;x++){
			data[x+y*width]=color;
		}
	}

}
