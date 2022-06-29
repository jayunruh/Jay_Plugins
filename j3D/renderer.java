/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import jalgs.jsort;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.image.MemoryImageSource;
import java.awt.image.PixelGrabber;

public class renderer implements Cloneable{
	// this class returns an awt image based on transformations of 3D spots and
	// lines
	public element3D[] elements;
	public int[] order;
	public int width,height,depth,centerx,centery,centerz;
	public double degx,degy,degz,horizon_dist;
	public Color background;

	public renderer(int width1,int height1){
		width=width1;
		centerx=width/2;
		height=height1;
		centery=height/2;
		centerz=height/2;
		degx=0.0;
		degy=0.0;
		degz=0.0;
		background=Color.white;
	}

	public void addelement(element3D element){
		// element.translate(width/2,width/2,height/2);
		if(elements==null){
			elements=new element3D[1];
			elements[0]=element;
		}else{
			element3D[] tempelements=new element3D[elements.length+1];
			System.arraycopy(elements,0,tempelements,0,elements.length);
			tempelements[elements.length]=element;
			elements=tempelements;
		}
	}
	
	public void insertelement(element3D element,int pos){
		//inserts an element after the element at pos
		// element.translate(width/2,width/2,height/2);
		if(pos>=(elements.length-1)) addelement(element);
		element3D[] tempelements=new element3D[elements.length+1];
		System.arraycopy(elements,0,tempelements,0,pos+1);
		tempelements[pos+1]=element;
		System.arraycopy(elements,pos+1,tempelements,pos+2,elements.length-(pos+1));
		elements=tempelements;
	}
	
	public void replaceelement(element3D element,int pos){
		elements[pos]=element;
	}
	
	public void removeelement(int remid){
		// element.translate(width/2,width/2,height/2);
		if(remid>=elements.length) return;
		if(elements!=null){
			element3D[] tempelements=new element3D[elements.length-1];
			for(int i=0;i<remid;i++) tempelements[i]=elements[i];
			for(int i=(remid+1);i<elements.length;i++){
				tempelements[i-1]=elements[i];
			}
			elements=tempelements;
		}
	}

	public void addText3D(String s,int x,int y,int z,Color color){
		addelement(new text3D(s,x,y,z,color));
	}

	public void addLine3D(int x1,int y1,int z1,int x2,int y2,int z2,Color color){
		addelement(new line3D(x1,y1,z1,x2,y2,z2,color));
	}

	public void addPoint3D(int x,int y,int z,int shape,Color color){
		addelement(new spot3D(x,y,z,shape,color));
	}

	public void addCube3D(int x,int y,int z,int size,Color color){
		addelement(new cube3D(x,y,z,size,color));
	}

	public void setelementarray(element3D[] elements1){
		elements=elements1;
	}

	public void addpointarray(int[] x,int[] y,int[] z,int shape,Color color){
		if(elements==null){
			elements=new element3D[x.length];
			for(int i=0;i<x.length;i++){
				elements[i]=new spot3D(x[i],y[i],z[i],shape,color);
			}
		}else{
			element3D[] tempelements=new element3D[elements.length+x.length];
			System.arraycopy(elements,0,tempelements,0,elements.length);
			for(int i=0;i<x.length;i++){
				tempelements[elements.length+i]=new spot3D(x[i],y[i],z[i],shape,color);
			}
			elements=tempelements;
		}
	}

	public void addpolyline(int[] x,int[] y,int[] z,Color color){
		if(elements==null){
			elements=new element3D[x.length-1];
			for(int i=0;i<x.length-1;i++){
				elements[i]=new line3D(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1],color);
			}
		}else{
			element3D[] tempelements=new element3D[elements.length+x.length-1];
			System.arraycopy(elements,0,tempelements,0,elements.length);
			for(int i=0;i<x.length-1;i++){
				tempelements[elements.length+i]=new line3D(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1],color);
			}
			elements=tempelements;
		}
	}

	public void setrotation(int degx1,int degy1,int degz1){
		degx=degx1;
		degy=degy1;
		degz=degz1;
		double cosx=Math.cos(Math.toRadians(degx));
		double cosy=Math.cos(Math.toRadians(degy));
		double cosz=Math.cos(Math.toRadians(-degz));
		double sinx=Math.sin(Math.toRadians(degx));
		double siny=Math.sin(Math.toRadians(degy));
		double sinz=Math.sin(Math.toRadians(-degz));
		for(int i=0;i<elements.length;i++){
			elements[i].rotatecossin(cosx,cosy,cosz,sinx,siny,sinz,centerx,centery,centerz);
		}
	}

	public void setrotation(float degx1,float degy1,float degz1){
		degx=degx1;
		degy=degy1;
		degz=degz1;
		double cosx=Math.cos(Math.toRadians(degx));
		double cosy=Math.cos(Math.toRadians(degy));
		double cosz=Math.cos(Math.toRadians(-degz));
		double sinx=Math.sin(Math.toRadians(degx));
		double siny=Math.sin(Math.toRadians(degy));
		double sinz=Math.sin(Math.toRadians(-degz));
		for(int i=0;i<elements.length;i++){
			elements[i].rotatecossin(cosx,cosy,cosz,sinx,siny,sinz,centerx,centery,centerz);
		}
	}
	
	public void set_obs_rot(double ux,double uy,double uz,double vx,double vy,double vz){
		//this rotates the object from 0,0,0 so that it is observed from vector, u
		//the v vector is orthogonal and denotes the up direction
		//rotate about the cross product between the oberver and the z axis
		if(ux==0.0 && uy==0.0 && uz==1.0){
			setrotation(0.0f,0.0f,0.0f);
			return;
		}
		double[] cross={uy,-ux,0.0};
		double norm=Math.sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
		cross[0]/=norm; cross[1]/=norm;
		double sintheta=norm;
		double costheta=uz;
		for(int i=0;i<elements.length;i++){
			elements[i].rotatecossinabout(costheta,-sintheta,cross[0],cross[1],cross[2],centerx,centery,centerz);
		}
		//rotate the observation vectors as well
		//point3D u=new point3D((int)(10000.0*ux),(int)(10000.0*uy),(int)(10000.0*uz));
		//u.set_rot_about_vector(costheta,-sintheta,cross[0],cross[1],cross[2],0,0,0);
		point3D v=new point3D((int)(10000.0*vx),(int)(10000.0*vy),(int)(10000.0*vz));
		v.set_rot_about_vector(costheta,-sintheta,cross[0],cross[1],cross[2],0,0,0);
		//now rotate about the obervation axis so that the vertical direction is preserved
		double newlength=Math.sqrt(v.rx*v.rx+v.ry*v.ry);
		double costheta2=v.ry/newlength;
		double sintheta2=v.rx/newlength;
		for(int i=0;i<elements.length;i++){
			elements[i].addrotatecossin(0,0,costheta2,0,0,-sintheta2,centerx,centery,centerz);
		}
	}

	public void set_perspective(double horizon_dist1){
		horizon_dist=horizon_dist1;
		for(int i=0;i<elements.length;i++){
			elements[i].transform_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void set_negative_perspective(double horizon_dist1){
		horizon_dist=horizon_dist1;
		for(int i=0;i<elements.length;i++){
			elements[i].transform_negative_perspective(horizon_dist,centerx,centery,centerz);
		}
	}

	public void rotate(int dx,int dy,int dz){
		degx+=dx;
		degy+=dy;
		degz+=dz;
		for(int i=0;i<elements.length;i++){
			elements[i].rotatedeg(degx,degy,degz,centerx,centery,centerz);
		}
	}

	public void addrotation(int dx,int dy,int dz){
		for(int i=0;i<elements.length;i++){
			elements[i].addrotatedeg(dx,dy,dz,centerx,centery,centerz);
		}
	}
	
	public void transform(double[][] transmat){
		for(int i=0;i<elements.length;i++){
			elements[i].transform(transmat,centerx,centery,centerz);
		}
	}
	
	public void translate(int dx,int dy,int dz){
		for(int i=0;i<elements.length;i++){
			elements[i].translate(dx,dy,dz);
		}
	}

	private void setyorder(){
		// here we sort the element list in order of decreasing z value
		// that way closer z values will be drawn last, on top of further z
		// values
		float[] zvals=new float[elements.length];
		for(int i=0;i<elements.length;i++){
			zvals[i]=elements[i].getzpos();
		}
		int[] temporder=jsort.javasort_order(zvals);
		order=new int[elements.length];
		for(int i=0;i<elements.length;i++){
			order[i]=temporder[elements.length-i-1];
		}
	}

	public Image renderimage(){
		Image retimage=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics g=retimage.getGraphics();
		Color tempcolor=g.getColor();
		g.setColor(background);
		g.fillRect(0,0,width,height);
		g.setColor(tempcolor);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g);
		}
		return retimage;
	}
	
	/*public Image[] render3D(){
		//this renders the entire 3D volume as an array of images
		if(depth==0) depth=32;
		Image[] stack=new Image[depth];
		for(int i=0;i<depth;i++) stack[i]=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics[] g=new Graphics[depth];
		for(int i=0;i<depth;i++){
			g[i]=stack[i].getGraphics();
			Color tempcolor=g[i].getColor();
			g[i].setColor(background);
			g[i].fillRect(0,0,width,height);
			g[i].setColor(tempcolor);
			setyorder();
		}
		for(int i=0;i<elements.length;i++){
			elements[order[i]].draw3Delement(g);
		}
		return stack;
	}*/

	public byte[] renderEMF(){
		RenderVector rv=new RenderVector(width,height,background,0);
		setyorder();
		for(int i=0;i<elements.length;i++){
			rv.drawElement(elements[order[i]]);
		}
		rv.endExport();
		return rv.getByteArray();
	}
	
	public byte[] renderPS(){
		RenderVector rv=new RenderVector(width,height,background,1);
		setyorder();
		for(int i=0;i<elements.length;i++){
			rv.drawElement(elements[order[i]]);
		}
		rv.endExport();
		return rv.getByteArray();
	}
	
	public byte[] renderPDF(){
		RenderVector rv=new RenderVector(width,height,background,2);
		setyorder();
		for(int i=0;i<elements.length;i++){
			rv.drawElement(elements[order[i]]);
		}
		rv.endExport();
		return rv.getByteArray();
	}
	
	public int[] renderrgb(){
		Image retimage=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics g=retimage.getGraphics();
		Color tempcolor=g.getColor();
		g.setColor(background);
		g.fillRect(0,0,width,height);
		g.setColor(tempcolor);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g);
		}
		int[] pixels=((BufferedImage)retimage).getRGB(0,0,width,height,null,0,width);
		return pixels;
	}

	public float[] renderfloat(float value){
		Image retimage=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics g=retimage.getGraphics();
		Color tempcolor=g.getColor();
		g.setColor(background);
		g.fillRect(0,0,width,height);
		g.setColor(tempcolor);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g);
		}
		int[] pixels=((BufferedImage)retimage).getRGB(0,0,width,height,null,0,width);
		float[] output=new float[width*height];
		for(int i=0;i<width*height;i++){
			if((pixels[i]&0xffffff)>0)
				output[i]=value;
		}
		return output;
	}

	public Image renderanalglyph(Color leftcolor,Color rightcolor,double angle){
		Color[] colors=new Color[elements.length];
		for(int i=0;i<elements.length;i++){
			colors[i]=copycolor(elements[i].color);
		}
		Color[] leftcolors=new Color[elements.length];
		for(int i=0;i<elements.length;i++){
			leftcolors[i]=scalecolor(leftcolor,colormagnitude(colors[i]));
		}
		Color[] rightcolors=new Color[elements.length];
		for(int i=0;i<elements.length;i++){
			rightcolors[i]=scalecolor(rightcolor,colormagnitude(colors[i]));
		}
		// set up the drawing
		Image retimage1=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Image retimage2=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics g1=retimage1.getGraphics();
		Color tempcolor=g1.getColor();
		g1.setColor(background);
		g1.fillRect(0,0,width,height);
		g1.setColor(tempcolor);
		Graphics g2=retimage2.getGraphics();
		tempcolor=g2.getColor();
		g2.setColor(background);
		g2.fillRect(0,0,width,height);
		g2.setColor(tempcolor);

		// start by drawing the left color
		// for(int
		// i=0;i<elements.length;i++){elements[i].color=scalecolor(leftcolor,colormagnitude(colors[i]));}
		for(int i=0;i<elements.length;i++){
			elements[i].color=leftcolor;
		}
		addrotation(0,-(int)(angle/2.0),0);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g1);
		}

		// now draw the right color
		// for(int
		// i=0;i<elements.length;i++){elements[i].color=scalecolor(rightcolor,colormagnitude(colors[i]));}
		for(int i=0;i<elements.length;i++){
			elements[i].color=rightcolor;
		}
		setrotation((int)degx,(int)degy,(int)degz);
		set_perspective(horizon_dist);
		addrotation(0,(int)(angle/2.0),0);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g2);
		}
		// now reset everything
		for(int i=0;i<elements.length;i++){
			elements[i].color=colors[i];
		}
		setrotation((int)degx,(int)degy,(int)degz);
		set_perspective(horizon_dist);
		return combinecoloranalglyphs(retimage1,retimage2);
	}

	public Image renderinterlacedanalglyph(double angle){
		return renderinterlacedanalglyph(angle,0);
	}

	public Image renderinterlacedanalglyph(double angle,int interlacetype){
		// note that interlace types are 0:diagonal,1:vertical,and 2:horizontal
		// set up the drawing
		Image retimage1=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Image retimage2=new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
		Graphics g1=retimage1.getGraphics();
		Color tempcolor=g1.getColor();
		g1.setColor(background);
		g1.fillRect(0,0,width,height);
		g1.setColor(tempcolor);
		Graphics g2=retimage2.getGraphics();
		tempcolor=g2.getColor();
		g2.setColor(background);
		g2.fillRect(0,0,width,height);
		g2.setColor(tempcolor);

		// start by drawing the left image
		double tempdegy=degy;
		addrotation(0,-(int)(angle/2.0),0);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g1);
		}

		// now draw the right image
		setrotation((int)degx,(int)degy,(int)degz);
		set_perspective(horizon_dist);
		addrotation(0,(int)(angle/2.0),0);
		setyorder();
		for(int i=0;i<elements.length;i++){
			elements[order[i]].drawelement(g2);
		}
		// now reset everything
		setrotation((int)degx,(int)tempdegy,(int)degz);
		set_perspective(horizon_dist);
		return interlaceanalglyphs(retimage1,retimage2,interlacetype);
	}

	public Image combinecoloranalglyphs(Image leftimage,Image rightimage){
		// here we add two images typically orthogonal color spaces should be
		// used
		int width=leftimage.getWidth(null);
		int height=leftimage.getHeight(null);
		int[] pix1=getimagepixels(leftimage);
		int[] pix2=getimagepixels(rightimage);
		int[] retpix=new int[width*height];
		for(int i=0;i<width;i++){
			// run the left image first
			for(int j=0;j<height;j++){
				retpix[i+j*width]=maxcolor(pix1[i+j*width],pix2[i+j*width]);
			}
		}
		return getimagefrompixels(retpix,width,height);
	}

	public Image interlaceanalglyphs(Image leftimage,Image rightimage,int interlace){
		int width=leftimage.getWidth(null);
		int height=leftimage.getHeight(null);
		int[] pix1=getimagepixels(leftimage);
		int[] pix2=getimagepixels(rightimage);
		int[] retpix=new int[width*height];
		// eliminate every other line and every other row
		if(interlace<2){
			boolean even=true;
			boolean vinterlace=(interlace==1);
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j+=2){
					int pixval1=0;
					int pixval2=0;
					if(vinterlace||even){
						pixval1=maxcolor(pix1[j+i*width],pix1[j+1+i*width]);
						pixval2=maxcolor(pix2[j+i*width],pix2[j+1+i*width]);
					}else{
						pixval1=maxcolor(pix2[j+i*width],pix2[j+1+i*width]);
						pixval2=maxcolor(pix1[j+i*width],pix1[j+1+i*width]);
					}
					retpix[j+i*width]=pixval1;
					retpix[j+1+i*width]=pixval2;
				}
				even=!even;
			}
		}else{
			for(int i=0;i<height;i+=2){
				for(int j=0;j<width;j++){
					int pixval1=0;
					int pixval2=0;
					pixval1=maxcolor(pix1[j+i*width],pix1[j+(i+1)*width]);
					pixval2=maxcolor(pix2[j+i*width],pix2[j+(i+1)*width]);
					retpix[j+i*width]=pixval1;
					retpix[j+1+i*width]=pixval2;
				}
			}
		}
		return getimagefrompixels(retpix,width,height);
	}

	public int maxcolor(int pix1,int pix2){
		boolean white1=(pix1==0xffffffff);
		boolean white2=(pix2==0xffffffff);
		int r1=(pix1&0xff0000)>>16;
		int g1=(pix1&0xff00)>>8;
		int b1=pix1&0xff;
		if(!white2){
			int r2=(pix2&0xff0000)>>16;
			int g2=(pix2&0xff00)>>8;
			int b2=(pix2&0xff);
			if(!white1){
				r1+=r2;
				g1+=g2;
				b1+=b2;
			}else{
				r1=r2;
				g1=g2;
				b1=b2;
			}
		}
		return 0xff000000+(r1<<16)+(g1<<8)+b1;
	}

	public int colormagnitude(Color colorval){
		int r1=colorval.getRed();
		int max=r1;
		int g1=colorval.getGreen();
		if(g1>max){
			max=g1;
		}
		int b1=colorval.getBlue();
		if(b1>max){
			max=b1;
		}
		// int mag=(int)Math.sqrt(r1*r1+g1*g1+b1*b1);
		// if(mag>255){mag=255;}
		// return (int)(0.33f*(r1+g1+b1));
		return max;
	}

	public Color scalecolor(Color colorval,int scale){
		float factor=scale/255.0f;
		int r1=(int)(factor*colorval.getRed());
		int g1=(int)(factor*colorval.getGreen());
		int b1=(int)(factor*colorval.getBlue());
		return new Color(r1,g1,b1);
	}

	public Color copycolor(Color colorval){
		return new Color(colorval.getRed(),colorval.getGreen(),colorval.getBlue());
	}

	public Image getimagefrompixels(int[] pixels,int width,int height){
		MemoryImageSource source=new MemoryImageSource(width,height,pixels,0,width);
		return Toolkit.getDefaultToolkit().createImage(source);
	}

	public int[] getimagepixels(Image input){
		int width=input.getWidth(null);
		int height=input.getHeight(null);
		int[] pixels=new int[width*height];
		PixelGrabber pg=new PixelGrabber(input,0,0,width,height,pixels,0,width);
		try{
			pg.grabPixels();
		}catch(InterruptedException e){}
		return pixels;
	}
	
	public renderer clone(){
		element3D[] temp=new element3D[elements.length];
		for(int i=0;i<temp.length;i++){
			temp[i]=elements[i].clone();
		}
		renderer rtemp=new renderer(width,height);
		rtemp.elements=temp;
		rtemp.centerx=centerx;
		rtemp.centery=centery;
		rtemp.centerz=centerz;
		rtemp.degx=degx;
		rtemp.degy=degy;
		rtemp.degz=degz;
		rtemp.order=order;
		return rtemp;
	}

}
