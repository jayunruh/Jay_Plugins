package jguis;

import org.freehep.graphics2d.*;
import org.freehep.graphicsio.emf.*;
import ij.gui.*;

import java.awt.*;
import java.awt.geom.*;
import java.io.*;
import java.text.DecimalFormat;

public class emf_converter {
	
	private static float magnification,magratio;
	public static int currx,curry;
	public static boolean logx,logy;
	public static float xMin,xMax,yMin,yMax,logxmin,logxmax,logymin,logymax,logxscale,logyscale;
	
	public static void export_plot(ImageWindow iw,String path){
		try{
			File file=new File(path);
			Plot4 p4=(Plot4)jutils.getPW4PlotCopy(iw);
			VectorGraphics graphics=new EMFGraphics2D(file,getsize(p4));
			graphics.startExport();
			//Graphics g=new Graphics();
			//we need to figure out how to create panel independent graphics here
			VectorGraphics g2=(VectorGraphics)graphics.create();
			drawPlot(g2,p4);
			g2.dispose();
			graphics.endExport();
		} catch(IOException e){
			
		}
	}
	
	private static Dimension getsize(Plot4 p4){
		magnification=p4.getmagnification();
		magratio=p4.getmagratio();
		float ymag=magnification*magratio/((float)Plot4.HEIGHT/(float)Plot4.WIDTH);
		int newleftmargin=(int)(magnification*Plot4.LEFT_MARGIN);
		int newtopmargin=(int)(ymag*Plot4.TOP_MARGIN);
		int newwidth=(int)(magnification*Plot4.WIDTH);
		int newheight=(int)(ymag*Plot4.HEIGHT);
		int newrightmargin=(int)(magnification*Plot4.RIGHT_MARGIN);
		int newbottommargin=(int)(ymag*Plot4.BOTTOM_MARGIN);
		int width=newwidth+newleftmargin+newrightmargin;
		int height=newheight+newtopmargin+newbottommargin;
		return new Dimension(width,height);
	}
	
	private static void drawPlot(VectorGraphics graphics,Plot4 p4){
		float magnification=p4.getmagnification();
		float magratio=p4.getmagratio();
		float ymag=magnification*magratio/((float)Plot4.HEIGHT/(float)Plot4.WIDTH);
		int newleftmargin=(int)(magnification*Plot4.LEFT_MARGIN);
		int newtopmargin=(int)(ymag*Plot4.TOP_MARGIN);
		int newwidth=(int)(magnification*Plot4.WIDTH);
		int newheight=(int)(ymag*Plot4.HEIGHT);
		int newrightmargin=(int)(magnification*Plot4.RIGHT_MARGIN);
		int newbottommargin=(int)(ymag*Plot4.BOTTOM_MARGIN);

		Rectangle frame = new Rectangle(newleftmargin,newtopmargin,newwidth,newheight);
		int width=newwidth+newleftmargin+newrightmargin;
		int height=newheight+newtopmargin+newbottommargin;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){temp[i]=0xffffffff;}

    	logxmin=0; logymin=0; logxscale=0; logyscale=0; logxmax=0; logymax=0;
    	boolean[] logs=p4.getLogAxes(); logx=logs[0]; logy=logs[1];
    	float[] limits=p4.getLimits();
    	xMin=limits[0]; xMax=limits[1]; yMin=limits[2]; yMax=limits[1];
    	float[][] xValues=p4.getXValues(); float[][] yValues=p4.getYValues();
    	int[] npts=p4.getNpts();
    	int nseries=p4.getNSeries();
    	int[] colors=p4.getColors();
    	int[] shapes=p4.getShapes();
    	int selected=p4.getSelected();
    	boolean showerrors=p4.getShowErrors();
    	float[][][] errors=p4.getErrors();
    	
    	if(logx){
    		if(xMin<=0.0f){
    			logxmin=(float)Math.log((double)findmingt0(xValues,npts,xMax));
    		} else {
    			logxmin=(float)Math.log((double)xMin);
    		}
    		logxmax=(float)Math.log((double)xMax);
    		logxscale=(float)newwidth/(logxmax-logxmin);
    	}
    	if(logy){
    		if(yMin<=0.0f){
    			logymin=(float)Math.log((double)findmingt0(yValues,npts,yMax));
    		} else {
    			logymin=(float)Math.log((double)yMin);
    		}
    		logymax=(float)Math.log((double)yMax);
    		logyscale=(float)newheight/(logymax-logymin);
    	}
       // IJ.showMessage("testdraw1");
    	
		drawAxisLabels(graphics,p4);
		graphics.setClip(frame);
        float xScale=(float)newwidth/(xMax-xMin);
        float yScale=(float)newheight/(yMax-yMin);
    	
        //IJ.showMessage("testdraw2");
        if(!logx && !logy){
	        for(int j=0;j<nseries;j++){
	        	graphics.setColor(getColor(colors[j]));
	        	int xpoints[] = new int[npts[j]];
	        	int ypoints[] = new int[npts[j]];
		        for (int i=0; i<npts[j]; i++) {
		            xpoints[i] = newleftmargin+ (int)((xValues[j][i]-xMin)*xScale);
		            ypoints[i] = newtopmargin + frame.height - (int)((yValues[j][i]-yMin)*yScale);
		        }
		        if(j!=selected){
			        if(shapes[j]==0){
			        	drawPolyline(graphics, xpoints, ypoints, npts[j]);
			        } else {
			        	drawPolyshape(graphics,xpoints,ypoints,npts[j],shapes[j]);
			        }
		        } else {
			        if(shapes[j]==0){
			        	drawPolyshape(graphics, xpoints, ypoints, npts[j],1);
			        } else {
			        	drawPolyline(graphics,xpoints,ypoints,npts[j]);
			        }
		        }
		        if(showerrors && errors!=null){
		        	int[] yerrptsu=new int[npts[j]];
		        	int[] yerrptsl=new int[npts[j]];
			        for (int i=0; i<npts[j]; i++) {
			            yerrptsu[i] = newtopmargin + frame.height - (int)((yValues[j][i]+errors[1][j][i]-yMin)*yScale);
			            yerrptsl[i] = newtopmargin + frame.height - (int)((yValues[j][i]-errors[0][j][i]-yMin)*yScale);
			        }
		        	drawPolyerrors(graphics,xpoints,yerrptsu,yerrptsl,npts[j]);
		        }
	        }
        } else {
	        for(int j=0;j<nseries;j++){
	        	graphics.setColor(getColor(colors[j]));
	        	int xpoints[] = new int[npts[j]];
	        	int ypoints[] = new int[npts[j]];
		        for (int i=0; i<npts[j]; i++) {
		        	if(logx){
		        		float xtemp;
		        		if(xValues[j][i]>0.0f){
		        			xtemp=(float)Math.log((double)xValues[j][i]);
		        		} else {
		        			xtemp=logxmin;
		        		}
		        		xpoints[i] = newleftmargin + (int)((xtemp-logxmin)*logxscale);
		        	} else {
		        		xpoints[i] = newleftmargin + (int)((xValues[j][i]-xMin)*xScale);
		        	}
		        }
		        for (int i=0; i<npts[j]; i++) {
		        	if(logy){
		        		float ytemp;
		        		if(yValues[j][i]>0.0f){
		        			ytemp=(float)Math.log((double)yValues[j][i]);
		        		} else {
		        			ytemp=logymin;
		        		}
			            ypoints[i] = newtopmargin + frame.height - (int)((ytemp-logymin)*logyscale);
		        	} else {
			            ypoints[i] = newtopmargin + frame.height - (int)((yValues[j][i]-yMin)*yScale);
		        	}
		        }
		        if(j!=selected){
			        if(shapes[j]==0){
			        	drawPolyline(graphics, xpoints, ypoints, npts[j]);
			        } else {
			        	drawPolyshape(graphics,xpoints,ypoints,npts[j],shapes[j]);
			        }
		        }else{
			        if(shapes[j]==0){
			        	drawPolyshape(graphics, xpoints, ypoints, npts[j],1);
			        } else {
			        	drawPolyline(graphics, xpoints, ypoints, npts[j]);
			        }
		        }
		        if(showerrors && errors!=null){
		        	int[] yerrptsu=new int[npts[j]];
		        	int[] yerrptsl=new int[npts[j]];
		        	if(logy){
		        		for (int i=0; i<npts[j]; i++) {
		        			float ytemp=yValues[j][i]+errors[1][j][i];
		        			if(ytemp>0.0f){
		        				ytemp=(float)Math.log((double)ytemp);
		        			} else {
		        				ytemp=logymin;
		        			}
		        			yerrptsu[i] = newtopmargin + frame.height - (int)((ytemp-logymin)*logyscale);
		        			ytemp=yValues[j][i]-errors[1][j][i];
		        			if(ytemp>0.0f){
		        				ytemp=(float)Math.log((double)ytemp);
		        			} else {
		        				ytemp=logymin;
		        			}
		        			yerrptsl[i] = newtopmargin + frame.height - (int)((ytemp-logymin)*yScale);
		        		}
		        	} else {
		        		for (int i=0; i<npts[j]; i++) {
		        			yerrptsu[i] = newtopmargin + frame.height - (int)((yValues[j][i]+errors[1][j][i]-yMin)*yScale);
		        			yerrptsl[i] = newtopmargin + frame.height - (int)((yValues[j][i]-errors[0][j][i]-yMin)*yScale);
		        		}
		        	}
		        	drawPolyerrors(graphics,xpoints,yerrptsu,yerrptsl,npts[j]);
		        }
	        }
        }
        graphics.setColor(Color.black);
        graphics.setClip(null);
        //return cp;
        //IJ.showMessage("testdraw3");
	}
	
    private static void drawAxisLabels(VectorGraphics g,Plot4 p4) {
    	//set up formating
    	DecimalFormat expformat=new DecimalFormat("0.00E0");
		float ymag=magnification*magratio/((float)Plot4.HEIGHT/(float)Plot4.WIDTH);
		int newleftmargin=(int)(magnification*Plot4.LEFT_MARGIN);
		int newtopmargin=(int)(ymag*Plot4.TOP_MARGIN);
		int newwidth=(int)(magnification*Plot4.WIDTH);
		int newheight=(int)(ymag*Plot4.HEIGHT);
		int newrightmargin=(int)(magnification*Plot4.RIGHT_MARGIN);
		int newbottommargin=(int)(ymag*Plot4.BOTTOM_MARGIN);
		int newticklength=(int)(magnification*Plot4.TICK_LENGTH);
		int newfontsize=(int)(magnification*Plot4.fontsize);
		
		//calculate the appropriate label numbers
		float[] xticklabels=new float[4];
		float[] yticklabels=new float[4];
		for(int i=0;i<4;i++){
			if(logx){
				float tempx=logxmin+((float)i/3.0f)*(logxmax-logxmin);
				xticklabels[i]=(float)Math.exp((double)tempx);
			} else {
				xticklabels[i]=xMin+((float)i/3.0f)*(xMax-xMin);
			}
			if(logy){
				float tempy=logymin+((float)i/3.0f)*(logymax-logymin);
				yticklabels[i]=(float)Math.exp((double)tempy);
			} else {
				yticklabels[i]=yMin+((float)i/3.0f)*(yMax-yMin);
			}
		}
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		Font tempfont=new Font("SansSerif", Font.PLAIN, newfontsize);
		g.setFont(tempfont);
		FontMetrics fm = g.getFontMetrics();
		int fontheight=fm.getAscent()+fm.getDescent();
		
    	//draw the y axis labels

		String s = jutils.formatted_string((double)yticklabels[0]);
		g.drawString(s, newleftmargin-getStringWidth(g,s)-newticklength-2, newtopmargin+newheight+fm.getDescent());
		s = jutils.formatted_string((double)yticklabels[1]);
		g.drawString(s,  newleftmargin-getStringWidth(g,s)-newticklength-2, newtopmargin+(int)((2*newheight)/3)+fontheight/2);
		s = jutils.formatted_string((double)yticklabels[2]);
		g.drawString(s,  newleftmargin-getStringWidth(g,s)-newticklength-2, newtopmargin+(int)(newheight/3)+fontheight/2);
    	s=jutils.formatted_string((double)yticklabels[3]);
		g.drawString(s,  newleftmargin-getStringWidth(g,s)-newticklength-2, newtopmargin+fontheight/2);
		
		//draw the x axis labels
		int y = newtopmargin + newheight + fontheight + newticklength+4;
		s = jutils.formatted_string((double)xticklabels[0]);
		g.drawString(s, newleftmargin-8, y);
		s = jutils.formatted_string((double)xticklabels[1]);
		g.drawString(s, newleftmargin+ (int)(newwidth/3-getStringWidth(g,s)/2), y);
		s = jutils.formatted_string((double)xticklabels[2]);
		g.drawString(s, newleftmargin + (int)((2*newwidth)/3-getStringWidth(g,s)/2), y);
		s = jutils.formatted_string((double)xticklabels[3]);
		g.drawString(s, newleftmargin + newwidth-getStringWidth(g,s)+8, y);
		
		//now draw the axis labels
		g.drawString(p4.getxLabel(), newleftmargin+(newwidth-getStringWidth(g,p4.getxLabel()))/2, newtopmargin+newheight+fm.getAscent()+newbottommargin/2);
		drawVerticalString(p4.getyLabel(),g,newleftmargin-(4*newleftmargin/5),newtopmargin+(int)(newheight/2));
		//now draw the tick marks and grid lines
		g.drawRect(newleftmargin, newtopmargin, newwidth+1, newheight+1);
		g.setColor(p4.gridColor);
		g.drawLine(newleftmargin+(int)(newwidth/3), newtopmargin, newleftmargin+(int)(newwidth/3), newtopmargin+newheight);
		g.drawLine(newleftmargin+(int)((2*newwidth)/3), newtopmargin, newleftmargin+(int)((2*newwidth)/3), newtopmargin+newheight);
		g.drawLine(newleftmargin, newtopmargin+(int)(newheight/3), newleftmargin+newwidth, newtopmargin+(int)(newheight/3));
		g.drawLine(newleftmargin, newtopmargin+(int)((2*newheight)/3), newleftmargin+newwidth, newtopmargin+(int)((2*newheight)/3));
		g.setColor(Color.black);
		g.drawLine(newleftmargin, newtopmargin, newleftmargin-newticklength, newtopmargin);
		g.drawLine(newleftmargin, newtopmargin+(int)(newheight/3), newleftmargin-newticklength, newtopmargin+(int)(newheight/3));
		g.drawLine(newleftmargin, newtopmargin+(int)((2*newheight)/3), newleftmargin-newticklength, newtopmargin+(int)((2*newheight)/3));
		g.drawLine(newleftmargin, newtopmargin+newheight, newleftmargin-newticklength, newtopmargin+newheight);
		
		g.drawLine(newleftmargin, newtopmargin+newheight, newleftmargin, newtopmargin+newheight+newticklength);
		g.drawLine(newleftmargin+(int)(newwidth/3), newtopmargin+newheight, newleftmargin+(int)(newwidth/3), newtopmargin+newheight+newticklength);
		g.drawLine(newleftmargin+(int)((2*newwidth)/3), newtopmargin+newheight, newleftmargin+(int)((2*newwidth)/3), newtopmargin+newheight+newticklength);
		g.drawLine(newleftmargin+newwidth, newtopmargin+newheight, newleftmargin+newwidth, newtopmargin+newheight+newticklength);
    }
    
    private static int getStringWidth(VectorGraphics g,String s){
    	return g.getFontMetrics().stringWidth(s);
    }
    
    private static void drawVerticalString(String text,VectorGraphics g,int x,int y){
    	Font currfont=g.getFont();
    	FontMetrics fm=g.getFontMetrics();
    	int ascent=fm.getAscent();
    	int descent=fm.getDescent();
    	int height=ascent+descent;
    	int width=getStringWidth(g,text);
    	AffineTransform at=new AffineTransform();
    	at.rotate(-Math.PI/2.0);
    	Font dfont=currfont.deriveFont(at);
    	g.setFont(dfont);
    	g.drawString(text,x-height/2,y-width/2);
    	g.setFont(currfont);
    }
	
    private static Color getColor(int index){
		int temp=index;
		if(temp>=8){temp=index%8;}
		Color[] temp2={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};
		return temp2[temp];
	}
	
	private static void drawPolyline(VectorGraphics g, int[] x, int[] y, int n) {
		g.drawPolyline(x,y,n);
    }
    
    private static void drawPolyshape(VectorGraphics g, int[] x, int[] y, int n, int shape1) {
    	int newshapesize=(int)(magnification*Plot4.shapesize);
        for (int i=0; i<n; i++){
            if(shape1==1){drawSquare(g,x[i],y[i],newshapesize);}
        	if(shape1==2){drawPlus(g,x[i],y[i],newshapesize);}
        	if(shape1==3){drawX(g,x[i],y[i],newshapesize);}
        	if(shape1==4){drawTriangle(g,x[i],y[i],newshapesize);}
        }
    }
    
    private static void drawPolyerrors(VectorGraphics g,int[] x,int[] yu,int[] yl,int n){
    	int newshapesize=(int)(magnification*Plot4.shapesize);
    	for(int i=0;i<n;i++){
    		if(yu[i]!=yl[i]){
    			g.drawLine(x[i]-newshapesize/2, yl[i], x[i]+newshapesize/2,yl[i]);
    			g.drawLine(x[i], yl[i],x[i], yu[i]);
    			g.drawLine(x[i]-newshapesize/2, yu[i],x[i]+newshapesize/2, yu[i]);
    		}
    	}
    }
    
    private static void drawSquare(VectorGraphics g,int x,int y,int size){
    	moveTo(x-size/2, y-size/2);
    	lineTo(g,x+size/2, y-size/2);
    	lineTo(g,x+size/2, y+size/2);
    	lineTo(g,x-size/2, y+size/2);
    	lineTo(g,x-size/2, y-size/2);
    }
    
    private static void drawPlus(VectorGraphics g,int x,int y,int size){
    	moveTo(x-size/2, y);
    	lineTo(g,x+size/2, y);
    	moveTo(x, y-size/2);
    	lineTo(g,x, y+size/2);
    }
    
    private static void drawX(VectorGraphics g,int x,int y,int size){
    	moveTo(x-size/2, y-size/2);
    	lineTo(g,x+size/2, y+size/2);
    	moveTo(x-size/2, y+size/2);
    	lineTo(g,x+size/2, y-size/2);
    }
    
    private static void drawTriangle(VectorGraphics g,int x,int y,int size){
    	moveTo(x, y-size/2);
    	lineTo(g,x-size/2, y+size/2);
    	lineTo(g,x+size/2, y+size/2);
    	lineTo(g,x, y-size/2);
    }
    
    private static void moveTo(int x,int y){
    	currx=x;
    	curry=y;
    }
    
    private static void lineTo(VectorGraphics g,int x,int y){
    	g.drawLine(currx,curry,x,y);
    	currx=x;
    	curry=y;
    }
    
	private static float[] findminmax(float[][] arr,int[] npts1){
		float[] temp=new float[2];
		temp[0]=arr[0][0];
		temp[1]=arr[0][0];
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[i];j++){
				if(arr[i][j]<temp[0]){temp[0]=arr[i][j];}
				if(arr[i][j]>temp[1]){temp[1]=arr[i][j];}
			}
		}
		return temp;
	}
	
	private static float findmingt0(float[][] arr,int[] npts1,float max){
		float temp=max;
		if(max<=0.0f){return 0.0f;}
		for(int i=0;i<arr.length;i++){
			for(int j=0;j<npts1[i];j++){
				if(arr[i][j]<temp && arr[i][j]>0.0f){temp=arr[i][j];}
			}
		}
		return temp;
	}
	
	

}
