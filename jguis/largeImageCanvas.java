package jguis;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

public class largeImageCanvas extends ImageCanvas{
	//this class attempts to draw a downsampled version of images to prevent display lock-up for large images
	//call the ImageWindow(imp,canvas) with this class to utilize
	private boolean painted;
	public int maxpixels;
	public boolean dynamicscale;

	public largeImageCanvas(ImagePlus imp){
		super(imp);
		IJ.register(this.getClass());
		maxpixels=512;
		dynamicscale=false;
		// TODO Auto-generated constructor stub
	}
	
    public void paint(Graphics g) {
        if (IJ.debugMode) IJ.log("ImageCanvas.paint: "+imp);
        painted = true;
        //IJ.log("setting painted");
        //jutils.setReflectionField(this,"painted",new Boolean(true));
        Roi roi = imp.getRoi(); 
        //eliminate double buffered painting--shouldn't need it
        /*if (roi!=null) {
            // Use double buffering to avoid flickering of ROIs and to work around
            // a Linux problem with large images not showing at low magnification.
            if (roi!=null) {
            	jutils.runReflectionMethod(roi,"updatePaste",null);
                //roi.updatePaste();
            }
        }*/
        try {
            if (imageUpdated) {
                imageUpdated = false;
                imp.updateImage();
            }
            //setInterpolation(g, Prefs.interpolateScaledImages);
            //not sure if we need to interpolate
            //jutils.runReflectionMethod(this,"setInterpolation",new Object[] {g,new Boolean(Prefs.interpolateScaledImages)});
            Image img = imp.getImage();
            if (img!=null) {
                //g.drawImage(img, 0, 0, (int)(srcRect.width*magnification+0.5), (int)(srcRect.height*magnification+0.5),srcRect.x, srcRect.y, srcRect.x+srcRect.width, srcRect.y+srcRect.height, null);
            	//params are img,dx1,dy1,dx2,dy2,sx1,sy1,sx2,sy2
                //here's where we need to change things
                //need to resample the img object at maxpixels pixels (along the large dimension)
                //then calculate the new srcRect information
            	//boolean dynamicscale=false;
            	int owidth=imp.getWidth();
            	int oheight=imp.getHeight();
            	if(dynamicscale) {
            		owidth=srcRect.width;
            		oheight=srcRect.height;
            	}
            	int newwidth=maxpixels; int newheight=maxpixels;
            	float stride=1.0f;
            	if(owidth>oheight) {
            		stride=(float)owidth/(float)newwidth;
            		newheight=(int)((float)oheight/stride);

            	} else {
            		stride=(float)oheight/(float)newheight;
            		newwidth=(int)((float)owidth/stride);
            	}
            	Rectangle newrect=new Rectangle();
        		newrect.x=(int)((float)srcRect.x/stride);
        		newrect.y=(int)((float)srcRect.y/stride);
        		newrect.width=(int)((float)srcRect.width/stride);
        		newrect.height=(int)((float)srcRect.height/stride);
        		Image newimg=null;
                if(!dynamicscale) newimg=resampleImage(img,imp.getWidth(),imp.getHeight(),newwidth,newheight,stride);
                else newimg=resampleImage(img,imp.getWidth(),imp.getHeight(),newrect.width,newrect.height,stride,srcRect.x,srcRect.y);
                if(!dynamicscale) g.drawImage(newimg,0,0,(int)(srcRect.width*magnification+0.5),(int)(srcRect.height*magnification+0.5),newrect.x,newrect.y,newrect.x+newrect.width,newrect.y+newrect.height,null);
                else g.drawImage(newimg,0,0,(int)(srcRect.width*magnification+0.5),(int)(srcRect.height*magnification+0.5),0,0,newrect.width,newrect.height,null);
            }
            Overlay overlay=super.getOverlay();
            Overlay showalloverlay=super.getShowAllList();
            /*
            //for now skip drawing the overlay--can add that later
            if (overlay!=null)
                drawOverlay(overlay, g);
            if (showAllOverlay!=null)
                drawOverlay(showAllOverlay, g);
            */
            if (roi!=null) {
            	//drawRoi(roi, g);
            	roi.draw(g);
            	//jutils.runReflectionMethod(this,"drawRoi",new Object[] {roi,g});
            }
            if (srcRect.width<imageWidth || srcRect.height<imageHeight) {
                drawZoomIndicator(g);
            	//jutils.runReflectionMethod(this,"drawZoomIndicator",new Object[] {g});
            }
            //if (IJ.debugMode) showFrameRate(g);
        }
        catch(OutOfMemoryError e) {IJ.outOfMemory("Paint");}
        setPaintPending(false);
    }
    
    /*private void drawOverlay(Overlay overlay, Graphics g) {
    	Overlay showAllOverlay=super.getShowAllList();
        if (imp!=null && imp.getHideOverlay() && overlay!=showAllOverlay)
            return;
        boolean flattening = imp!=null && ImagePlus.flattenTitle.equals(imp.getTitle());
        if (imp!=null && showAllOverlay!=null && overlay!=showAllOverlay)
            overlay.drawLabels(false);
        Color labelColor = overlay.getLabelColor();
        if (labelColor==null) labelColor = Color.white;
        initGraphics(overlay, g, labelColor, Roi.getColor());
        int n = overlay.size();
        if (IJ.debugMode) IJ.log("drawOverlay: "+n);
        int currentImage = imp!=null?imp.getCurrentSlice():-1;
        int stackSize = imp.getStackSize();
        if (stackSize==1)
            currentImage = -1;
        int channel=0, slice=0, frame=0;
        boolean hyperstack = imp.isHyperStack();
        if (hyperstack) {
            channel = imp.getChannel();
            slice = imp.getSlice();
            frame = imp.getFrame();
        }
        drawNames = overlay.getDrawNames() && overlay.getDrawLabels();
        boolean drawLabels = drawNames || overlay.getDrawLabels();
        if (drawLabels)
            labelRects = new Rectangle[n];
        else
            labelRects = null;
        font = overlay.getLabelFont();
        if (overlay.scalableLabels() && font!=null) {
            double mag = getMagnification();
            if (mag!=1.0)
                font = font.deriveFont((float)(font.getSize()*mag));
        }
        Roi activeRoi = imp.getRoi();
        boolean roiManagerShowAllMode = overlay==showAllOverlay && !Prefs.showAllSliceOnly;
        for (int i=0; i<n; i++) {
            if (overlay==null) break;
            Roi roi = overlay.get(i);
            if (roi==null) break;
            if (hyperstack) {
                int c = roi.getCPosition();
                int z = roi.getZPosition();
                int t = roi.getTPosition();
                int position = roi.getPosition();
                //IJ.log(c+" "+z+" "+t+"  "+position+" "+roiManagerShowAllMode);
                if (position>0) {
                    if (z==0 && imp.getNSlices()>1)
                        z = position;
                    else if (t==0)
                        t = position;
                }
                if (((c==0||c==channel) && (z==0||z==slice) && (t==0||t==frame)) || roiManagerShowAllMode)
                    drawRoi(g, roi, drawLabels?i+LIST_OFFSET:-1);
            } else {
                int position =  stackSize>1?roi.getPosition():0;
                if (position==0 && stackSize>1)
                    position = getSliceNumber(roi.getName());
                if (position>0 && imp.getCompositeMode()==IJ.COMPOSITE)
                    position = 0;
                //IJ.log(position+"  "+currentImage+" "+roiManagerShowAllMode);
                if (position==0 || position==currentImage || roiManagerShowAllMode)
                    drawRoi(g, roi, drawLabels?i+LIST_OFFSET:-1);
            }
        }
        ((Graphics2D)g).setStroke(Roi.onePixelWide);
        drawNames = false;
        font = null;
    }
    
    private void initGraphics(Overlay overlay, Graphics g, Color textColor, Color defaultColor) {
        if (smallFont==null) {
            smallFont = new Font("SansSerif", Font.PLAIN, 9);
            largeFont = ImageJ.SansSerif12;
        }
        if (textColor!=null) {
            labelColor = textColor;
            if (overlay!=null && overlay.getDrawBackgrounds())
                bgColor = new Color(255-labelColor.getRed(), 255-labelColor.getGreen(), 255-labelColor.getBlue());
            else
                bgColor = null;
        } else {
            int red = defaultColor.getRed();
            int green = defaultColor.getGreen();
            int blue = defaultColor.getBlue();
            if ((red+green+blue)/3<128)
                labelColor = Color.white;
            else
                labelColor = Color.black;
            bgColor = defaultColor;
        }
        this.defaultColor = defaultColor;
        g.setColor(defaultColor);
    }*/
    
    /*private void drawRoi(Roi roi, Graphics g) {
    	//getting the currentRoi
    	Roi currentRoi=(Roi)jutils.getReflectionField(this,"currentRoi");
        if (roi==currentRoi) {
            Color lineColor = roi.getStrokeColor();
            Color fillColor = roi.getFillColor();
            float lineWidth = roi.getStrokeWidth();
            roi.setStrokeColor(null);
            roi.setFillColor(null);
            boolean strokeSet = roi.getStroke()!=null;
            if (strokeSet)
                roi.setStrokeWidth(1);
            roi.draw(g);
            roi.setStrokeColor(lineColor);
            if (strokeSet)
                roi.setStrokeWidth(lineWidth);
            roi.setFillColor(fillColor);
            currentRoi = null;
        } else
            roi.draw(g);
    }*/
    
    private void drawZoomIndicator(Graphics g) {
        int x1 = 10;
        int y1 = 10;
        double aspectRatio = (double)imageHeight/imageWidth;
        int w1 = 64;
        if (aspectRatio>1.0)
            w1 = (int)(w1/aspectRatio);
        int h1 = (int)(w1*aspectRatio);
        if (w1<4) w1 = 4;
        if (h1<4) h1 = 4;
        int w2 = (int)(w1*((double)srcRect.width/imageWidth));
        int h2 = (int)(h1*((double)srcRect.height/imageHeight));
        if (w2<1) w2 = 1;
        if (h2<1) h2 = 1;
        int x2 = (int)(w1*((double)srcRect.x/imageWidth));
        int y2 = (int)(h1*((double)srcRect.y/imageHeight));
        g.setColor(new Color(128, 128, 255));
        ((Graphics2D)g).setStroke(Roi.onePixelWide);
        g.drawRect(x1, y1, w1, h1);
        if (w2*h2<=200 || w2<10 || h2<10)
            g.fillRect(x1+x2, y1+y2, w2, h2);
        else
            g.drawRect(x1+x2, y1+y2, w2, h2);
    }
    
    /**************
     * here we resample an entire image at a reduced resolution
     * @param img
     * @param width
     * @param height
     * @param newwidth
     * @param newheight
     * @param scale
     * @return
     */
    public Image resampleImage(Image img,int width,int height,int newwidth,int newheight,float scale) {
    	ColorProcessor cp=new ColorProcessor(img);
    	int[] pix=(int[])cp.getPixels();
    	int[] newpix=new int[newwidth*newheight];
    	for(int i=0;i<newheight;i++) {
    		int ypos=(int)((float)i*scale);
    		if(ypos>=height) ypos=height-1;
    		for(int j=0;j<newwidth;j++) {
    			int xpos=(int)((float)j*scale);
        		if(xpos>=width) xpos=width-1;
        		newpix[j+i*newwidth]=pix[xpos+ypos*width];
    		}
    	}
    	return (new ColorProcessor(newwidth,newheight,newpix)).createImage();
    }
    
    /***********
     * here we resample a selected crop rectangle at a desired resolution
     * @param img
     * @param width
     * @param height
     * @param newwidth
     * @param newheight
     * @param scale
     * @param startx
     * @param starty
     * @return
     */
    public Image resampleImage(Image img,int width,int height,int newwidth,int newheight,float scale,int startx,int starty) {
    	ColorProcessor cp=new ColorProcessor(img);
    	int[] pix=(int[])cp.getPixels();
    	int[] newpix=new int[newwidth*newheight];
    	for(int i=0;i<newheight;i++) {
    		int ypos=starty+(int)((float)i*scale);
    		if(ypos>=height) ypos=height-1;
    		for(int j=0;j<newwidth;j++) {
    			int xpos=startx+(int)((float)j*scale);
        		if(xpos>=width) xpos=width-1;
        		newpix[j+i*newwidth]=pix[xpos+ypos*width];
    		}
    	}
    	return (new ColorProcessor(newwidth,newheight,newpix)).createImage();
    }

}
