package j3D;

import java.awt.Color;
import java.awt.Dimension;
import java.io.ByteArrayOutputStream;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.emf.EMFGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.freehep.graphicsio.ps.PSGraphics2D;

public class RenderVector{
	public ByteArrayOutputStream os;
	public VectorGraphics g;
	
	public RenderVector(int width,int height,Color background,int type){
		//types are EMF, PS (EPS), and PDF
		try{
			os=new ByteArrayOutputStream();
			if(type==0){
				g=new EMFGraphics2D(os,new Dimension(width,height));	
			} else if(type==1){
				g=new PSGraphics2D(os,new Dimension(width,height));
			} else {
				g=new PDFGraphics2D(os,new Dimension(width,height));
			}
			g.setDeviceIndependent(true);
			g.startExport();
			Color tempcolor=g.getColor();
			g.setColor(background);
			g.fillRect(0,0,width,height);
			g.setColor(tempcolor);
		}catch(Throwable e){
			g=null;
		}
	}
	
	public void endExport(){
		if(g!=null) g.endExport();
		g=null;
	}
	
	public void drawElement(element3D element){
		if(g!=null) element.drawelement(g);
	}
	
	public byte[] getByteArray(){
		try{
			if(g!=null) return os.toByteArray();
		} catch(Throwable e){
			return null;
		}
		return null;
	}

}
