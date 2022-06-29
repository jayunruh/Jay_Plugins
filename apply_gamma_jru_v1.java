import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;

public class apply_gamma_jru_v1 implements PlugIn, DialogListener {
	ImagePlus imp;
	ImageProcessor firstip;
	Object temppix;
	boolean canceled;
	float min,max,gamma;

	public void run(String arg) {
		imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int channel=imp.getChannel();
		firstip=imp.getProcessor();
		min=(float)firstip.getMin();
		max=(float)firstip.getMax();
		//IJ.log(""+min+" , "+max);
		gamma=1.0f;
		firstip.snapshot();
		showDialog(firstip);
		firstip.reset();
		if(canceled){
			imp.updateAndDraw();
			return;
		}
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				ImageProcessor ip=stack.getProcessor(channel+j*channels+i*channels*slices);
				apply_gamma(ip.getPixels());
			}
		}
		imp.updateAndDraw();
	}

	public void apply_gamma(Object image){
		if(image instanceof float[]){
			float[] temp=(float[])image;
			for(int i=0;i<temp.length;i++){
				float temp2=(float)Math.pow(((temp[i]-min)/(max-min)),gamma);
				temp[i]=temp2*(max-min)+min;
			}
		} else {
			if(image instanceof short[]){
				short[] temp=(short[])image;
				for(int i=0;i<temp.length;i++){
					float temp2=(float)(temp[i]&0xffff);
					float temp3=(float)Math.pow(((temp2-min)/(max-min)),gamma);
					float temp4=temp3*(max-min)+min;
					temp[i]=(short)((int)temp4);
				}
			} else {
				byte[] temp=(byte[])image;
				for(int i=0;i<temp.length;i++){
					float temp2=(float)(temp[i]&0xff);
					float temp3=(float)Math.pow(((temp2-min)/(max-min)),gamma);
					float temp4=temp3*(max-min)+min;
					temp[i]=(byte)((int)temp4);
				}
			}
		}
	}

	void showDialog(ImageProcessor ip) {
		GenericDialog gd = new GenericDialog("Gamma_Adjuster");
		gd.addNumericField("Gamma",gamma,5,15,null);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled())
			{canceled = true; return;}
		gamma=(float)gd.getNextNumber();
	}

	public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		firstip.reset();
		firstip.snapshot();
		gamma=(float)gd.getNextNumber();
		apply_gamma(firstip.getPixels());
		imp.updateAndDraw();
		return true;
	}

}
