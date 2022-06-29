import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;

public class Time_Stamper_jru_v1 implements PlugIn, DialogListener {
	ImagePlus imp;
	ImageProcessor firstip;
	double time;
	static int x = 2;
	static int y = 15;
	static int size = 12;
	int maxWidth;
	Font font;
	static double start = 0;
	static double interval = 1;
	static String suffix = "sec";
	static int decimalPlaces = 0;
	boolean canceled;
	int frame, first, last,channels,slices,frames;
	Object temppix;

	public void run(String arg) {
		imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		channels=imp.getNChannels();
		if(frames==1){
			frames=slices*channels;
			slices=1;
			channels=1;
		}
		imp.setPosition(1,1,1);
		firstip=imp.getProcessor();
		firstip.snapshot();
		first=1;
		last=frames;
		showDialog(firstip);
		firstip.reset();
		if(canceled){
			imp.updateAndDraw();
			return;
		}
		time=start+(first-1)*interval;
		for(int i=(first-1);i<last;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					ImageProcessor ip=stack.getProcessor(1+k+j*channels+i*channels*slices);
					ip.setAntialiasedText(true);
					ip.setFont(font);
					ip.setColor(Toolbar.getForegroundColor());
					String s = getString(time);
					ip.moveTo(x+maxWidth-ip.getStringWidth(s), y);
					ip.drawString(s);
				}
			}
			time+=interval;
		}
		imp.updateAndDraw();
	}
	
	String getString(double time) {
		if (interval==0.0)
			return suffix;
		else
			return (decimalPlaces==0?""+(int)time:IJ.d2s(time, decimalPlaces))+" "+suffix;
	}

	void showDialog(ImageProcessor ip) {
		Rectangle roi = ip.getRoi();
		if (roi.width<ip.getWidth() || roi.height<ip.getHeight()) {
			x = roi.x;
			y = roi.y+roi.height;
			size = (int) ((roi.height - 1.10526)/0.934211);	
			if (size<7) size = 7;
			if (size>80) size = 80;
		}
		GenericDialog gd = new GenericDialog("Time Stamper");
		gd.addNumericField("Starting Time:", start, 2);
		gd.addNumericField("Time Btw Frames:", interval, 2);
		gd.addNumericField("X Location:", x, 0);
		gd.addNumericField("Y Location:", y, 0);
		gd.addNumericField("Font Size:", size, 0);
		gd.addNumericField("Decimal Places:", decimalPlaces, 0);
		gd.addNumericField("First Frame:", first, 0);
		gd.addNumericField("Last Frame:", last, 0);
		gd.addStringField("Suffix:", suffix);
		gd.addCheckbox("Bold",false);
		gd.addCheckbox("Only_Suffix",false);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled())
			{canceled = true; return;}
		start = gd.getNextNumber();
 		interval = gd.getNextNumber();
		x = (int)gd.getNextNumber();
		y = (int)gd.getNextNumber();
		size = (int)gd.getNextNumber();
		decimalPlaces = (int)gd.getNextNumber();
		first = (int)gd.getNextNumber();
		last = (int)gd.getNextNumber();
		suffix = gd.getNextString();
		boolean bold=gd.getNextBoolean();
		if(gd.getNextBoolean()){interval=0.0;}
		if(bold){
			font = new Font("SansSerif", Font.BOLD, size);
		} else {
			font = new Font("SansSerif", Font.PLAIN, size);
		}
		ip.setFont(font);
		time = start;
		maxWidth = ip.getStringWidth(getString(start+interval*frames));
	}

	public boolean dialogItemChanged(GenericDialog gd,AWTEvent e){
		firstip.reset();
		firstip.snapshot();
		start=gd.getNextNumber();
		interval=gd.getNextNumber();
		x=(int)gd.getNextNumber();
		y=(int)gd.getNextNumber();
		size=(int)gd.getNextNumber();
		decimalPlaces=(int)gd.getNextNumber();
		first=(int)gd.getNextNumber();
		last=(int)gd.getNextNumber();
		suffix=gd.getNextString();
		boolean bold=gd.getNextBoolean();
		if(gd.getNextBoolean()){interval=0.0;}
		if(bold){
			font = new Font("SansSerif", Font.BOLD, size);
		} else {
			font = new Font("SansSerif", Font.PLAIN, size);
		}
		firstip.setFont(font);
		time = start;
		maxWidth = firstip.getStringWidth(getString(start+interval*frames));
		firstip.setAntialiasedText(true);
		firstip.setFont(font);
		firstip.setColor(Toolbar.getForegroundColor());
		String s = getString(time);
		firstip.moveTo(x+maxWidth-firstip.getStringWidth(s), y);
		firstip.drawString(s);
		imp.updateAndDraw();
		//IJ.log(s);
		return true;
	}

}
