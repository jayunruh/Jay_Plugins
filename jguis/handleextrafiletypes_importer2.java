package jguis;

import ij.IJ;
import ij.io.OpenDialog;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.net.URLDecoder;
import java.util.Iterator;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.TransferHandler;

public class handleextrafiletypes_importer2 extends JPanel implements ActionListener{
	// this class is an ugly hack that provides a drag and drop window from
	// which
	// fiji users can access the old "HandleExtraFileTypes" plugin
	private JButton openbutton;
	public boolean usecustomloci;

	public static Frame launch_frame(handleextrafiletypes_importer2 panel){
		final Frame f=new Frame("HEFT Importer");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				Component[] comps=f.getComponents();
				for(int i=0;i<comps.length;i++){
					comps[i].setVisible(false);
				}
				f.dispose();
			}
		});

		f.setLayout(null);
		panel.setBounds(10,40,220,400);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(new Dimension(250,450));
		f.setVisible(true);
		panel.requestFocus();
		return f;
	}

	public void init(){
		setLayout(null);
		openbutton=new JButton("Open/Drop");
		openbutton.setBounds(10,10,120,30);
		openbutton.addActionListener(this);
		add(openbutton);
		//DropTarget dp=new DropTarget(this,this);
		setTransferHandler(new ShortTransHandlerJ(this));
		IJ.register(this.getClass());
	}

	public void actionPerformed(ActionEvent e){
		if(e.getSource()==openbutton){
			open();
		}
	}

	public void open(){
		OpenDialog od=new OpenDialog("Open","");
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name!=null){
			String path=directory+name;
			openWithHEFT(path);
		}
	}

	public void openWithHEFT(String path){
		try{
			if(path.endsWith(".pw")||path.endsWith(".pw2")){
				Class heft=Class.forName("import plot jru v1");
				IJ.runPlugIn("import plot jru v1",path);
			}else{
				if(usecustomloci){
					Class customloci=Class.forName("LOCI_file_reader_jru_v1");
					IJ.runPlugIn("LOCI_file_reader_jru_v1",path);
				} else {
					Class heft=Class.forName("HandleExtraFileTypesjru");
					IJ.runPlugIn("HandleExtraFileTypesjru",path);
				}
			}
		}catch(ClassNotFoundException e){
			IJ.runPlugIn("HandleExtraFileTypes",path);
		}
		/*
		 * Object thePlugIn=null; try { Class c = Class.forName(className);
		 * thePlugIn = c.newInstance(); if (thePlugIn instanceof PlugIn)
		 * ((PlugIn)thePlugIn).run(arg); else new PlugInFilterRunner(thePlugIn,
		 * commandName, arg); } catch (ClassNotFoundException e) { if
		 * (IJ.getApplet()==null) log("Plugin or class not found: \"" +
		 * className + "\"\n(" + e+")"); } catch (InstantiationException e)
		 * {log("Unable to load plugin (ins)");} catch (IllegalAccessException
		 * e)
		 * {log("Unable to load plugin, possibly \nbecause it is not public.");}
		 */
	}


}

class ShortTransHandlerJ extends TransferHandler{
	
	protected handleextrafiletypes_importer2 spanel;
	
	public ShortTransHandlerJ(handleextrafiletypes_importer2 spanel){
		this.spanel=spanel;
	}
	
	public boolean canImport(JComponent comp,DataFlavor[] transferFlavors){
		return true;
	}
	
	public boolean importData(JComponent comp,Transferable t){
		//DataFlavor[] flavors=t.getTransferDataFlavors();
		try{
			if(t.isDataFlavorSupported(DataFlavor.javaFileListFlavor)){
				Object data=t.getTransferData(DataFlavor.javaFileListFlavor);
				Iterator iterator=((List)data).iterator();
				while(iterator.hasNext()){
					File file=(File)iterator.next();
					if(!file.isDirectory()) // don't handle directories for now
						spanel.openWithHEFT(file.getAbsolutePath());
				}
			} else if(t.isDataFlavorSupported(DataFlavor.stringFlavor)){
				String s=(String)t.getTransferData(DataFlavor.stringFlavor);
				IJ.log("transfer string = "+s);
				String[] fnames=s.split("[ \t\n\r\f]");
				for(int i=0;i<fnames.length;i++){
					//IJ.log(fnames[i]);
					fnames[i]=fnames[i].replaceAll("^file:/*", "/");
					if (!new File(fnames[i]).exists()) {
			            fnames[i] = URLDecoder.decode(fnames[i], "UTF-8");
			        }
					File file=new File(fnames[i]);
					if(!file.isDirectory()) // don't handle directories for now
						spanel.openWithHEFT(file.getAbsolutePath());
				}
			} else {return false;}
		}catch(Exception e){
			e.printStackTrace();
			IJ.error("Error in drop() method. Reason: "+e.getMessage());
			return false;
		}
		return true;
	}
}
