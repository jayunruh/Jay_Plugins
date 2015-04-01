package jguis;

import ij.IJ;
import ij.io.OpenDialog;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Iterator;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;

public class handleextrafiletypes_importer extends JPanel implements ActionListener,DropTargetListener{
	// this class is an ugly hack that provides a drag and drop window from
	// which
	// fiji users can access the old "HandleExtraFileTypes" plugin
	private JButton openbutton;
	public boolean usecustomloci;

	public static Frame launch_frame(handleextrafiletypes_importer panel){
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
		DropTarget dp=new DropTarget(this,this);
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

	public void dragEnter(DropTargetDragEvent arg0){
		// check if there is at least one File Type in the list
		DataFlavor[] flavors=arg0.getCurrentDataFlavors();
		for(int i=0;i<flavors.length;i++){
			String mimeType=flavors[i].getMimeType();
			// System.out.println(flavors[i].toString());
			// mimeType= application/x-java-file-list
			if(flavors[i].isFlavorJavaFileListType()){
				arg0.acceptDrag(DnDConstants.ACTION_COPY);
				return;
			}
		}
		arg0.rejectDrag();
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

	public void dragExit(DropTargetEvent arg0){
	}

	public void dragOver(DropTargetDragEvent arg0){
	}

	public void drop(DropTargetDropEvent arg0){
		arg0.acceptDrop(DnDConstants.ACTION_COPY);
		try{
			Transferable t=arg0.getTransferable();
			if(t.isDataFlavorSupported(DataFlavor.javaFileListFlavor)){
				Object data=t.getTransferData(DataFlavor.javaFileListFlavor);
				Iterator iterator=((List)data).iterator();
				while(iterator.hasNext()){
					File file=(File)iterator.next();
					if(!file.isDirectory()) // don't handle directories for now
						openWithHEFT(file.getAbsolutePath());
				}
			}
		}catch(Exception e){
			e.printStackTrace();
			IJ.error("Error in drop() method. Reason: "+e.getMessage());
			arg0.dropComplete(false);
			return;
		}

		arg0.dropComplete(true);
	}

	public void dropActionChanged(DropTargetDragEvent arg0){
	}

}
