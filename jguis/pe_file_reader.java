/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import loci.common.ByteArrayHandle;
import loci.common.DataTools;
import loci.common.IRandomAccess;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.common.services.DependencyException;
import loci.common.services.ServiceFactory;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MissingLibraryException;
import loci.formats.codec.LZOCodec;
import loci.formats.services.MetakitService;
import ome.metakit.MetakitException;
import ome.metakit.MetakitReader;

public class pe_file_reader extends FormatReader{
	
	  private static final String DATA_DIR = "Data";
	  private static final String EMBEDDED_STREAM = "embedded-stream.raw";

	  private static final int SIGNATURE_SIZE = 13;

	  // -- Fields --

	  private ArrayList<Stack> stacks;
	  private ArrayList<String> extraFiles;
	  private Object[][] sampleTable, stringTable;
	  private Location dir = null;

	  private ArrayList<Double[]> timestamps = new ArrayList<Double[]>();

	public void open(String dir,String fname){
		try{
			MetakitReader mkr=new MetakitReader(dir+fname);
			int tablecount=mkr.getTableCount();
			String[] tablenames=mkr.getTableNames();
			for(int i=0;i<tablecount;i++){
				IJ.log(tablenames[i]);
				Object[][] tdata=mkr.getTableData(i);
				//data are stored in row arrays
				if(tdata==null || tdata.length<1) continue;
				IJ.log("obj dimensions = "+tdata.length+" , "+tdata[0].length);
				List<List<String>> listtable=new ArrayList<List<String>>();
				//int nrows=mkr.getRowCount(i);
				//IJ.log("nrows = "+nrows);
				String[] colnames=mkr.getColumnNames(i);
				Class[] coltypes=mkr.getColumnTypes(i);
				for(int j=0;j<colnames.length;j++){
					IJ.log(colnames[j]+" , "+coltypes[j].getName());
				}
				for(int j=0;j<tdata.length;j++){
					if(tdata[j]==null) continue;
					List<String> row=new ArrayList<String>();
					for(int k=0;k<tdata[j].length;k++){
						if(tdata[j][k]==null){
							row.add("NULL");
							continue;
						}
						if(tdata[j][k] instanceof Integer){
							row.add(((Integer)tdata[j][k]).toString());
						} else if(tdata[j][k] instanceof String) {
							row.add(((String)tdata[j][k]));
						} else {
							row.add(tdata[j][k].toString());
						}
					}
					listtable.add(row);
				}
				if(listtable.size()>0) table_tools.create_table(tablenames[i],listtable,colnames);
			}
			mkr.close();
		}catch(IOException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}catch(MetakitException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	  // -- Constructor --

	  /** Constructs a new Volocity reader. */
	  public pe_file_reader() {
	    //super("pe_file_reader",new String[] {"mvd2", "aisf", "aiix", "dat", "atsf"});
		  super("pe_file_reader","");
	    domains = new String[] {FormatTools.UNKNOWN_DOMAIN};
	    hasCompanionFiles = true;
	    datasetDescription = "One .mvd2 file plus a 'Data' directory";
	  }

	  // -- IFormatReader API methods --

	  /* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
	  @Override
	  public String[] getSeriesUsedFiles(boolean noPixels) {
	    FormatTools.assertId(currentId, true, 1);

	    ArrayList<String> files = new ArrayList<String>();
	    files.addAll(extraFiles);
	    if (!noPixels) {
	      Stack stack = stacks.get(getSeries());
	      for (int c=0; c<getEffectiveSizeC(); c++) {
	        files.add(stack.pixelsFiles[c]);
	      }
	      if (stack.timestampFile != null) {
	        files.add(stack.timestampFile);
	      }
	    }
	    return files.toArray(new String[files.size()]);
	  }

	  /* @see loci.formats.IFormatReader#isThisType(String, boolean) */
	  @Override
	  public boolean isThisType(String name, boolean open) {
	    if (checkSuffix(name, "mvd2")) {
	      return super.isThisType(name, open);
	    }

	    if (open && checkSuffix(name, suffixes)) {
	      Location file = new Location(name).getAbsoluteFile();
	      Location parent = file.getParentFile();
	      parent = parent.getParentFile();
	      if (parent != null) {
	        parent = parent.getParentFile();
	        if (parent != null) {
	          Location mvd2 = new Location(parent, parent.getName() + ".mvd2");
	          return mvd2.exists() && super.isThisType(mvd2.getAbsolutePath());
	        }
	      }
	    }
	    return false;
	  }

	  /* @see loci.formats.IFormatReader#isThisType(RandomAccessInputStream) */
	  @Override
	  public boolean isThisType(RandomAccessInputStream stream) throws IOException {
	    final int blockLen = 2;
	    if (!FormatTools.validStream(stream, blockLen, false)) return false;
	    String check = stream.readString(blockLen);
	    return check.equals("JL") || check.equals("LJ");
	  }

	  /**
	   * @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
	   */
	  @Override
	  public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
	    throws FormatException, IOException
	  {
	    FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);

	    int[] zct = getZCTCoords(no);

	    Stack stack = stacks.get(getSeries());

	    if (!new Location(stack.pixelsFiles[zct[1]]).exists()) {
	      Arrays.fill(buf, (byte) 0);
	      return buf;
	    }

	    RandomAccessInputStream pix =
	      new RandomAccessInputStream(stack.pixelsFiles[zct[1]]);

	    int padding = zct[2] * stack.planePadding;

	    long planeSize = FormatTools.getPlaneSize(this);
	    int planesInFile = (int) (pix.length() / planeSize);
	    int planeIndex = no / getEffectiveSizeC();
	    if (planesInFile == getSizeT()) {
	      planeIndex = zct[2];

	      int block = stack.blockSize;
	      padding = block - (int) (planeSize % block);
	      if (padding == block) {
	        padding = 0;
	      }
	      padding *= zct[2];
	    }

	    long offset = (long) stack.blockSize + planeIndex * planeSize + padding;
	    if (offset >= pix.length()) {
	      pix.close();
	      return buf;
	    }
	    pix.seek(offset);

	    if (stack.clippingData) {
	      pix.seek(offset - 3);
	      ByteArrayHandle v = new ByteArrayHandle();
	      while (v.length() < FormatTools.getPlaneSize(this) &&
	        pix.getFilePointer() < pix.length())
	      {
	        try {
	          byte[] b = new LZOCodec().decompress(pix, null);
	          pix.skipBytes(4);
	          v.write(b);
	        }
	        catch (IOException e) { }
	      }
	      RandomAccessInputStream s = new RandomAccessInputStream(v);
	      s.seek(0);
	      readPlane(s, x, y, w, h, buf);
	      s.close();
	    }
	    else {
	      if (pix.getFilePointer() + planeSize > pix.length()) {
	        pix.close();
	        return buf;
	      }
	      readPlane(pix, x, y, w, h, buf);
	    }
	    pix.close();

	    if (getRGBChannelCount() == 4) {
	      // stored as ARGB, need to swap to RGBA
	      for (int i=0; i<buf.length/4; i++) {
	        byte a = buf[i * 4];
	        buf[i * 4] = buf[i * 4 + 1];
	        buf[i * 4 + 1] = buf[i * 4 + 2];
	        buf[i * 4 + 2] = buf[i * 4 + 3];
	        buf[i * 4 + 3] = a;
	      }
	    }

	    return buf;
	  }

	  /* @see loci.formats.IFormatReader#close(boolean) */
	  @Override
	  public void close(boolean fileOnly) throws IOException {
	    super.close(fileOnly);
	    if (!fileOnly) {
	      stacks = null;
	      sampleTable = null;
	      stringTable = null;
	      dir = null;
	      timestamps.clear();
	      Location.mapFile(EMBEDDED_STREAM, null);
	    }
	  }

	  // -- Internal FormatReader API methods --

	  /* @see loci.formats.FormatReader#initFile(String) */
	  public void initFile(String id) throws FormatException, IOException {
		    if (!checkSuffix(id, "mvd2")) {
		      Location file = new Location(id).getAbsoluteFile();
		      Location parent = file.getParentFile().getParentFile();
		      String[] files = parent.list(true);
		      for (String f : files) {
		        if (checkSuffix(f, "mvd2")) {
		          id = new Location(parent, f).getAbsolutePath();
		          break;
		        }
		      }
		    }

		    super.initFile(id);

		    stacks = new ArrayList<Stack>();
		    extraFiles = new ArrayList<String>();

		    Location file = new Location(id).getAbsoluteFile();
		    extraFiles.add(file.getAbsolutePath());

		    Location parentDir = file.getParentFile();
		    dir = new Location(parentDir, DATA_DIR);

		    if (dir.exists()) {
		      String[] files = dir.list(true);
		      for (String f : files) {
		        if (!checkSuffix(f, "aisf") && !checkSuffix(f, "atsf")) {
		          extraFiles.add(new Location(dir, f).getAbsolutePath());
		        }
		      }
		    }

		    try {
		      ServiceFactory factory = new ServiceFactory();
		      MetakitService reader = factory.getInstance(MetakitService.class);
		      reader.initialize(id);
		      sampleTable = reader.getTableData(1);
		      stringTable = reader.getTableData(2);

		      reader.close();
		    }
		    catch (DependencyException e) {
		      throw new MissingLibraryException("Could not find Metakit library", e);
		    }

		    ArrayList<String> stackNames = new ArrayList<String>();
		    ArrayList<Integer> parentIDs = new ArrayList<Integer>();

		    for (int i=0; i<sampleTable.length; i++) {
		      Integer stringID = (Integer) sampleTable[i][11];
		      String name = getString(stringID);

		      int channelIndex = getChildIndex((Integer) sampleTable[i][0], "Channels");

		      if (i > 0 && (Integer) sampleTable[i][2] == 1 && (channelIndex >= 0 ||
		        (sampleTable[i][14] != null && !sampleTable[i][14].equals(0)) ||
		        ((byte[]) sampleTable[i][13]).length > 21))
		      {
		        if (channelIndex < 0) {
		          RandomAccessInputStream s = getStream(i);
		          s.seek(0);
		          if (s.read() != 'I') {
		            s.order(false);
		          }
		          s.seek(22);
		          int x = s.readInt();
		          int y = s.readInt();
		          int z = s.readInt();
		          if (x * y * z > 0 && x * y * z < (s.length() * 3)) {
		            stackNames.add(name);
		            parentIDs.add((Integer) sampleTable[i][0]);
		          }
		          s.close();
		        }
		        else {
		          stackNames.add(name);
		          parentIDs.add((Integer) sampleTable[i][0]);
		        }
		      }
		    }

		    for (int i=0; i<parentIDs.size(); i++) {
		      Stack stack = new Stack();
		      stack.core = new CoreMetadata();
		      Integer parent = parentIDs.get(i);

		      int channelIndex = getChildIndex(parent, "Channels");
		      if (channelIndex >= 0) {
		        Integer[] channels =
		          getAllChildren((Integer) sampleTable[channelIndex][0]);
		        stack.core.sizeC = channels.length;
		        stack.pixelsFiles = new String[stack.core.sizeC];

		        stack.channelNames = new String[channels.length];
		        for (int c=0; c<channels.length; c++) {
		          stack.channelNames[c] =
		            getString((Integer) sampleTable[channels[c]][11]);

		          RandomAccessInputStream data = getStream(channels[c]);
		          if (data.length() > 22) {
		            data.seek(22);
		            int stackID = data.readInt();
		            Location f = new Location(dir, stackID + ".aisf");
		            if (!f.exists()) {
		              f = new Location(dir, DataTools.swap(stackID) + ".aisf");
		            }
		            stack.pixelsFiles[c] = f.getAbsolutePath();
		          }
		          else {
		            Integer child =
		              getAllChildren((Integer) sampleTable[channels[c]][0])[0];
		            stack.pixelsFiles[c] =
		              getFile((Integer) sampleTable[child][0], dir);
		          }
		          data.close();
		        }
		      }
		      else {
		        stack.pixelsFiles = new String[1];
		        stack.pixelsFiles[0] = getFile(parent, dir);

		        if (stack.pixelsFiles[0] == null ||
		          !new Location(stack.pixelsFiles[0]).exists())
		        {
		          int row = -1;
		          for (int r=0; r<sampleTable.length; r++) {
		            if (sampleTable[r][0].equals(parent)) {
		              row = r;
		              break;
		            }
		          }

		          stack.pixelsFiles[0] = EMBEDDED_STREAM;
		          IRandomAccess data =
		            new ByteArrayHandle((byte[]) sampleTable[row][13]);
		          Location.mapFile(stack.pixelsFiles[0], data);
		        }
		      }

		      RandomAccessInputStream data = null;

		      int timestampIndex = getChildIndex(parent, "Timepoint times stream");
		      if (timestampIndex >= 0) {
		        data = getStream(timestampIndex);
		        data.seek(22);
		        int timestampID = data.readInt();
		        Location f = new Location(dir, timestampID + ".atsf");
		        if (!f.exists()) {
		          f = new Location(dir, DataTools.swap(timestampID) + ".atsf");
		        }
		        stack.timestampFile = f.getAbsolutePath();
		        data.close();
		      }

		      int xIndex = getChildIndex(parent, "um/pixel (X)");
		      if (xIndex >= 0) {
		        data = getStream(xIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.physicalX = data.readDouble();
		        data.close();
		      }

		      int yIndex = getChildIndex(parent, "um/pixel (Y)");
		      if (yIndex >= 0) {
		        data = getStream(yIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.physicalY = data.readDouble();
		        data.close();
		      }

		      int zIndex = getChildIndex(parent, "um/pixel (Z)");
		      if (zIndex >= 0) {
		        data = getStream(zIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.physicalZ = data.readDouble();
		        data.close();
		      }

		      timestampIndex = getChildIndex(parent, "TimepointTimes");
		      if (timestampIndex >= 0) {
		        data = getStream(timestampIndex);
		        data.seek(SIGNATURE_SIZE);
		        data.close();
		      }

		      int objectiveIndex = getChildIndex(parent, "Microscope Objective");
		      if (objectiveIndex >= 0) {
		        data = getStream(objectiveIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.magnification = data.readDouble();
		        data.close();
		      }

		      int detectorIndex = getChildIndex(parent, "Camera/Detector");
		      if (detectorIndex >= 0) {
		        data = getStream(detectorIndex);
		        data.seek(SIGNATURE_SIZE);
		        int len = data.readInt();
		        stack.detector = data.readString(len);
		        data.close();
		      }

		      int descriptionIndex = getChildIndex(parent, "Experiment Description");
		      if (descriptionIndex >= 0) {
		        data = getStream(descriptionIndex);
		        data.seek(SIGNATURE_SIZE);
		        int len = data.readInt();
		        stack.description = data.readString(len);
		        data.close();
		      }
		      
		      int exptlogIndex = getChildIndex(parent, "Experiment Event Log");
		      if (exptlogIndex >= 0) {
		        data = getStream(exptlogIndex);
		        data.seek(SIGNATURE_SIZE);
		        int len = data.readInt();
		        stack.eventLog = data.readString(len);
		        data.close();
		        IJ.log(stack.eventLog);
		      }

		      int xLocationIndex = getChildIndex(parent, "X Location");
		      if (xLocationIndex >= 0) {
		        data = getStream(xLocationIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.xLocation = data.readDouble();
		        data.close();
		      }

		      int yLocationIndex = getChildIndex(parent, "Y Location");
		      if (yLocationIndex >= 0) {
		        data = getStream(yLocationIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.yLocation = data.readDouble();
		        data.close();
		      }

		      int zLocationIndex = getChildIndex(parent, "Z Location");
		      if (zLocationIndex >= 0) {
		        data = getStream(zLocationIndex);
		        data.seek(SIGNATURE_SIZE);
		        stack.zLocation = data.readDouble();
		        data.close();
		      }

		      stacks.add(stack);
		    }

	}
	  
	  private String getString(Integer stringID) {
		    for (int row=0; row<stringTable.length; row++) {
		      if (stringID.equals(stringTable[row][0])) {
		        String s = (String) stringTable[row][1];
		        if (s != null) {
		          s = s.trim();
		        }
		        return s;
		      }
		    }
		    return null;
		  }

		  private int getChildIndex(Integer parentID, String childName) {
		    for (int row=0; row<sampleTable.length; row++) {
		      if (parentID.equals(sampleTable[row][1])) {
		        String name = getString((Integer) sampleTable[row][11]);
		        if (childName.equals(name)) {
		          return row;
		        }
		      }
		    }
		    return -1;
		  }

		  private Integer[] getAllChildren(Integer parentID) {
		    ArrayList<Integer> children = new ArrayList<Integer>();
		    for (int row=0; row<sampleTable.length; row++) {
		      if (parentID.equals(sampleTable[row][1])) {
		        children.add(row);
		      }
		    }
		    return children.toArray(new Integer[children.size()]);
		  }

		  private RandomAccessInputStream getStream(int row) throws IOException {
		    Object o = sampleTable[row][14];
		    String fileLink = o == null ? "0" : o.toString().trim();
		    RandomAccessInputStream data = null;
		    if (fileLink.equals("0")) {
		      data = new RandomAccessInputStream((byte[]) sampleTable[row][13]);
		    }
		    else {
		      fileLink = new Location(dir, fileLink + ".dat").getAbsolutePath();
		      data = new RandomAccessInputStream(fileLink);
		    }

		    data.order(true);
		    return data;
		  }

		  private String getFile(Integer parent, Location dir) {
		    for (int row=0; row<sampleTable.length; row++) {
		      if (parent.equals(sampleTable[row][0])) {
		        Object o = sampleTable[row][14];
		        if (o != null) {
		          String fileLink = o.toString().trim() + ".dat";
		          return new Location(dir, fileLink).getAbsolutePath();
		        }
		      }
		    }
		    return null;
		  }

		  // -- Helper class --

		  class Stack {
		    public String[] pixelsFiles;
		    public String timestampFile;
		    public int planePadding;
		    public int blockSize;
		    public boolean clippingData;

		    public CoreMetadata core;
		    public String[] channelNames;
		    public Double physicalX;
		    public Double physicalY;
		    public Double physicalZ;
		    public Double magnification;
		    public String detector;
		    public String description;
		    public String eventLog;
		    public double xLocation;
		    public double yLocation;
		    public double zLocation;
		  }

}
