package jalgs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class file_walker2 extends SimpleFileVisitor<Path>{
	//this version uses the java file walker functionality
	//long totsize=0;
	//need to make a hash table holding parents, sizes, file numbers, creationdates
	HashMap<Path,Long[]> dirmap;
	String[] keepexts;
	int maxdocdepth;
	
	public file_walker2(String[] keepexts,int maxdocdepth) {
		dirmap=new HashMap<>();
		this.keepexts=keepexts;
		//keepexts=new String[]{".xlsx",".pdf",".docx",".html",".pptx",".yaml",".yml"};
		this.maxdocdepth=maxdocdepth;
	}
	
	public FileVisitResult visitFile(Path file, BasicFileAttributes attr) {
		if(attr.isSymbolicLink()) {
			//do nothing for links
		} else if (attr.isRegularFile()) {
			//System.out.println(file.toString()+","+attr.size());
			Path parent=file.getParent();
			if(dirmap.containsKey(parent)) {
				Long[] temp=dirmap.get(parent);
				temp[0]+=attr.size();
				temp[1]++;
				temp[2]=Math.max(temp[2], attr.creationTime().to(TimeUnit.SECONDS));
				dirmap.put(parent,temp);
			} else {
				dirmap.put(parent,new Long[] {attr.size(),1L,attr.creationTime().to(TimeUnit.SECONDS)});
				//System.out.println("missing key");
			}
			if(file.getNameCount()<=maxdocdepth) {
				String tname=file.toString().toLowerCase();
				for(int i=0;i<keepexts.length;i++) {
					if(tname.endsWith(keepexts[i])) {
						dirmap.put(file,new Long[] {-1L,-1L,attr.creationTime().to(TimeUnit.SECONDS)});
					}
				}
			}
			//totsize+=attr.size();
		} else {
			//not sure what falls in this category
		}
		return FileVisitResult.CONTINUE;
	}
	
	public FileVisitResult preVisitDirectory(Path dir,IOException e) {
		//System.out.println(dir.toString()+",");
		//dirmap.put(dir,0L);
		return FileVisitResult.CONTINUE;
	}
	
	public FileVisitResult postVisitDirectory(Path dir,IOException e) {
		//System.out.println(dir.toString()+","+((float)dirmap.get(dir))/1048576.0f);
		//dirmap.remove(dir);
		return FileVisitResult.CONTINUE;
	}
	
	public FileVisitResult visitFileFailed(Path file,IOException e) {
		System.err.println(e);
		return FileVisitResult.CONTINUE;
	}
	
	public int[] done(String outdirname,String outdocname) {
		StringBuffer dirsb=new StringBuffer();
		//System.out.println("Directories:name,size(MB),nfiles,creation");
		//dirsb.append("name,size,nfiles,creationdate\n");
		int ndirs=0;
		for(Path key:dirmap.keySet()) {
			Long[] temp=dirmap.get(key);
			if(temp[0]>=0L) {
				//this is a directory
				String outstr=key.toString()+","+((float)temp[0])/1048576.0f+","+temp[1]+","+temp[2];
				//System.out.println(outstr);
				dirsb.append(outstr+"\n");
				ndirs++;
			}
		}
		if(outdirname.equals("sys")) {
			System.out.print(dirsb.toString());
		} else {
			try {
				File outdirfile=new File(outdirname);
				BufferedWriter bw=new BufferedWriter(new FileWriter(outdirfile));
				bw.write(dirsb.toString());
				bw.close();
			} catch(IOException e) {
				System.out.println("Error writing file");
				System.err.println(e.toString());
			}
		}
		//System.out.println("Wrote "+ndirs+" directories");
		dirsb=null;
		//System.out.println("Document list:name,creation");
		int ndocs=0;
		StringBuffer docsb=new StringBuffer();
		//docsb.append("name,creation\n");
		for(Path key:dirmap.keySet()) {
			Long[] temp=dirmap.get(key);
			if(temp[0]<0L) {
				//this is a document
				String outstr=key.toString()+","+temp[2];
				//System.out.println(outstr+"\n");
				docsb.append(outstr+"\n");
				ndocs++;
			}
		}
		if(outdocname.equals("sys")) {
			System.out.print(docsb.toString());
		} else {
			try {
				File outdocfile=new File(outdocname);
				BufferedWriter bw=new BufferedWriter(new FileWriter(outdocfile));
				bw.write(docsb.toString());
				bw.close();
			} catch(IOException e) {
				System.out.println("Error writing file");
				System.err.println(e.toString());
			}
		}
		//System.out.println("Wrote "+ndocs+" documents");
		return new int[] {ndirs,ndocs};
	}
	
	public static void main(String[] args) throws IOException {
		if(args[0].equals("--help")){
			System.out.println("usage: file_walker2 indir maxreldepth");
			return;
		}
		Path startdir=Paths.get(args[0]);
		int relmaxdepth=Integer.parseInt(args[1]);
		int startdepth=(startdir).getNameCount();
		String[] keepexts={".xlsx",".pdf",".docx",".html",".pptx",".yaml",".yml"};
		/*file_walker2 fw=new file_walker2(keepexts,startdepth+relmaxdepth);
		Files.walkFileTree(startdir,fw);
		fw.done(args[1],args[2]);*/
		//run multithreaded file walkers over top level folders
		File[] flist=(new File(startdir.toString())).listFiles();
		int nthreads=10;
		ExecutorService executor=Executors.newFixedThreadPool(10);
		for(int i=0;i<flist.length;i++) {
			if(flist[i].isDirectory()) {
				Runnable worker=new DirWalker(flist[i].getPath(),startdepth+relmaxdepth-1,keepexts);
				executor.execute(worker);
			}
		}
		executor.shutdown();
		return;
	}
}

class DirWalker implements Runnable{
	//here is a class to run a directory in a thread
	String startdir;
	file_walker2 fw;
	int[] nfiles;
	
	public DirWalker(String startdir,int maxdepth,String[] keepexts) {
		this.startdir=startdir;
		fw=new file_walker2(keepexts,maxdepth);
	}
	
	public void run() {
		try {
			Files.walkFileTree(Paths.get(startdir),fw);
		} catch(IOException e) {
			System.err.println(e.toString());
		}
		nfiles=fw.done("sys", "sys");
	}
}
