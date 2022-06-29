package jalgs;

import java.io.File;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Random;

public class file_walker {
	//this class has methods to walk the file paths multithreaded and return lists of files and sizes
	//if a folder has more than maxdirsize files, it is assumed to contain no directories
	//and it's size is estimated from 25 randomly selected files
	
	public static void main(String[] args) {
		if(args[0].contains("help")) {
			System.out.println("Usage: java file_walker dir [maxdirsize(1000) nthreads(10) docexts(commalist)]");
			System.out.println("Output: path,files_size(MB),epoch_date_modified,nfiles");
			return;
		}
		String[] keepexts= {".xlsx",".pdf",".docx",".html",".pptx",".yaml",".yml"};
		int maxdirsize=1000;
		int nthreads=10;
		if(args.length>1) maxdirsize=Integer.parseInt(args[1]);
		if(args.length>2) nthreads=Integer.parseInt(args[2]);
		if(args.length>3) keepexts=args[3].split(",");
		walkFolder(args[0],keepexts,nthreads,maxdirsize);
	}
	
	public static void walkFolder(String parentdir,String[] keepexts,int nthreads,int maxdirsize) {
		ExecutorService executor=Executors.newFixedThreadPool(nthreads);
		File parent=new File(parentdir);
		File[] firstlist=parent.listFiles();
		//List<Object[]> biglist=new ArrayList<Object[]>();
		int firstcount=0;
		long firstctime=0;
		long firstsize=0;
		for(int i=0;i<firstlist.length;i++) {
			if(firstlist[i].isDirectory()) {
				File[] secondlist=firstlist[i].listFiles();
				int secondcount=0;
				long secondctime=0;
				long secondsize=0;
				for(int j=0;j<secondlist.length;j++) {
					if(secondlist[j].isDirectory()) {
						//at this level start a thread
						Runnable worker=new DirectoryWalker(secondlist[j],keepexts,maxdirsize);
						executor.execute(worker);
					} else {
						secondcount++;
						secondsize+=secondlist[j].length();
						long tempctime=secondlist[j].lastModified();
						secondctime=Math.max(secondctime, tempctime);
						String tname=secondlist[j].getName().toLowerCase();
						if(hasKeepExt(tname,keepexts)) System.out.println(secondlist[j].getPath()+",-1,"+tempctime+",-1");
					}
				}
				System.out.println("\""+firstlist[i].getPath()+"\","+(float)secondsize/1048576.0f+","+secondctime/1000L+","+secondcount);
			} else {
				firstcount++;
				firstsize+=firstlist[i].length();
				long tempctime=firstlist[i].lastModified();
				firstctime=Math.max(firstctime, tempctime);
				String tname=firstlist[i].getName().toLowerCase();
				if(hasKeepExt(tname,keepexts)) System.out.println(firstlist[i].getPath()+",-1,"+tempctime+",-1");
			}
			System.out.println("\""+parent.getPath()+"\","+(float)firstsize/1048576.0f+","+firstctime/1000L+","+firstcount);
		}
		executor.shutdown();
		while(!executor.isTerminated()) {}
	}
	
	public static boolean hasKeepExt(String name,String[] keepexts) {
		for(int j=0;j<keepexts.length;j++) {
			if(name.endsWith(keepexts[j])) {
				return true;
			}
		}
		return false;
	}
	
}

class DirectoryWalker implements Runnable{
	File parent;
	String[] keepexts;
	long size;
	int fcount;
	int maxdirsize;
	
	public DirectoryWalker(File parent,String[] keepexts,int maxdirsize) {
		this.parent=parent;
		this.keepexts=keepexts;
		this.size=0L;
		this.fcount=0;
		this.maxdirsize=maxdirsize;
	}
	
	public void run() {
		//run through this directory, outputting names of sub-directories and document files
		//if the directory has >maxdirsize files in it ignore it and it's documents and estimate it's size from the first 20
		File[] flist=parent.listFiles();
		if(flist==null) return;
		long mtime=0;
		if(flist.length<maxdirsize) {
			for(int i=0;i<flist.length;i++) {
				if(flist[i].isDirectory()) {
					Runnable worker=new DirectoryWalker(flist[i],keepexts,maxdirsize);
					worker.run();
					//executor.execute(worker);
				} else {
					fcount++;
					size+=flist[i].length();
					long tempmtime=flist[i].lastModified();
					mtime=Math.max(mtime, tempmtime);
					String tname=flist[i].getName().toLowerCase();
					if(file_walker.hasKeepExt(tname,keepexts)) {
						System.out.println("\""+flist[i].getPath()+"\",-1,"+tempmtime/1000L+",-1");
						break;
					}
				}
			}
		} else {
			long sumsize=0;
			fcount+=flist.length;
			int tcount=0;
			//would be better to randomly sample here rather than use first 25
			Random jrnd=new Random();
			for(int i=0;i<25;i++) {
				int idx=(int)(jrnd.nextDouble()*(double)flist.length);
				//these are probably not directories but check anyway just in case
				if(flist[idx].isDirectory()) {
					Runnable worker=new DirectoryWalker(flist[idx],keepexts,maxdirsize);
					worker.run();
					//executor.execute(worker);
				} else {
					sumsize+=flist[idx].length();
					mtime=Math.max(mtime, flist[idx].lastModified());
					tcount++;
				}
			}
			size=(long)(((float)sumsize*(float)flist.length)/(float)tcount);
		}
		System.out.println("\""+parent.getPath()+"\","+(float)size/1048576.0f+","+mtime/1000L+","+fcount);
		return;
	}
}
