package jalgs.jseq;

import jalgs.algutils;
import jalgs.gui_interface;
import jalgs.jdataio;
import jalgs.jsort;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class SlideSeqDecoding{
	
	public String bridgecode="TCTTCAGCGTTCCCGAGA";
	public gui_interface gi;
	
	public SlideSeqDecoding(gui_interface gi){
		if(gi!=null) this.gi=gi;
	}

	//here we have static methods to decode slide-seq data
	public Object[] getBarCodes(String fastqpath){
		//String[] flines=null;
		String[][] barcodes=null;
		int[] indices=null; //this keeps track of the sequence indices for paired end matching of other data sets
		int nseqs=0;
		int nvalid=0;
		System.out.println("Getting Number of Reads");
		try{
    		BufferedReader b=new BufferedReader(new FileReader(fastqpath));
    		int nlines=0;
    		while(b.readLine()!=null){
    			nlines++;
    		}
    		b.close();
    		System.out.println("Reading "+nlines+" Lines");
    		b=new BufferedReader(new FileReader(fastqpath));
    		//flines=(new jdataio()).readstringfilelines(b);
    		//flines=new String[nlines/2]; //assume that 1/2 of the lines are headers or spacers
    		barcodes=new String[nlines/4][2]; //1/4 of the sequences don't have bridges
    		indices=new int[nlines/4];
    		for(int i=0;i<nlines/4;i++){
    			b.readLine();
    			String temp=b.readLine();
    			b.readLine();
    			b.readLine();
    			//if(temp.startsWith("@")) continue;
    			//if(temp.startsWith("#")) continue;
    			//if(temp.startsWith("+")) continue;
    			nseqs++;
    			//int bridgeindex=temp.indexOf(bridgecode);
    			//try finding the the bridge with mismatches
    			int bridgeindex=indexOfMismatch(temp,bridgecode,8,12,1);
    			int bridgeend=bridgeindex+bridgecode.length()+6;
    			if(bridgeindex>7 && bridgeend<=(temp.length()-1)){
    				String barcode1=temp.substring((bridgeindex-8),bridgeindex);
    				String barcode2=temp.substring((bridgeindex+bridgecode.length()),bridgeend);
    				String seq2=temp.substring(bridgeend,temp.length());
    				String barcode=barcode1+barcode2;
    				barcodes[nvalid][0]=barcode;
    				barcodes[nvalid][1]=seq2;
    				indices[nvalid]=i;
    				nvalid++;
    				//List<String> entry=new ArrayList<String>();
    				//entry.add(barcode); entry.add(seq2);
    				//barcodes.add(entry);
    			}
    			if(i%1000000==0) showProgress(i,nlines/4);
    		}
    		b.close();
		} catch(Exception e){
			System.out.println("error reading file on seq "+nseqs);
			System.out.println(e.getMessage());
		}
		//now get all of the sequences
		//they don't start with @ or # or +
		//boolean[] valid=new boolean[flines.length];
		//List<String> seqs=new ArrayList<String>();
		/*System.out.println("Finding Valid Reads");
		int nseqs=0;
		for(int i=0;i<flines.length;i++){
			if(flines[i].startsWith("@")) flines[i]=null;
			if(flines[i].startsWith("#")) flines[i]=null;
			if(flines[i].startsWith("+")) flines[i]=null;
			if(flines[i]!=null) nseqs++;
			//seqs.add(flines[i]);
			//valid[i]=(flines[i]!=null);
		}*/
		//now search for the bridge string code in each barcode
		/*System.gc();
		System.out.println("Found "+nseqs+" Sequences");
		System.out.println("Searching For Barcodes");
		//List<List<String>> barcodes=new ArrayList<List<String>>();
		String[][] barcodes=new String[nseqs][2];
		int nvalid=0;
		//String bridgecode="TCTTCAGCGTTCCCGAGA";
		for(int i=0;i<flines.length;i++){
			if(flines[i]==null) continue;
			String temp=flines[i];
			int bridgeindex=temp.indexOf(bridgecode);
			if(bridgeindex>7){
				String barcode1=temp.substring((bridgeindex-6),bridgeindex);
				String barcode2=temp.substring((bridgeindex+bridgecode.length()),(bridgeindex+bridgecode.length()+7));
				String seq2=temp.substring((bridgeindex+bridgecode.length()+7),temp.length());
				String barcode=barcode1+barcode2;
				barcodes[nvalid][0]=barcode;
				barcodes[nvalid][1]=seq2;
				nvalid++;
				//List<String> entry=new ArrayList<String>();
				//entry.add(barcode); entry.add(seq2);
				//barcodes.add(entry);
			}
			if(i%100000==0) gi.showProgress(i,flines.length);
		}*/
		System.out.println("finished decoding "+nvalid+" barcodes");
		String[][] validbarcodes=new String[nvalid][2];
		int[] validindices=new int[nvalid];
		for(int i=0;i<nvalid;i++){
			validbarcodes[i]=barcodes[i];
			validindices[i]=indices[i];
		}
		return new Object[]{validbarcodes,validindices};
	}
	
	public void showProgress(int pos,int finalpos){
		if(gi!=null){
			gi.showProgress(pos,finalpos);
		} else {
			System.out.print("position "+pos+" of "+finalpos+"\r");
		}
	}
	
	public void showMessage(String message){
		if(gi!=null){
			gi.showMessage(message);
		} else {
			System.out.println(message);
		}
	}
	
	public Object[] matchPositions(String[][] barcodes1,int[] indices,List<String> beadcodes,float[][] positions,boolean dedup){
		//this matches the bead codes that go with each barcode and their positions
		//start by collecting the identical barcodes together
		//first sort the table
		System.out.println("Sorting the barcode table");
		List<Object[]> barcodes=get_sorted_listtable(barcodes1,indices,0);
		//at this point the barcodes have the original string list (sorted on the first string) plus the integer position
		//now get the unique list of barcodes
		System.out.println("Getting the unique barcode list");
		Object[] uniqueoutput=get_cell_list(barcodes,0);
		List<String> uniquebarcodes=(List<String>)uniqueoutput[0];
		int[] uniquecounts=(int[])uniqueoutput[1];
		int[] uniquestarts=(int[])uniqueoutput[2];
		System.out.println("Found "+uniquebarcodes.size()+" unique barcodes");
		//now sort the beadcodes and search for our unique barcodes
		int[] beadorder=jsort.get_javasort_order(beadcodes);
		List<String> sortbeadcodes=new ArrayList<String>();
		float[][] sortpositions=new float[positions.length][];
		for(int i=0;i<positions.length;i++){
			sortbeadcodes.add(beadcodes.get(beadorder[i]));
			sortpositions[i]=positions[beadorder[i]];
		}
		//need to make a big index with all positions singly substituted so that we can map single mismatches
		List<Object[]> subbedbarcodes=makeSortedSubbedList(sortbeadcodes);
		System.out.println("Matching the barcodes");
		int[] matchindex=new int[uniquebarcodes.size()];
		int[] cellcounts=new int[beadcodes.size()];
		//int[] cellstarts=new int[beadcodes.size()];
		int[][] cellindices=new int[beadcodes.size()][];
		int matchcount=0;
		for(int i=0;i<uniquebarcodes.size();i++){
			String temp=uniquebarcodes.get(i);
			matchindex[i]=find_sorted_listtable_string(sortbeadcodes,temp);
			if(matchindex[i]>=0){
				cellcounts[matchindex[i]]+=uniquecounts[i];
				//cellstarts[matchindex[i]]=uniquestarts[i];
				cellindices[matchindex[i]]=addIndicesToArray(cellindices[matchindex[i]],uniquestarts[i],uniquecounts[i]);
				matchcount++;
			} else {
				//search for singly mismatched barcodes
				//here's the problem, this may end up merging some of our unique codes which will be out of order in our original sorted list
				//to do this, we may need to rehash the way we do accounting
				//perhaps need to keep explicit indices for every barcode?
				/*matchindex[i]=find_sorted_listtable_string_mismatch(subbedbarcodes,temp);
				if(matchindex[i]>=0){
					cellcounts[matchindex[i]]+=uniquecounts[i];
					//cellstarts[matchindex[i]]=uniquestarts[i];
					cellindices[matchindex[i]]=addIndicesToArray(cellindices[matchindex[i]],uniquestarts[i],uniquecounts[i]);
					matchcount++;
				}*/
			}
			//if(i%100000==0) showProgress(i,uniquebarcodes.size());
		}
		System.out.println("Matched "+matchcount+" barcodes");
		System.out.println("Assembling Matched Statistics");
		//now output a table of barcodes and counts and positions
		ArrayList<List<String>> outpositions=new ArrayList<List<String>>();
		//at the same time, filter the barcode/umi list to only contain the matches
		List<Object[]> matchedlist=new ArrayList<Object[]>();
		int pos=-1;
		for(int i=0;i<uniquebarcodes.size();i++){
			if(matchindex[i]>=0){
				//output is barcode, beadcode, xpos, ypos, count
				int codepos=matchindex[i];
				List<String> row=new ArrayList<String>();
				String tbarcode=uniquebarcodes.get(i);
				row.add(tbarcode);
				row.add(sortbeadcodes.get(codepos));
				String xpos=""+sortpositions[codepos][0];
				String ypos=""+sortpositions[codepos][1];
				row.add(xpos);
				row.add(ypos);
				//need to wait and add counts when UMI's are deduplicated
				//int codecount=cellcounts[codepos];
				//row.add(""+codecount);
				//outpositions.add(row);
				//need to deduplicate UMI's within the barcode list and add to the matched list
				List<Object[]> barcodelist=new ArrayList<Object[]>();
				for(int j=0;j<cellcounts[codepos];j++){
					//barcodelist.add(barcodes.get(j+cellstarts[codepos]));
					if(cellindices[codepos]==null) System.out.println("null on "+i);
					else barcodelist.add(barcodes.get(cellindices[codepos][j]));
				}
				List<Object[]> tlist=null;
				if(dedup){
					tlist=deDuplicateUMIs(barcodelist,sortpositions[codepos]);
				} else {
					tlist=addPosition(barcodelist,sortpositions[codepos]);
				}
				matchedlist.addAll(tlist);
				row.add(""+tlist.size());
				outpositions.add(row);
			}
		}
		//on output the matched list has the original string list, the integer index, and the float[] position
		return new Object[]{outpositions,uniquebarcodes,matchedlist};
	}
	
	public static int[] addIndicesToArray(int[] arr,int start,int count){
		int[] temp=new int[count];
		for(int i=0;i<count;i++) temp[i]=start+i;
		if(arr==null) arr=temp;
		else arr=(int[])algutils.combine_arrays(arr,temp);
		return arr;
	}
	
	public static List<Object[]> deDuplicateUMIs(List<Object[]> barcodelist,float[] pos){
		//this deduplicates the umi's and adds the position to the object array
		//first sort by the UMI's
		sort_listtable2(barcodelist,1);
		List<Object[]> deduplist=new ArrayList<Object[]>();
		//now find unique ones and delete the duplicates
		Object[] currcello=barcodelist.get(0);
		String currumi=((String[])currcello[0])[1];
		for(int i=1;i<barcodelist.size();i++){
			Object[] tempo=barcodelist.get(i);
			String temp=((String[])tempo[0])[1];
			if(!temp.equals(currumi)){
				Object[] tempo2={tempo[0],tempo[1],pos};
				deduplist.add(tempo2);
				currumi=temp;
			}
		}
		return deduplist;
	}
	
	public static List<Object[]> addPosition(List<Object[]> barcodelist,float[] pos){
		//this just adds the position without deduplication
		//sort_listtable2(barcodelist,1);
		List<Object[]> tlist=new ArrayList<Object[]>();
		//now find unique ones and delete the duplicates
		for(int i=0;i<barcodelist.size();i++){
			Object[] tempo=barcodelist.get(i);
			Object[] tempo2={tempo[0],tempo[1],pos};
			tlist.add(tempo2);
		}
		return tlist;
	}
	
	public static void sort_listtable(List<List<String>> list,int sortcolumn){
		final int tempcolumn=sortcolumn;
		Collections.sort(list,new Comparator<List<String>>(){
			public int compare(List<String> o1,List<String> o2){
				return o1.get(tempcolumn).compareTo(o2.get(tempcolumn));
				//String temp1=o1.get(tempcolumn).trim().toUpperCase();
				//String temp2=o2.get(tempcolumn).trim().toUpperCase();
				//return temp1.compareTo(temp2);
			}
		});
	}
	
	public static void sort_listtable2(List<Object[]> list,int sortcolumn){
		//in this version each list row contains a string array and an integer
		final int tempcolumn=sortcolumn;
		Collections.sort(list,new Comparator<Object[]>(){
			public int compare(Object[] o1,Object[] o2){
				String s1=((String[])o1[0])[tempcolumn];
				String s2=((String[])o2[0])[tempcolumn];
				return s1.compareTo(s2);
			}
		});
	}
	
	public static void sort_listtable2(List<Object[]> list){
		//in this version each list row contains a string and an integer
		Collections.sort(list,new Comparator<Object[]>(){
			public int compare(Object[] o1,Object[] o2){
				String s1=((String)o1[0]);
				String s2=((String)o2[0]);
				return s1.compareTo(s2);
			}
		});
	}
	
	public static List<Object[]> get_sorted_listtable(String[][] list,int[] indices,int sortcolumn){
		List<Object[]> listtable=new ArrayList<Object[]>();
		for(int i=0;i<list.length;i++){
			Object[] row={list[i],Integer.valueOf(indices[i])};
			listtable.add(row);
		}
		sort_listtable2(listtable,sortcolumn);
		return listtable;
	}
	
	public static int find_sorted_listtable_string(List<String> list,String target){
		// note that the listtable must be sorted for this to work
		//searches for exact match but right justified
		int index=Collections.binarySearch(list,target,new Comparator<String>(){
			public int compare(String o1,String o2){
				int lendiff=o1.length()-o2.length();
				if(lendiff==0) return o1.compareTo(o2);
				//if(lendiff==0) return compareToMismatch(o1,o2,0);
				if(lendiff>0){
					String temp=o1.substring(lendiff,o1.length());
					return temp.compareTo(o2);
					//return compareToMismatch(temp,o2,0);
				} else {
					String temp=o2.substring(-lendiff,o2.length());
					return temp.compareTo(o1);
					//return compareToMismatch(temp,o1,0);
				}
			}
		});
		return index;
	}
	
	public static int find_sorted_listtable_string_mismatch(List<Object[]> subbedlist,String target){
		//brute force searching for string with same length and single mismatch
		//takes waaaaay too long with brute force
		//how about we sub each nucleotide in turn for the target nucleotide (already there in subbedlist)
		//assume we've already failed at binary search
		//to get the original index, the subbed list carries an index with it
		Object[] targetobj=new Object[]{target,Integer.valueOf(0)};
		int index=Collections.binarySearch(subbedlist,targetobj,new Comparator<Object[]>(){
			public int compare(Object[] o1,Object[] o2){
				String s1=(String)o1[0];
				String s2=(String)o2[0];
				return s1.compareTo(s2);
			}
		});
		if(index>=0){
			Object[] temp=subbedlist.get(index);
			return ((Integer)temp[1]).intValue();
		} else {
			return -1;
		}
	}
	
	public static List<Object[]> makeSortedSubbedList(List<String> list){
		//here we make a big long list from an existing list with every character subbed out by every alternative character
		//characters are in rotating list ATCG
		//add in the original indexes to keep things in order
		List<Object[]> subbed=new ArrayList<Object[]>();
		int nchars=list.get(0).length();
		for(int i=0;i<list.size();i++){
			char[] temp=list.get(i).toCharArray();
    		for(int j=0;j<nchars;j++){
        		for(int k=0;k<3;k++){
        			char[] subbed2=temp.clone();
        			subbed2[j]=getShiftedChar(subbed2[j],k+1);
        			subbed.add(new Object[]{new String(subbed2),Integer.valueOf(i)});
        		}
    		}
		}
		sort_listtable2(subbed);
		return subbed;
	}
	
	public static char getShiftedChar(char currchar,int nshift){
		if(nshift==1){
			if(currchar=='A') return 'T';
			if(currchar=='T') return 'C';
			if(currchar=='C') return 'G';
			if(currchar=='G') return 'A';
		} else if(nshift==2){
			if(currchar=='A') return 'C';
			if(currchar=='T') return 'G';
			if(currchar=='C') return 'A';
			if(currchar=='G') return 'T';
		} else {
			if(currchar=='A') return 'G';
			if(currchar=='T') return 'A';
			if(currchar=='C') return 'T';
			if(currchar=='G') return 'C';
		}
		return 'A';
	}
	
	public static int indexOfMismatch(String bait,String query,int start,int end,int maxmismatch){
		//this is a brute force finder with n mismatches
		//make sure start and end are close together for reasonable performance
		//returns -1 if query isn't found
		//need to make sure that start+querylength and end+querylength are contained
		int minmismatch=query.length();
		int minpos=start;
		for(int i=start;i<end;i++){
			int nmismatch=0;
			for(int j=0;j<query.length();j++){
				if(bait.charAt(i+j)!=query.charAt(j)) nmismatch++;
			}
			if(nmismatch<minmismatch){
				minmismatch=nmismatch;
				minpos=i;
			}
		}
		if(minmismatch<=maxmismatch) return minpos;
		else return -1;
	}
	
	public static int compareToMismatch(String bait,String query,int maxmismatch){
		//this works but doesn't support binary search
		int nmismatch=countMismatches(bait,query);
		if(nmismatch<=maxmismatch) return 0;
		else return bait.compareTo(query);
	}
	
	public static int countMismatches(String bait,String query){
		//counts mismatches in strings of same length
		int nmismatch=0;
		for(int j=0;j<query.length();j++){
			if(bait.charAt(j)!=query.charAt(j)) nmismatch++;
		}
		return nmismatch;
	}
	
	public static Object[] get_cell_list(List<Object[]> list,int cellcolumn){
		List<String> celllist=new ArrayList<String>();
		int[] counts=new int[list.size()];
		int[] starts=new int[list.size()];
		int pos=0;
		Object[] currcello=list.get(0);
		String currcell=((String[])currcello[0])[cellcolumn];
		celllist.add(currcell);
		counts[pos]++;
		for(int i=1;i<list.size();i++){
			Object[] tempo=list.get(i);
			String temp=((String[])tempo[0])[cellcolumn];
			if(!temp.equals(currcell)){
				celllist.add(temp);
				currcell=temp;
				pos++;
				starts[pos]=i;
			}
			counts[pos]++;
		}
		return new Object[]{celllist,algutils.get_subarray(counts,0,pos+1),algutils.get_subarray(starts,0,pos+1)};
	}
	
	public String[][] getSecondarySeqs(List<Object[]> barcodes,String path){
		//the barcode objects have the original string array, the integer index, and the float[] position
		int[] indices=new int[barcodes.size()];
		for(int i=0;i<indices.length;i++){
			Object[] tempo=barcodes.get(i);
			indices[i]=(Integer)tempo[1];
		}
		int[] sortorder=jsort.get_javasort_order(indices);
		System.out.println("min index = "+indices[sortorder[0]]);
		System.out.println("max index = "+indices[sortorder[indices.length-1]]);
		/*int checktemp=indices[sortorder[0]];
		for(int i=1;i<indices.length;i++){
			int checktemp2=indices[sortorder[i]];
			if(checktemp>checktemp2){
				System.out.println("check "+checktemp+" , "+checktemp2);
			}
			checktemp=checktemp2;
		}*/
		String[][] seqs=new String[3][indices.length];
		//here we loop through the secondary file and find the sequences that correspond to our indices
		int nseqs=-1;
		try{
    		BufferedReader b=new BufferedReader(new FileReader(path));
    		int nlines=0;
    		while(b.readLine()!=null){
    			nlines++;
    		}
    		b.close();
    		b=new BufferedReader(new FileReader(path));
    		int pos=0;
    		int currindex=indices[sortorder[pos]];
    		for(int i=0;i<nlines/4;i++){
        		String header=b.readLine();
        		String temp=b.readLine();
        		b.readLine();
        		String qual=b.readLine();
        		nseqs++;
        		if(i>=currindex){
        			//lets add a bunch of stuff to the header separated by underscores
        			//add _barcode_umi_xpos_ypos
        			Object[] bdetails=barcodes.get(sortorder[pos]);
        			String[] btemp=(String[])bdetails[0];
        			float[] bpositions=(float[])bdetails[2];
        			header=header.replace(' ','_');
        			String theader=header+"_"+btemp[0]+"_"+btemp[1]+"_"+bpositions[0]+"_"+bpositions[1];
        			seqs[0][sortorder[pos]]=temp;
        			seqs[1][sortorder[pos]]=theader;
        			seqs[2][sortorder[pos]]=qual;
        			pos++;
        			if(pos>=sortorder.length) break;
        			/*if(currindex>indices[sortorder[pos]]){
        				System.out.println("Sort failed at "+currindex+" and "+indices[sortorder[pos]]);
        			}*/
        			currindex=indices[sortorder[pos]];
        			//System.out.println(currindex);
        		}
    		}
    		b.close();
		} catch(Exception e){
			System.out.println("error reading file on seq "+nseqs);
			System.out.println(e.getMessage());
		}
		return seqs;
	}
	
	public static void main(String[] args){
		//args should be barcode fastq file, secondary (3') fastq file, and "dedupumis" if deduplication is desired
		String inname=args[0];
		SlideSeqDecoding decoder=new SlideSeqDecoding(null);
		Object[] barcodeoutput=decoder.getBarCodes(inname);
		String[][] barcodes=(String[][])barcodeoutput[0];
		int[] indices=(int[])barcodeoutput[1];
		String dir=(new File(inname)).getParent();
		String prefix=(new File(inname)).getName();
		prefix=prefix.substring(0,prefix.length()-6);
		List<String> beadcodes=null;
		try{
			BufferedReader b=new BufferedReader(new FileReader(dir+File.separator+"BeadBarcodes_nocommas.txt"));
			beadcodes=(new jdataio()).readstringfilelines2(b);
			b.close();
			for(int i=0;i<beadcodes.size();i++){
				beadcodes.set(i,beadcodes.get(i).trim());
			}
		} catch(Exception e){
			System.out.println(e.getMessage());
			return;
		}
		System.out.println("read "+beadcodes.size()+" beadcodes");
		float[][] positions=null;
		try{
			BufferedReader b=new BufferedReader(new FileReader(dir+File.separator+"BeadLocations.txt"));
			String temp=(new jdataio()).readstringfile(b);
			b.close();
			String[] xystrings=temp.split("\n");
			String[] xposstrings=xystrings[0].split(",");
			String[] yposstrings=xystrings[1].split(",");
			System.out.println("read "+xposstrings.length+" positions");
			int npos=xposstrings.length;
			if(npos!=beadcodes.size()){
				System.out.println("wrong number of positions");
				npos=Math.min(npos,beadcodes.size());
			}
			positions=new float[npos][2];
			for(int i=0;i<npos;i++){
				positions[i][0]=Float.parseFloat(xposstrings[i].trim());
				positions[i][1]=Float.parseFloat(yposstrings[i].trim());
			}
		} catch(Exception e){
			System.out.println(e.getMessage());
			return;
		}
		System.gc();
		boolean dedupumis=(args.length>2 && args[2].equals("dedupumis"));
		Object[] matchout=decoder.matchPositions(barcodes,indices,beadcodes,positions,dedupumis);
		List<List<String>> matches=(List<List<String>>)matchout[0];
		List<String> uniquebarcodes=(List<String>)matchout[1];
		List<Object[]> matchedbarcodes=(List<Object[]>)matchout[2];
		String outname=prefix+"_decoded.csv";
		System.out.println("writing output to "+outname);
		try{
			BufferedWriter b=new BufferedWriter(new FileWriter(dir+File.separator+outname));
			//output is barcode, beadcode, xpos, ypos, count
			b.write("Barcode,Beadcode,xpos,ypos,count\n");
			(new jdataio()).writestringfile(b,matches);
			b.close();
		} catch(Exception e){
			System.out.println(e.getMessage());
		}
		//now want to pull out the sequences from the 3' file that match our barcodes and UMI's
		System.out.println("finding 3' sequences in "+args[1]);
		String[][] secondaryseqs=decoder.getSecondarySeqs(matchedbarcodes,args[1]);
		String outname3=prefix+"_matchedbarcodes.csv";
		System.out.println("writing output to "+outname3);
		try{
			BufferedWriter b=new BufferedWriter(new FileWriter(dir+File.separator+outname3));
			b.write("barcode,UMI,readIndex,secondarySeq,xpos,ypos\n");
			for(int i=0;i<matchedbarcodes.size();i++){
				Object[] templineo=matchedbarcodes.get(i);
				String[] templine2=(String[])templineo[0];
				String templine=templine2[0];
				for(int j=1;j<templine2.length;j++){
					templine+=","+templine2[j];
				}
				templine+=","+((Integer)templineo[1]).toString();
				templine+=","+secondaryseqs[0][i];
				float[] temppos=(float[])templineo[2];
				templine+=","+temppos[0];
				templine+=","+temppos[1];
				b.write(templine+"\n");
			}
			//(new jdataio()).writestringfile3(b,matchedbarcodes);
			b.close();
		} catch(Exception e){
			System.out.println(e.getMessage());
		}
		String prefix2=(new File(args[1])).getName();
		prefix2=prefix2.substring(0,prefix2.length()-6);
		System.out.println("outputting 3' fastq file");
		String outname4=prefix2+"_matchedbarcodes.fastq";
		try{
			BufferedWriter b=new BufferedWriter(new FileWriter(dir+File.separator+outname4));
			//b.write("barcode,UMI,readIndex,secondarySeq\n");
			for(int i=0;i<matchedbarcodes.size();i++){
				b.write(secondaryseqs[1][i]+"\n");
				b.write(secondaryseqs[0][i]+"\n");
				b.write("+\n");
				b.write(secondaryseqs[2][i]+"\n");
			}
			//(new jdataio()).writestringfile3(b,matchedbarcodes);
			b.close();
		} catch(Exception e){
			System.out.println(e.getMessage());
		}
		System.out.println("finished");
	}

}

class MatchBarCodes implements Runnable{
	private String query;
	private List<String> beadcodes;
	int[] matchindices;
	int index;
	
	public MatchBarCodes(String query,List<String> beadcodes,int[] matchindices,int index){
		this.query=query;
		this.beadcodes=beadcodes;
		this.index=index;
		this.matchindices=matchindices;
	}
	
	public void run(){
		// note that the listtable must be sorted for this to work
		//searches for exact match but right justified
		int temp=Collections.binarySearch(beadcodes,query,new Comparator<String>(){
			public int compare(String o1,String o2){
				int lendiff=o1.length()-o2.length();
				if(lendiff==0) return o1.compareTo(o2);
				if(lendiff>0){
					String temp=o1.substring(lendiff,o1.length());
					return temp.compareTo(o2);
				} else {
					String temp=o2.substring(-lendiff,o2.length());
					return temp.compareTo(o1);
				}
			}
		});
		matchindices[index]=temp;
	}
	
}
