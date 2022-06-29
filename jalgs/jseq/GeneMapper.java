package jalgs.jseq;

import jalgs.table_tools2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GeneMapper{
	
	/******************
	 * here we map gene alignments from a sam file to gene interval
	 * @param sampath: path to sam file that we are mapping
	 * @param gipath: path to gene interval file
	 */
	public List<List<String>> mapAlignedGeneLocations(String sampath,String gipath){
		System.out.println("reading sam file");
		String[][] samdata=readSAMFile(new File(sampath));
		String[] sam_labels={"header","flag","chr","pos","qual","cigar","rnext","pnext","tlen","seq","oqual","xa","md","nm"};
		//lets remove the unmapped reads (flag==4) and those with atypical chromosomes
		//at the same time, we can delete the unnecessary columns (for now everything past pos)
		System.out.println("filtering sam file");
		List<List<String>> samfiltered=new ArrayList<List<String>>();
		for(int i=0;i<samdata.length;i++){
			String[] row=samdata[i];
			int flag=Integer.parseInt(row[1]);
			if(flag!=4 && row[2].length()<6){
				String[] sublist={row[0],row[1],row[2],row[3]};
				samfiltered.add(table_tools2.stringarray2list(sublist));
			}
			row=null;
		}
		samdata=null;
		System.gc();
		System.out.println("sorting sam file");
		table_tools2.sort_listtable(samfiltered,3,true); //sort on position
		table_tools2.sort_listtable(samfiltered,1,true); //sort on flag (strand)
		table_tools2.sort_listtable(samfiltered,2); //sort on chromosome
		//delete qual, cigar, rnext, pnext, tlen, xa, md, and nm columns
		/*table_tools2.delete_listtable_column(samfiltered,4);
		table_tools2.delete_listtable_column(samfiltered,4);
		table_tools2.delete_listtable_column(samfiltered,4);
		table_tools2.delete_listtable_column(samfiltered,4);
		table_tools2.delete_listtable_column(samfiltered,4);
		table_tools2.delete_listtable_column(samfiltered,6);
		table_tools2.delete_listtable_column(samfiltered,6);
		table_tools2.delete_listtable_column(samfiltered,6);*/
		System.out.println("reading and sorting gi file");
		String[] samfiltered_labels={"header","flag","chr","pos","matchedname"};
		List<List<String>> gidf1=table_tools2.getTableFromFile(new File(gipath),",",false);
		String[] gilabels=gidf1.get(0).toArray(new String[0]);
		//these should be chr, start, end, strand (+ or -), gene_id
		gidf1.remove(0);
		List<List<String>> gidf=new ArrayList<List<String>>();
		//going to remove the atypical chromosomes (don't have many genes anyway)
		for(int i=0;i<gidf1.size();i++){
			if(gidf1.get(i).get(0).length()<6) gidf.add(gidf1.get(i));
		}
		gidf1=null;
		//now sort and filter the gene interval table
		table_tools2.sort_listtable(gidf,1,true); //sort on start
		table_tools2.sort_listtable(gidf,3); //sort on strand
		table_tools2.sort_listtable(gidf,0); //sort on chromosome
		//String gidir=(new File(gipath)).getParent();
		//table_tools2.writeTableToFile(gidir+File.separator+"sorted.csv",gilabels,gidf,1);
		Object[] giindices=getGIIndices(gidf);
		System.out.println("finding intervals");
		int prevpos=0;
		String prevchr="chr0";
		int chrindex=0;
		String prevstrand="+";
		String prevname="";
		int nfiltered=samfiltered.size();
		for(int i=0;i<samfiltered.size();i++){
			String chr=samfiltered.get(i).get(2);
			if(!chr.equals(prevchr)){
				chrindex=((List<String>)giindices[0]).indexOf(chr);
			}
			String strand=(samfiltered.get(i).get(1).equals("0")?"+":"-");
			int genepos=table_tools2.get_integer(samfiltered,i,3);
			String matchedname=prevname;
			if(prevpos!=genepos || !strand.equals(prevstrand) || !chr.equals(prevchr)){
				matchedname=geneFinder(chrindex,strand,genepos,gidf,giindices);
			}
			samfiltered.get(i).add(matchedname);
			prevname=matchedname; prevpos=genepos; prevchr=chr; prevstrand=strand;
			if((i%100000)==0) System.out.println("line "+i+" of "+nfiltered+"\r");
		}
		System.out.println("\n");
		samfiltered.add(0,table_tools2.stringarray2list(samfiltered_labels));
		return samfiltered;
	}
	
	public Object[] getGIIndices(List<List<String>> gidf){
		//now create the gene interval indices (start indices for each unique chromosome and its reverse strand)
		//needs to be sorted by start, strand, and chromsome
		List<String> uniquechroms=table_tools2.get_cell_list(gidf,0);
		int[] chrindex=new int[uniquechroms.size()];
		int[] minusstrandindex=new int[chrindex.length];
		String currchrom=uniquechroms.get(0);
		String currstrand=gidf.get(0).get(3);
		boolean strandmode=true;
		int pos=0;
		chrindex[pos]=0;
		for(int i=1;i<gidf.size();i++){
			String tchrom=gidf.get(i).get(0);
			String tstrand=gidf.get(i).get(3);
			if(strandmode){
    			if(!tstrand.equals(currstrand)){
    				minusstrandindex[pos]=i;
    				currstrand=tstrand;
    				strandmode=false;
    			}
			} else {
				if(!tchrom.equals(currchrom)){
					pos++;
					chrindex[pos]=i;
					currchrom=tchrom;
					currstrand=tstrand;
					strandmode=true;
				}
			}
		}
		//System.out.println(table_tools.print_string_array(uniquechroms));
		//System.out.println(table_tools.print_int_array(chrindex));
		//System.out.println(table_tools.print_int_array(minusstrandindex));
		return new Object[]{uniquechroms,chrindex,minusstrandindex};
	}
	
	public String geneFinder(String chr,String strand,int genepos,List<List<String>> gidf,Object[] giindices){
		int chromindex=((List<String>)giindices[0]).indexOf(chr);
		return geneFinder(chromindex,strand,genepos,gidf,giindices);
	}
	
	public String geneFinder(int chromindex,String strand,int genepos,List<List<String>> gidf,Object[] giindices){
		//here we use the gene interval indices to help find which gene interval our gene position belongs to
		int nchroms=((List<String>)giindices[0]).size();
		//int chromindex=((List<String>)giindices[0]).indexOf(chr);
		int startindex=((int[])giindices[1])[chromindex];
		int endindex=((int[])giindices[2])[chromindex];
		if(strand.equals("-")){
			startindex=endindex;
			if(chromindex<(nchroms-1)){
				endindex=((int[])giindices[1])[chromindex+1];
			} else {
				endindex=gidf.size();
			}
		}
		for(int i=startindex;i<endindex;i++){
			int intstart=table_tools2.get_integer(gidf,i,1);
			int intend=table_tools2.get_integer(gidf,i,2);
			if(genepos<=intend && genepos>=intstart){
				return gidf.get(i).get(4);
			}
			if(intstart>genepos) break;
		}
		return "genenotfound";
	}
	
	public String[][] readSAMFile(File infile){
		try{
			BufferedReader b=new BufferedReader(new FileReader(infile));
			int nlines=0;
			String temp=b.readLine();
			while(b.readLine()!=null){
				nlines++;
			}
			b.close();
			b=new BufferedReader(new FileReader(infile));
			//List<String> lines=new ArrayList<String>();
			//start by reading through the @lines that start the file
			int headlines=0;
			temp=b.readLine();
			while(temp.startsWith("@")){
				headlines++;
				temp=b.readLine();
			}
			String[][] lines=new String[nlines-headlines+1][];
			int pos=0;
			//now read the rest of the lines
			while(temp!=null && temp.length()>1){
				String[] split=table_tools2.split(temp,"\t");
				lines[pos]=split;
				pos++;
				temp=b.readLine();
			}
			//List<List<String>> listtable=table_tools2.table2listtable(lines,"\t",false);
			//lines=null;
			b.close();
			return lines;
		}
		catch(IOException e){
			System.out.println(e.getMessage());
			return null;
		}
	}
	
	public List<List<List<String>>> filterAnnotations(List<List<String>> annotated){
		//once we have our annotated list, we can do some post processing
		//want to filter out genenotfound
		//want to combine the UMI's (eliminate the ones with genenotfound)
		//finally build an expression list for each bead
		//start by eliminating genenotfound and parsing the header for barcode,umi,xpos,ypos
		List<List<String>> filtered=new ArrayList<List<String>>();
		int ncols=annotated.get(0).size();
		for(int i=1;i<annotated.size();i++){
			List<String> row=annotated.get(i);
			if(!row.get(4).equals("genenotfound")){
				String[] headinfo=table_tools2.split(row.get(0),"_",false);
				int temp1=headinfo.length;
				List<String> newrow=row.subList(1,ncols);
				newrow.add(headinfo[temp1-4]);
				newrow.add(headinfo[temp1-3]);
				newrow.add(headinfo[temp1-2]);
				newrow.add(headinfo[temp1-1]);
				filtered.add(newrow);
			}
		}
		//now the columns are 0flag, 1chr, 2pos, 3matchedname,4barcode,5umi,xpos,ypos
		//now deduplicate the UMI's, sort first
		table_tools2.sort_listtable(filtered,5);
		table_tools2.sort_listtable(filtered,4);
		for(int i=0;i<10;i++){
			System.out.println(table_tools2.print_string_array(filtered.get(i),1));
		}
		List<List<List<String>>> deduped=new ArrayList<List<List<String>>>(); //order is bead, unique transcript, details with transcript count
		List<List<String>> templist=new ArrayList<List<String>>(); //this will hold the bead codes as we accumulate them
		String currumi=filtered.get(0).get(5);
		String currcode=filtered.get(0).get(4);
		List<List<String>> umilist=new ArrayList<List<String>>(); //this will hold each umi as we accumulate it
		umilist.add(filtered.get(0));
		for(int i=1;i<filtered.size();i++){
			String tumi=filtered.get(i).get(5);
			String tcode=filtered.get(i).get(4);
			if(tcode.equals(currcode)){
				if(tumi.equals(currumi)){
					umilist.add(filtered.get(i)); //add to the current umi list
				} else {
					//we reached the end of this umi list, project it
					List<String> bestumi=projectGene(umilist);
					if(bestumi!=null) templist.add(bestumi);
					else System.out.println("empty umi at line "+i);
					currumi=tumi;
					umilist=new ArrayList<List<String>>();
					umilist.add(filtered.get(i));
				}
			} else {
				//we reached the end of this umi list, project it
				List<String> bestumi=projectGene(umilist); //shouldn't be null but sometimes is
				if(bestumi!=null) templist.add(bestumi);
				else System.out.println("empty umi at line "+i);
				//sort the templist by gene id and count 
				List<List<String>> uniquelist=accumulateBeadGenes(templist);
				deduped.add(uniquelist);
				currcode=tcode;
				templist=new ArrayList<List<String>>();
				currumi=tumi;
				umilist=new ArrayList<List<String>>();
				umilist.add(filtered.get(i));
			}
		}
		//at the end, we have one last umi to project
		List<String> bestumi=projectGene(umilist); //shouldn't be null but sometimes is
		if(bestumi!=null) templist.add(bestumi);
		//sort the templist by gene id and count 
		List<List<String>> uniquelist=accumulateBeadGenes(templist);
		deduped.add(uniquelist);
		return deduped;
	}
	
	public List<String> projectGene(List<List<String>> umilist){
		//finds the best (most frequent?) gene from a umilist
		//but this needs to be fast, so just take the most common among the first 3
		if(umilist.size()<1){
			return null;
		}
		if(umilist.size()<3) return umilist.get(0); //if size is 1 or 2 we don't have a way to decide, just output the first
		if(umilist.get(0).get(3).equals(umilist.get(1).get(3))){
			//here 1 and 2 match
			return umilist.get(0);
		} else {
			//see if 2 and 3 match
			if(umilist.get(1).get(3).equals(umilist.get(2).get(3))){
				return umilist.get(1);
			} else {
				//no matches, just return the first one
				return umilist.get(0);
			}
		}
	}
	
	public List<List<String>> accumulateBeadGenes(List<List<String>> beadreads){
		//this takes all bead reads and returns a unique list with counts
		//first sort by geneid
		table_tools2.sort_listtable(beadreads,3);
		List<List<String>> unique=new ArrayList<List<String>>();
		List<String> currread=beadreads.get(0);
		String currgene=currread.get(3);
		int currcounts=1;
		for(int i=1;i<beadreads.size();i++){
			String tgene=beadreads.get(i).get(3);
			if(tgene.equals(currgene)){
				currcounts++;
			} else {
				currread.add(""+currcounts);
				unique.add(currread);
				currread=beadreads.get(i);
				currgene=currread.get(3);
				currcounts=1;
			}
		}
		return unique;
	}
	
	public static List<List<String>> flattenTables(List<List<List<String>>> nested,String[] collabels){
		//here we flatten a list of tables into a single table, adding an index as the last column
		//put column labels in first row
		List<List<String>> flattened=new ArrayList<List<String>>();
		for(int i=0;i<nested.size();i++){
			List<List<String>> ttable=nested.get(i);
			for(int j=0;j<ttable.size();j++){
				List<String> row=ttable.get(j);
				row.add(""+i);
				flattened.add(row);
			}
		}
		List<String> labels=table_tools2.stringarray2list(collabels);
		labels.add("index");
		flattened.add(0,labels);
		return flattened;
	}

	public static void main(String[] args){
		System.out.println("mapping alignments");
		List<List<String>> annotated=(new GeneMapper()).mapAlignedGeneLocations(args[0],args[1]);
		String dir=(new File(args[0])).getParent();
		String name=(new File(args[0])).getName();
		String prefix=name.substring(0,name.length()-4);
		System.out.println("writing output to "+prefix+"_mapped.csv");
		table_tools2.writeTableToFile(dir+File.separator+prefix+"_mapped.csv",annotated,1);
		System.out.println("filtering and counting reads");
		List<List<List<String>>> beadGenes=(new GeneMapper()).filterAnnotations(annotated);
		//need to rearray this to write it out
		String[] collabels2={"flag","chr","pos","matchedname","barcode","umi","xpos","ypos","readcount"};
		List<List<String>> flatGenes=GeneMapper.flattenTables(beadGenes,collabels2);
		List<String> head=flatGenes.get(0);
		flatGenes.remove(0);
		table_tools2.sort_listtable(flatGenes,3);
		flatGenes.add(0,head);
		System.out.println("writing output to "+prefix+"_mapped_counted.csv");
		table_tools2.writeTableToFile(dir+File.separator+prefix+"_mapped_counted.csv",flatGenes,1);
		System.out.println("finished");
	}

}
