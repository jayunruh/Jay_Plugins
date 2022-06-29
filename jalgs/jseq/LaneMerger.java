package jalgs.jseq;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import jalgs.table_tools2;

public class LaneMerger{

	public static Object[] getInflatedFile(String lane1path){
		List<List<String>> lane1=table_tools2.getTableFromFile(new File(lane1path),",",false);
		String[] collabels1=table_tools2.list2stringarray(lane1.get(0));
		lane1.remove(0);
		return new Object[]{inflateTables(lane1),collabels1};
	}
	
	public List<List<List<String>>> mergeLanes(String lane1path,String lane2path){
		Object[] lane1obj=getInflatedFile(lane1path);
		List<List<List<String>>> lane1inf=(List<List<List<String>>>)lane1obj[0];
		Object[] lane2obj=getInflatedFile(lane2path);
		List<List<List<String>>> lane2inf=(List<List<List<String>>>)lane2obj[0];
		mergeLanes(lane1inf,lane2inf);
		return lane1inf;
	}
	
	public void mergeLanes(List<List<List<String>>> lane1inf,List<List<List<String>>> lane2inf){
		//merges lane2 (inflated) into lane1
		//umi's have already been accumulated, just sum the genes together
		//if we are saturated, some genes will be overcounted (same umi in different lanes)
		List<String> barcodelist1=new ArrayList<String>();
		//lets make a bead barcode list here
		for(int i=0;i<lane1inf.size();i++){
			barcodelist1.add(lane1inf.get(i).get(0).get(4));
		}
		List<String> barcodelist2=new ArrayList<String>();
		//lets make a bead barcode list here
		for(int i=0;i<lane2inf.size();i++){
			barcodelist2.add(lane2inf.get(i).get(0).get(4));
		}
		boolean[] matched=new boolean[barcodelist2.size()];
		for(int i=0;i<lane1inf.size();i++){
			int matchindex=barcodelist2.indexOf(barcodelist1.get(i));
			if(matchindex>=0){
				//we have a match, lets merge beads
				mergeBeads(lane1inf.get(i),lane2inf.get(matchindex));
				matched[matchindex]=true;
			}
		}
		//need to add all of the unmatched beads to our list
		for(int i=0;i<matched.length;i++){
			if(!matched[i]){
				lane1inf.add(lane2inf.get(i));
			}
		}
		return;
	}
	
	public void mergeBeads(List<List<String>> bead1,List<List<String>> bead2){
		//add bead2 reads into bead1, or add them to the list if they don't exist
		//sort our beads by genename to make matching easier
		table_tools2.sort_listtable(bead1,3);
		table_tools2.sort_listtable(bead2,3);
		List<String> geneList1=getColumn(bead1,3);
		List<String> geneList2=getColumn(bead2,3);
		boolean[] matched=new boolean[bead2.size()];
		for(int i=0;i<bead1.size();i++){
			int matchindex=geneList2.indexOf(geneList1.get(i));
			if(matchindex>=0){
				List<String> row=bead1.get(i);
				String added=addStringInts(row.get(8),bead2.get(matchindex).get(8));
				row.set(8,added);
				matched[matchindex]=true;
			}
		}
		for(int i=0;i<matched.length;i++){
			if(!matched[i]){
				bead1.add(bead2.get(i));
			}
		}
		return;
	}
	
	public static String addStringInts(String int1,String int2){
		return ""+(Integer.parseInt(int1)+Integer.parseInt(int2));
	}
	
	public static List<String> getColumn(List<List<String>> listtable,int selcol){
		List<String> col=new ArrayList<String>();
		for(int i=0;i<listtable.size();i++){
			col.add(listtable.get(i).get(selcol));
		}
		return col;
	}
	
	public static List<List<List<String>>> inflateTables(List<List<String>> flattened){
		//the last index of flattened dictates which subtable each row belongs to
		int ncols=flattened.get(0).size();
		table_tools2.sort_listtable(flattened,ncols-1,true);
		List<List<List<String>>> reinflated=new ArrayList<List<List<String>>>();
		String currindex=flattened.get(0).get(ncols-1);
		List<List<String>> subtable=new ArrayList<List<String>>();
		subtable.add(flattened.get(0).subList(0,ncols-1));
		for(int i=1;i<flattened.size();i++){
			String tindex=flattened.get(i).get(ncols-1);
			if(tindex.equals(currindex)){
				subtable.add(flattened.get(i).subList(0,ncols-1));
			} else {
				currindex=tindex;
				reinflated.add(subtable);
				subtable=new ArrayList<List<String>>();
				subtable.add(flattened.get(i).subList(0,ncols-1));
			}
		}
		reinflated.add(subtable);
		return reinflated;
	}

	public static void main(String[] args){
		//merges all of the merged_counted lanes in args
		if(args.length<2){
			System.out.println("need more than one file to merge");
		}
		String dir=(new File(args[0])).getParent();
		Object[] lane1obj=LaneMerger.getInflatedFile(args[0]);
		String[] collabels=(String[])lane1obj[1];
		List<List<List<String>>> lane1inf=(List<List<List<String>>>)lane1obj[0];
		for(int i=1;i<args.length;i++){
			Object[] lanexobj=LaneMerger.getInflatedFile(args[i]);
			List<List<List<String>>> lanexinf=(List<List<List<String>>>)lanexobj[0];
			(new LaneMerger()).mergeLanes(lane1inf,lanexinf);
		}
		List<List<String>> flatGenes=GeneMapper.flattenTables(lane1inf,collabels);
		List<String> head=flatGenes.get(0);
		flatGenes.remove(0);
		table_tools2.sort_listtable(flatGenes,3);
		flatGenes.add(0,head);
		System.out.println("writing output to merged_mapped_counted.csv");
		table_tools2.writeTableToFile(dir+File.separator+"merged_mapped_counted.csv",flatGenes,1);
		System.out.println("finished");
	}

}
