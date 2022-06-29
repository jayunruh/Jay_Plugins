package jalgs.jseq;

import jalgs.table_tools2;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class GTFReader{

	public List<Object[]> getGTFFile(String path){
		//gtf files are tab delimited with the following fields:
		//0: chr, 1: source (ensemble), 2: type, 3: start, 4: end, 5: score, 6:strand (- or +), 7: frame (0, 1, or 2), 8: semicolon delimited attributes
		//types are stop_codon, exon, transcript (spliced?), gene (everything?), CDS (no stop or utr), five_prime_UTR, three_prime_UTR
		//the attributes are just key value pairs separated by a space
		//not all types have all attributes
		//some example attributes: gene_id, gene_version, transcript_id, gene_source, gene_biotype, transcript_source, transcript_biotype, exon_id, exon_version, tss_id, p_id
		List<List<String>> listtable=table_tools2.getTableFromFile((new File(path)),"\t",false);
		String[] collabels={"scaffold","source","type","start","end","score","strand","frame"};
		List<Object[]> outtable=new ArrayList<Object[]>();
		for(int i=0;i<listtable.size();i++){
			List<String> row=listtable.get(i);
			String[][] decoded=decodeAttr(row.get(8));
			row.remove(8);
			String[] temp=row.toArray(new String[0]);
			Object[] temp2={temp,decoded};
			outtable.add(temp2);
		}
		return outtable;
	}
	
	public String[][] decodeAttr(String attr){
		//decodes the semicolon deliminated key value attribute pairs
		String[] pairs=attr.split(";");
		String[][] pairs2=new String[pairs.length][2];
		for(int i=0;i<pairs.length;i++){
			//pair might start with a space
			int firstspace=pairs[i].indexOf(" ",1);
			pairs2[i][0]=pairs[i].substring(0,firstspace).trim();
			pairs2[i][1]=pairs[i].substring(firstspace+1,pairs[i].length()).trim();
		}
		return pairs2;
	}

}
