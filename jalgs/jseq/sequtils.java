package jalgs.jseq;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.DNASequence;

public class sequtils{
	//this class has a bunch of static methods to deal with DNA and protein sequences
	public static String nucnames="ATCG";
	public static String aanames="AVLIMFYWCGPSTNQRHKDE";
	public static String start="ATG";
	public static String amber="TAG";
	public static String opal="TGA";
	public static String ochre="TAA";
	public static String[][] codons={{"A","GCT","GCC","GCA","GCG"},{"V","GTT","GTC","GTA","GTG"},{"L","CTT","CTC","CTA","CTG","TTA","TTG"},
		{"I","ATT","ATC","ATA"},{"M","ATG"},{"F","TTT","TTC"},{"Y","TAT","TAC"},{"W","TGG"},{"C","TGT","TGC"},{"G","GGT","GGC","GGA","GGG"},
		{"P","CCT","CCC","CCA","CCG"},{"S","TCT","TCC","TCA","TCG","AGT","AGC"},{"T","ACG","ACA","ACC","ACT"},{"N","AAC","AAT"},{"Q","CAG","CAA"},
		{"R","AGG","AGA","CGT","CGC","CGA","CGG"},{"H","CAT","CAC"},{"K","AAG","AAA"},{"D","GAC","GAT"},{"E","GAG","GAA"},{"*",amber,opal,ochre}};

	/********************
	 * this function finds orfs longer than minlen in seq
	 * @param seq
	 * @param minlen
	 * @return
	 */
	public static int[][] orfFinder(String seq1,int minlen){
		//make sure all of our sequence is in DNA format and upper case
		String seq=makeDNA(seq1.toUpperCase());
    	int seqlen=seq.length();
    	int maxorfs=(int)(seqlen/minlen);
    	//start by searching for start codons
    	int tempstart=seq.indexOf("ATG");
    	int[][] orfbounds =new int[maxorfs][];
    	int norfs=0;
    	while(tempstart>0 && tempstart <(seqlen-minlen-1)){ 
    		//we found a valid start, search for the nearest stop 
    		int stop=findStop(seq,tempstart); 
    		if(stop>=0){
    			//if the product is long enough, record it
    			if((stop-tempstart)>(minlen-3)){
    				orfbounds[norfs]=new int[]{tempstart,stop};
    				norfs++;
    			}
    		} else {
    			//if no products are found, break
    			break;
    		}
    		tempstart=seq.indexOf("ATG",stop+3);
    	}
    	if(norfs>0) return subarray(orfbounds,norfs); 
    	else return null; 
	}
	
	public static String makeDNA(String seq){
		//here we convert all U's to T's
		return seq.replace('U','T');
	}

	/***********************
	 * this function searches for the nearest stop codon in frame with a start codon
	 * @param seq
	 * @param start
	 * @param minlen
	 * @return
	 */
	public static int findStop(String seq,int start){
		//search 3 at a time for stop codons
		int pos=start;
		int seqlen=seq.length();
		int maxpos=seqlen-4;
		while(pos<=maxpos && !isStop(seq.substring(pos,pos+3))){
			pos+=3;
		}
		if(pos<=maxpos) return pos;
		else return -1;
	}
	
	public static boolean isStop(String codon){
		if(codon.equals(amber)) return true;
		if(codon.equals(opal)) return true;
		if(codon.equals(ochre)) return true;
		return false;
	}

	public static int[][] subarray(int[][] arr,int len){
		int[][] newarr =new int[len][]; 
		for(int i=0;i<len;i++){
			newarr[i]=arr[i]; 
		} 
		return newarr; 
	}
	
	public static String makeComplement(String seq1){
		String seq=makeDNA(seq1.toUpperCase());
    	char[] seq2=seq.toCharArray();
    	for(int i=0;i<seq2.length;i++){
    		if(seq2[i]=='A') seq2[i]='T';
    		else if(seq2[i]=='T') seq2[i]='A';
    		else if(seq2[i]=='G') seq2[i]='C';
    		else if(seq2[i]=='C') seq2[i]='G';
    	}
    	return new String(seq2);
	}
	
	public static String reverse(String seq){
		char[] seq2=seq.toCharArray();
		char[] rev=new char[seq2.length];
		for(int i=0;i<seq2.length;i++){
			rev[seq2.length-1-i]=seq2[i];
		}
		return new String(rev);
	}
	
	public static Object[] translate(String seq){
		//translates an orf and creates a codon usage table
		//need to handle N values
		int naas=aanames.length();
		int[][] codonfreq=new int[naas+1][6];
		int ncodons=seq.length()/3;
		StringBuffer translation=new StringBuffer();
		for(int i=0;i<ncodons;i++){
			int[] searchres=findCodon(seq.substring(i*3,i*3+3));
			if(searchres!=null){
				translation.append(codons[searchres[0]][0]);
				codonfreq[searchres[0]][searchres[1]]++;
			} else {
				translation.append("*");
			}
		}
		return new Object[]{translation.toString(),codonfreq};
	}
	
	public static int[] findCodon(String query){
		for(int i=0;i<codons.length;i++){
			for(int j=1;j<codons[i].length;j++){
				if(query.equals(codons[i][j])){
					return new int[]{i,j-1};
				}
			}
		}
		return null;
	}
	
	public static Object[] pairwiseDNA(String seq1,String seq2,boolean global){
		try{
			DNASequence s1=new DNASequence(seq1);
			DNASequence s2=new DNASequence(seq2);
			SubstitutionMatrix<NucleotideCompound> matrix=SubstitutionMatrixHelper.getNuc4_4();
			PairwiseSequenceAlignerType atype=PairwiseSequenceAlignerType.LOCAL;
			if(global) atype=PairwiseSequenceAlignerType.GLOBAL;
			SequencePair<DNASequence,NucleotideCompound> pair=Alignments.getPairwiseAlignment(s1,s2,atype,new SimpleGapPenalty(),matrix);
			double fraction=pair.getPercentageOfIdentity(true);
			return new Object[]{pair.toString(),fraction};
		}catch(CompoundNotFoundException e){
			e.printStackTrace();
		}
		return null;
	}
	
	public static Object[] pairwiseProtein(String seq1,String seq2,boolean global){
		//ProteinSequence s1=new ProteinSequence(seq1);
		try{
			ProteinSequence s1=new ProteinSequence(seq1);
			ProteinSequence s2=new ProteinSequence(seq2);
			SubstitutionMatrix<AminoAcidCompound> matrix=SubstitutionMatrixHelper.getBlosum62();
			PairwiseSequenceAlignerType atype=PairwiseSequenceAlignerType.LOCAL;
			if(global) atype=PairwiseSequenceAlignerType.GLOBAL;
			SequencePair<ProteinSequence,AminoAcidCompound> pair=Alignments.getPairwiseAlignment(s1,s2,atype,new SimpleGapPenalty(),matrix);
			double fraction=pair.getPercentageOfIdentity(true);
			return new Object[]{pair.toString(),fraction};
		}catch(CompoundNotFoundException e){
			e.printStackTrace();
		}
		return null;
	}

}
