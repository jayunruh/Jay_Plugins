package jalgs.jseq;

import jalgs.jsim.rngs;

public class sim_sequences{
	//here we simulate DNA or protein sequences
	//for DNA, 0,1,2,3 are equivalent to ATCG
	//for protein, 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19
	//are equivalent to AVLIMFYW CGP STNQ RHK DE
	rngs random;
	int basisSize;
	public static String nucnames="ATCG";
	public static String aanames="AVLIMFYWCGPSTNQRHKDE";
	
	public sim_sequences(int basisSize) {
		//basis size is 4 for oligos and 20 for peptides
		this.basisSize=basisSize;
		random=new rngs();
	}
	
	/*********************
	 * here we simulate a random oligo or peptide (depending on basisSize)
	 * @param avglength: average length
	 * @param lengthstdev: length stdev
	 * @param maxlen: max sequence length
	 * @param padflag: 0 for zeros, 1 for random, 2 for mirror
	 * @return: peptide array and length
	 */
	public Object[] simNegative(float avglength,float lengthstdev,int maxlen,int padflag) {
		//the padflag is 0 for zeros, 1 for random, and 2 for mirrored
		int seqlen=(int)Math.round(random.gasdev(avglength,lengthstdev));
		int minlen=(int)Math.round(avglength-2.0f*lengthstdev);
		if(seqlen<minlen) seqlen=minlen;
		if(seqlen>maxlen) seqlen=maxlen;
		byte[] peparr=new byte[maxlen];
		int offset=(maxlen-seqlen)/2;
		for(int i=0;i<seqlen;i++) {
			peparr[i+offset]=(byte)random.unidev((double)basisSize,0.0);
		}
		//now do the padding
		if(padflag==1) {
			//this is random padding
			for(int i=0;i<offset;i++) {
				peparr[i]=(byte)random.unidev((double)basisSize,0.0);
				peparr[maxlen-i-1]=(byte)random.unidev((double)basisSize,0.0);
			}
		}
		if(padflag==2) {
			//this is mirror image padding
			for(int i=0;i<offset;i++) {
				peparr[i]=peparr[2*offset-i];
				peparr[maxlen-i-1]=peparr[seqlen+i];
			}
		}
		return new Object[] {peparr,Integer.valueOf(seqlen)};
	}
	
	public byte[] padSequence(byte[] seq,int padlen,int padflag){
		byte[] peparr=new byte[padlen];
		int offset=(padlen-seq.length)/2;
		System.arraycopy(seq,0,peparr,offset,seq.length);
		if(padflag==1){
			//this is random padding
			for(int i=0;i<offset;i++) {
				peparr[i]=(byte)random.unidev((double)basisSize,0.0);
				peparr[padlen-i-1]=(byte)random.unidev((double)basisSize,0.0);
			}
		}
		if(padflag==2) {
			//this is mirror image padding
			for(int i=0;i<offset;i++) {
				peparr[i]=peparr[2*offset-i];
				peparr[padlen-i-1]=peparr[seq.length+i];
			}
		}
		return peparr;
	}
	
	/***************
	 * here we simulate many sequences
	 * @param avglength
	 * @param lengthstdev
	 * @param nsims
	 * @return
	 */
	public Object[] simNegatives(float avglength,float lengthstdev,int padflag,int nsims){
		int maxlen=(int)Math.round(avglength+2.0f*lengthstdev);
		byte[][] peparr=new byte[nsims][maxlen];
		int[] lengths=new int[nsims];
		for(int i=0;i<nsims;i++) {
			Object[] temp=simNegative(avglength,lengthstdev,maxlen,padflag);
			peparr[i]=(byte[])temp[0];
			lengths[i]=((Integer)temp[1]).intValue();
		}
		return new Object[] {peparr,lengths};
	}
	
	/**********************
	 * here we simulate many sequences with positive signals inserted
	 * @param avglength
	 * @param lengthstdev
	 * @param padflag
	 * @param nsims
	 * @param signal
	 * @param border
	 * @return
	 */
	public Object[] simPositives(float avglength,float lengthstdev,int padflag,int nsims,byte[] signal,int border) {
		Object[] temp=simNegatives(avglength,lengthstdev,padflag,nsims);
		byte[][] peparr=(byte[][])temp[0];
		int[] lengths=(int[])temp[1];
		int maxlen=peparr[0].length;
		byte[][] locations=new byte[nsims][maxlen];
		for(int i=0;i<nsims;i++) {
			//need to make sure our negative simulation didn't have accidental signals in it
			fixAccidentalSignals(peparr[i],signal);
			int maxpos=lengths[i]-border-signal.length;
			int pos=(int)random.unidev(maxpos,border);
			int offset=(maxlen-lengths[i])/2;
			for(int j=0;j<signal.length;j++) {
				peparr[i][j+offset+pos]=signal[j];
				locations[i][j+offset+pos]=(byte)1;
			}
		}
		return new Object[] {peparr,lengths,locations};
	}
	
	/****************
	 * here we simulate a mixture of positive and negative sequences sorted randomly
	 * @param avglength
	 * @param lengthstdev
	 * @param padflag
	 * @param nsims
	 * @param signal
	 * @param border
	 * @return
	 */
	public Object[] simMixture(float avglength,float lengthstdev,int padflag,int nsims,byte[] signal,int border) {
		int npos=nsims/2;
		Object[] possims=simPositives(avglength,lengthstdev,padflag,npos,signal,border);
		Object[] negsims=simNegatives(avglength,lengthstdev,padflag,nsims-npos);
		byte[][] posarr=(byte[][])possims[0];
		byte[][] posloc=(byte[][])possims[2];
		int[] poslengths=(int[])possims[1];
		byte[][] negarr=(byte[][])negsims[0];
		int maxlen=posloc[0].length;
		int[] neglengths=(int[])negsims[1];
		for(int i=0;i<(nsims-npos);i++) {
			//fix accidental signals in the negative sims
			fixAccidentalSignals(negarr[i],signal);
		}
		//get a random order for all of these
		int[] order=random.random_order(nsims);
		byte[][] combarr=new byte[nsims][];
		byte[][] combloc=new byte[nsims][];
		int[] comblengths=new int[nsims];
		for(int i=0;i<nsims;i++) {
			if(order[i]>=npos) {
				int temp=order[i]-npos;
				combarr[i]=negarr[temp];
				combloc[i]=new byte[maxlen];
				comblengths[i]=neglengths[temp];
			} else {
				int temp=order[i];
				combarr[i]=posarr[temp];
				combloc[i]=posloc[temp];
				comblengths[i]=poslengths[temp];
			}
		}
		return new Object[] {combarr,comblengths,combloc};
	}
	
	/*****************
	 * converts a simulation into an nsims x maxlength x basisSize binary array
	 * @param seqarr
	 * @return
	 */
	public byte[][][] binarizeSeq(byte[][] seqarr){
		int nsims=seqarr.length;
		int maxlen=seqarr[0].length;
		byte[][][] binarized=new byte[nsims][maxlen][basisSize];
		for(int i=0;i<nsims;i++) {
			for(int j=0;j<maxlen;j++) {
				binarized[i][j][seqarr[i][j]&0xff]=(byte)1;
			}
		}
		return binarized;
	}
	
	/****************
	 * this version keeps the last two dimensions in raster format for image purposes
	 * @param seqarr
	 * @return
	 */
	public byte[][] binarizeSeq2(byte[][] seqarr){
		int nsims=seqarr.length;
		int maxlen=seqarr[0].length;
		byte[][] binarized=new byte[nsims][maxlen*basisSize];
		for(int i=0;i<nsims;i++) {
			for(int j=0;j<maxlen;j++) {
				int id=seqarr[i][j]&0xff;
				binarized[i][j*basisSize+id]=(byte)1;
			}
		}
		return binarized;
	}
	
	public int findSignal(byte[] sim,byte[] signal){
		//here we search for the presence of the signal sequence in our simulation
		for(int i=0;i<(sim.length-signal.length);i++){
			boolean found=true;
			for(int j=0;j<signal.length;j++){
				if(sim[i+j]!=signal[j]){
					found=false;
					break;
				}
			}
			if(found){
				return i;
			}
		}
		return -1;
	}

	public void fixAccidentalSignals(byte[] sim,byte[] signal){
		//here we fix signals that appear by accident in our sequence
		int sigpos=findSignal(sim,signal);
		while(sigpos>=0){
			for(int i=0;i<signal.length;i++){
				sim[i+sigpos]=(byte)random.unidev((double)basisSize,0.0);
			}
			sigpos=findSignal(sim,signal);
		}
		return;
	}

	public int random_nucleotide(float[] probdist){
		//generates a random integer according to a probability distribution
		double temp=random.unidev(1.0,0.0);
		float[] probdistcum=new float[probdist.length];
		probdistcum[0]=probdist[0];
		for(int i=1;i<probdist.length;i++){
			probdistcum[i]=probdistcum[i-1]+probdist[i];
		}
		for(int i=0;i<probdist.length;i++){
			if(temp<=probdistcum[i]) return i;
		}
		return (probdist.length-1);
	}
	
	public static String getAlphaSeq(byte[] seq,String alphacode) {
		StringBuffer sb=new StringBuffer();
		int tbasisSize=alphacode.length();
		for(int i=0;i<seq.length;i++) {
			int index=seq[i]&0xff;
			if(index<tbasisSize && index>=0) {
				sb.append(alphacode.charAt(index));
			} else {
				return null;
			}
		}
		return sb.toString();
	}
	
	public static byte[] getByteSeq(String alphaseq,String alphacode) {
		byte[] temp=alphaseq.getBytes();
		byte[] code=new byte[temp.length];
		for(int i=0;i<temp.length;i++) {
			int index=alphacode.indexOf(temp[i]);
			if(index<0) return null;
			code[i]=(byte)index;
		}
		return code;
	}

}
