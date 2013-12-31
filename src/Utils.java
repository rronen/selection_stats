import java.util.ArrayList;
import java.util.HashMap;

public class Utils{
	
	/************** avgDiffBetweenPairwise ***************/
	public double avgDiffBetweenPairwise(SnpMat caseMat, int caseStartCol, int caseEndCol, int caseNumSamp,
								  		 SnpMat contMat, int contStartCol, int contEndCol, int contNumSamp){
		/* calc average pairwise differences between case & control populations */
		
		ArrayList<ArrayList<Integer>> matCase = caseMat.m_mat ;
		ArrayList<ArrayList<Integer>> matCont = contMat.m_mat ;
		
		// get SNP to column mapping for case & control windows 
		HashMap<Integer,Integer> caseSNPs2cols = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> contSNPs2cols = new HashMap<Integer,Integer>();
		for(int i=caseStartCol; i<caseEndCol; i++)
			caseSNPs2cols.put(caseMat.m_col2pos.get(i), i);
		for(int i=contStartCol; i<contEndCol; i++)
			contSNPs2cols.put(contMat.m_col2pos.get(i), i);
		
		// get column indices of SNPs shared by case & control
		HashMap<Integer,Integer> commonSNPs2colsCase = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> commonSNPs2colsCont = new HashMap<Integer,Integer>();
		for(Integer caseSNP : caseSNPs2cols.keySet()){
			if(contSNPs2cols.containsKey(caseSNP)){ 
				commonSNPs2colsCase.put(caseSNP, caseSNPs2cols.get(caseSNP));
				commonSNPs2colsCont.put(caseSNP, contSNPs2cols.get(caseSNP));
			}
		}
		
		// finally, get average number of pairwise differences
		double numPairs = 0.0;
		double sumDiffs = 0.0;
		for(int caseN=0; caseN < caseNumSamp; caseN++){
			for(int contN=0; contN < contNumSamp; contN++){
				numPairs += 1.0;
				
				// add common SNP differences
				for(Integer commonSNP : commonSNPs2colsCase.keySet()){
					int caseCol = commonSNPs2colsCase.get(commonSNP);
					int contCol = commonSNPs2colsCont.get(commonSNP);
					if(matCase.get(caseN).get(caseCol) != matCont.get(contN).get(contCol))
						sumDiffs += 1.0;
				}
				
				// add unique SNP differences from case
				for(int caseCol=caseStartCol; caseCol<caseEndCol; caseCol++)
					sumDiffs += (double)matCase.get(caseN).get(caseCol);
				
				// add unique SNP differences from control
				for(int contCol=contStartCol; contCol<contEndCol; contCol++)
					sumDiffs += (double)matCont.get(contN).get(contCol);
				
			}
		}
		return sumDiffs / numPairs;
	}
	
	/***************** mask_fixed_snps *******************/	
	public ArrayList<ArrayList<Integer>> mask_fixed_snps(ArrayList<ArrayList<Integer>> mat, ArrayList<Double> snpFreq, int startCol, int endCol) {
		/* create a sub matrix over the desired columns, and remove fixed SNPs */
		
		// initialized masked matrix
		ArrayList<ArrayList<Integer>> subMatNoFixed = new ArrayList<ArrayList<Integer>>();
		for (int j = 0; j < mat.size(); j++)
			subMatNoFixed.add(j, new ArrayList<Integer>());
		
		
		// fill'er up!
		int non_fixed_count = 0;
		for (int snp = startCol; snp < endCol; snp++){
			if(snpFreq.get(snp) != 1.0){
				// add SNP-column only if not fixed
				for (int hap = 0; hap < mat.size(); hap++){	
					subMatNoFixed.get(hap).add(non_fixed_count, mat.get(hap).get(snp)); 	
				}
				non_fixed_count++;
			}
		}
		return subMatNoFixed;
	}
	
	/******************* tajimas_theta ********************/
	public double tajimas_theta(ArrayList<ArrayList<Integer>> mat, int numSamp) {
		/* calculate average number of pair-wise differences on SNP matrix */
		double sum_dist = 0.0;
		double num_pairs = 0.0;
		for(int i=0; i<numSamp-1; i++){
			for(int j=i+1; j<numSamp; j++){
				sum_dist += hamming_dist(mat.get(i), mat.get(j));
				num_pairs += 1.0;
			}
		}
		return sum_dist / num_pairs;
	}
	
	/******************* hamming_dist *******************/
	private double hamming_dist(ArrayList<Integer> list1, ArrayList<Integer> list2) {
		// check input sane
		if(list1.size() != list2.size()){
			System.out.println("[hamming_dist] Lists of unequal size. Quitting...");
			System.exit(1);
		}
		
		// get Hamming distance
		double dist = 0;
		for(int i = 0; i < list1.size(); i++){
			if(list1.get(i) != list2.get(i))
				dist += 1.0;
		}
		return dist;
	}
	
	/******************** exit_print ********************/
	public static void exit_print(String s){
		System.out.println("\n" + s + "\n");
		System.exit(1);
	}
	
	/*********************** usage ***********************/
	public static void exit_usage(){
		System.out.println();
		System.out.println("USAGE:");
		System.out.println("\tjava -Xmx4g -jar Selection.jar [OPTIONS...] -case case.vcf -cont cont.vcf\n");
		System.out.println("INPUT:");
		System.out.println("\tcase.vcf         case population variants (VCF)");
		System.out.println("\tcont.vcf         control population variants (VCF)");
		System.out.println();
		System.out.println("OPTIONS:");
		System.out.println("\t-o outfile       output file name");
		System.out.println("\t-outg og.vcf     out-group population variants (VCF), for PBS statistic");
		System.out.println("\t-b reg.bed       regions (e.g. exons) to perform tests in (BED), specify -f");
		System.out.println("\t-f type          filter, use \"full\" for many short regions (e.g. exons) or \"region\" for long regions (e.g. genes)");
		System.out.println();
		System.exit(1);
	}
	
} // class
