import java.util.ArrayList;
import java.util.HashMap;

public class Tests{
	
	// TODO Weir & Cockrahm Fst
	
	int caseSize, contSize, outgSize; // sample sizes (e.g., number of haplotypes)
	SnpMat matCase, matCont, matOutg; // SNP matrices
	
	/********************** C'tor ***********************/
	public Tests(int case_size, int cont_size, int outg_size, SnpMat mat_case, SnpMat mat_cont, SnpMat mat_outg){
		caseSize = case_size; contSize = cont_size; outgSize = outg_size;
		matCase = mat_case; matCont = mat_cont; matOutg = mat_outg;
	}
	
	/******************** testWindow *********************/
	public String testWindow(int caseScol, int caseEcol, int contScol, int contEcol, int outgScol, int outgEcol, int winStart, int winEnd){
		/* apply tests on given window for case, control & outgroup */
		
		// Theta-estimates & statistics
		String ans = "";
		int caseNsites = 0, contNsites = 0, outgNsites = 0;
		int caseNSitesNoFix = 0, contNSitesNoFix = 0, outgNSitesNoFix = 0;
		double caseWT = 0.0, contWT = 0.0;
		double caseH = 0.0, contH = 0.0;
		double caseTT = 0.0, contTT = 0.0, outgTT = 0.0; 
		double caseSF = 0.0, contSF = 0.0, caseTD = 0.0, contTD = 0.0;
		double caseContFst = 0.0, caseOutgFst = 0.0, contOutgFst = 0.0;
		double caseContFst_t = 0.0, caseOutgFst_t = 0.0, contOutgFst_t = 0.0;
		double case_cont_logTT = 0.0, case_cont_logSF = 0.0, PBS = 0.0, case_cont_logH = 0.0;
		
		// CASE population window Thetas & statistics
		if(caseScol != -1 && caseEcol != -1 && caseScol < caseEcol){			
			caseNsites = caseEcol - caseScol;
			caseNSitesNoFix = count_non_fixed_snps(matCase.m_freqs_est, caseScol, caseEcol);
			if(caseNSitesNoFix > 0){
				// non-empty matrix after masking fixed SNPs
				caseTD = tajimas_d(matCase.m_freqs_est, caseSize, caseScol, caseEcol, caseNSitesNoFix);
				caseTT = tajimas_theta(matCase.m_freqs_est, caseSize, caseScol, caseEcol);
				caseWT = wattersons_theta(caseSize, caseNSitesNoFix);
				caseSF = sum_frequencies_theta(matCase.m_freqs_est, caseScol, caseEcol);
				caseH = fay_and_wus(matCase.m_freqs_est, caseSize, caseScol, caseEcol);
			}
		}
		
		// CONTROL population window Thetas & statistics 		
		if(contScol != -1 && contEcol != -1 && contScol < contEcol){
			contNsites = contEcol - contScol;
			contNSitesNoFix = count_non_fixed_snps(matCont.m_freqs_est, contScol, contEcol);
			if(contNSitesNoFix > 0){
				// non-empty matrix after masking fixed SNPs
				contTD = tajimas_d(matCont.m_freqs_est, contSize, contScol, contEcol, contNSitesNoFix);	
				contTT = tajimas_theta(matCont.m_freqs_est, contSize, contScol, contEcol);
				contWT = wattersons_theta(contSize, contNSitesNoFix);
				contSF = sum_frequencies_theta(matCont.m_freqs_est, contScol, contEcol);
				contH = fay_and_wus(matCont.m_freqs_est, contSize, contScol, contEcol);
			}
		}
		
		// OUT-GROUP population window Thetas & statistics 		
		if(outgScol != -1 && outgEcol != -1 && outgScol < outgEcol){
			outgNsites = outgEcol - outgScol;
			outgNSitesNoFix = count_non_fixed_snps(matOutg.m_freqs_est, outgScol, outgEcol);
			if(outgNSitesNoFix > 0){
				// non-empty matrix after masking fixed SNPs
				outgTT = tajimas_theta(matOutg.m_freqs_est, outgSize, outgScol, outgEcol);
			}
		}
		
		// calc Fst & PBS statistics
		caseContFst = Fst(caseNsites, contNsites, caseSize, contSize, caseTT, contTT, matCase, matCont, caseScol, caseEcol, contScol, contEcol);
		if(caseContFst != 0.0 && caseContFst != 1.0) caseContFst_t = -Math.log(1.0 - caseContFst);
		
		caseOutgFst = Fst(caseNsites, outgNsites, caseSize, outgSize, caseTT, outgTT, matCase, matOutg, caseScol, caseEcol, outgScol, outgEcol);
		if(caseOutgFst != 0.0 && caseOutgFst != 1.0) caseOutgFst_t = -Math.log(1.0 - caseOutgFst);
		
		contOutgFst = Fst(contNsites, outgNsites, contSize, outgSize, contTT, outgTT, matCont, matOutg, contScol, contEcol, outgScol, outgEcol);
		if(contOutgFst != 0.0 && caseOutgFst != 1.0) contOutgFst_t = -Math.log(1.0 - contOutgFst);
		
		if(matOutg != null)
			PBS = (caseContFst_t + caseOutgFst_t - contOutgFst_t) / 2.0;
		
		// calc LR statistics
		case_cont_logSF = Math.log( (contSF + 0.1) / (caseSF + 0.1) );
		case_cont_logTT = Math.log( (contTT + 0.1) / (caseTT + 0.1) );
		case_cont_logH = Math.log( (contH + 0.1) / (caseH + 0.1) );
		
		// window information
		ans += matCase.m_chr + "\t" + winStart + "\t" +  winEnd + "\t" + caseNsites + "\t" + caseNSitesNoFix + "\t" + 
			   contNsites + "\t" + contNSitesNoFix + "\t" + outgNsites + "\t" + outgNSitesNoFix + "\t";
		
		// window statistics
		ans += caseWT + "\t" + contWT + "\t" + caseTT + "\t" + contTT + "\t" + case_cont_logTT + "\t" + 
			   caseTD + "\t" + contTD + "\t" + caseSF + "\t" + contSF + "\t" + case_cont_logSF + "\t" + 
			   caseContFst_t + "\t" + contOutgFst_t + "\t" + caseOutgFst_t + "\t" + PBS + "\t" + case_cont_logH + "\n";
		
		return ans;
	}
	
	/*********************** Fst ************************/
	public static double Fst(int nSites1, int nSites2, int size1, int size2, double tt1, double tt2, 
						   	 SnpMat mat1, SnpMat mat2, int sCol1, int eCol1, int sCol2, int eCol2){
		
		double Fst = 0.0;
		if(nSites1 > 0 && nSites2 > 0){
			double pi_within = ( (tt1 * size1 * (size1-1) ) + (tt2 * size2 * (size2-1)) ) / (double)( size1*(size1-1) + size2*(size2-1) );
			double pi_between = piBetween(mat1, sCol1, eCol1, mat2, sCol2, eCol2);
			double pi_total   = piTotal(size1, size2, tt1, tt2, pi_between);
			
			Fst =  1.0 - ( pi_within / pi_total );
			if(Fst < 0.0 || Double.isNaN(Fst) || Double.isInfinite(Fst)) Fst = 0.0; // not cheating!
			
		}
		
		return Fst;
	}
	
	/********************* piTotal **********************/
	public static double piTotal(int nn1, int nn2, double within1, double within2, double between) {
		/* returns pi-Total given pi-Within of two populations */
		
		double n1 = (double)nn1, n2 = (double)nn2;
		double totalPai = ( (2.0*n1*n2 * between) + (n1 * (n1 - 1.0) * within1) + (n2 * (n2 - 1.0) * within2) ) /
						  ( (n1 + n2) * (n1 + n2 - 1.0) );
		return totalPai;
	}

	/*************** count_non_fixed_snps ****************/
	public static int count_non_fixed_snps(ArrayList<Double> snpFreq, int startCol, int endCol){
		/* returns the number of non-fixed SNPs in a window */
		int non_fixed = 0;
		for (int snp = startCol; snp < endCol; snp++)
			if(snpFreq.get(snp) != 1.0) non_fixed++;
		return non_fixed;
	}
	
	/******************** piBetween *********************/
	public static double piBetween(SnpMat mat1, int startCol1, int endCol1, SnpMat mat2, int startCol2, int endCol2){
		
		/* calc average heterozygosity between case & control populations */
		double piBetween = 0.0;
		
		// get SNP to column mapping for case & control windows 
		HashMap<Integer,Integer> caseSNPs2cols = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> contSNPs2cols = new HashMap<Integer,Integer>();
		for(int i=startCol1; i<endCol1; i++)
			caseSNPs2cols.put(mat1.m_col2pos.get(i), i);
		for(int i=startCol2; i<endCol2; i++)
			contSNPs2cols.put(mat2.m_col2pos.get(i), i);
		
		// add between heterozygosity of SNPs in case (unique + shared)
		for(Integer caseSNP : caseSNPs2cols.keySet()){
			double fcase = mat1.m_freqs_est.get(caseSNPs2cols.get(caseSNP));
			
			if(contSNPs2cols.containsKey(caseSNP)){ 
				double fcont = mat2.m_freqs_est.get(contSNPs2cols.get(caseSNP));
				piBetween += fcase * (1.0 - fcont) + fcont * (1.0 - fcase); // SNP in case & control
			}
			else
				piBetween += fcase; // SNP only in case
		}
		
		// add between heterozygosity of SNPs unique to control
		for(Integer contSNP : contSNPs2cols.keySet()){
			double fcont = mat2.m_freqs_est.get(contSNPs2cols.get(contSNP));
			if(!caseSNPs2cols.containsKey(contSNP))
				piBetween += fcont; // SNP only in control
		}
		
		return piBetween;
	}

	/*************** sum_frequencies_theta ****************/	
	public static double sum_frequencies_theta(ArrayList<Double> snp_freq, int start_col, int end_col) {
		/* return the sum of SNP frequencies in the given matrix window */
		double sumF = 0.0;
		for (int i=start_col; i < end_col; i++){
			if(snp_freq.get(i) != 1.0)
				sumF += snp_freq.get(i); // count SNP only if not fixed	
		}
		return sumF;
	}

	/***************** wattersons_theta ******************/
	public static double wattersons_theta(int nSamp, int nSites) {
		/* calculate Watterson's theta */
		double num_sites = (double)nSites;
		double harmonic = 0.0;
		for(int i = 1; i < nSamp; i++) 
			harmonic += 1 /(double)i; 
		
		return num_sites / harmonic;
	}
	
	/******************* tajimas_theta ********************/
	public static double tajimas_theta(ArrayList<Double> freqs, int numSamp, int startCol, int endCol) {
		/* calculate average heterozygosity from SNP frequencies */
		double avg_heterozigosity = 0.0;
		for(int i = startCol; i < endCol; i++){
			double f = freqs.get(i);
			if(f < 1.0) // only non-fixed SNPs 
				avg_heterozigosity += f * (1 - f);
		}
		double scaling_factor = 2.0 * ( (double)numSamp / ( (double)numSamp - 1.0 ) );
		return scaling_factor * avg_heterozigosity;
	}

	/******************* tajimas_theta ********************/
	public static double fay_and_wus(ArrayList<Double> freqs, int numSamp, int startCol, int endCol) {
		/* calculate Fay & Wu's theta from SNP frequencies */
		double theta = 0.0;
		for(int i = startCol; i < endCol; i++){
			double f = freqs.get(i);
			if(f < 1.0) // only non-fixed SNPs 
				theta += f * f;
		}
		double scaling_factor = 2.0 * ( (double)numSamp / ( (double)numSamp - 1.0 ) );
		return scaling_factor * theta;
	}
	
	/********************* tajimas_d **********************/
	public static double tajimas_d(ArrayList<Double> freqs, int nSamp, int startCol, int endCol, int numNonFixedSites) {
		/* calculate Tajima's D (k^ - S/a1 ) / standard deviation */

		// get Tajima's theta
		double thetaTaj = tajimas_theta(freqs, nSamp, startCol, endCol);
		
		// get Watterson's theta (not via method, to get Squared Harmonic number)
		double numSites = (double)numNonFixedSites;
		double harmonic = 0.0;
		double square_harmonic = 0.0;
		for(int i=1; i < nSamp; i++){
			harmonic += 1/(double)i;
			square_harmonic += 1/(double)(i*i); 
		}
		double thetaWaterson = numSites / harmonic;
		if(thetaTaj - thetaWaterson == 0) 
			return 0.0; // Tajima's D is 0 (no need for std-dev)
		
		// calculate standard-deviation (from variance) of difference
		double n = (double)nSamp;
		double a1 = harmonic;
		double b1 = (n + 1) / ( 3*(n -1) );
		double c1 = b1 - (1 / a1);
		double e1 = c1 / a1;
		double a2 = square_harmonic;
		double b2 = ( 2 * (n*n + n + 3) ) / ( 9 * n * (n-1) );
		double c2 = b2 - (n+2)/(a1*n) + a2/(a1*a1);
		double e2 = c2 / ( (a1*a1) + a2 );
		double var = (e1 * numSites) + (e2 * numSites * (numSites - 1) );
		double stddev = Math.sqrt(var);
		
		return (thetaTaj - thetaWaterson) / stddev; // (non-zero) Tajima's D
	
	}
}
