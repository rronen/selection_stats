import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class SnpMat {
	
	/****************** class variables ******************/

	protected File m_snp_file;
	protected ArrayList<ArrayList<Integer>> m_mat;
	protected HashMap<Integer, String> m_row2SampNames;
	protected HashMap<Integer, Integer> m_col2pos;
	protected HashMap<Integer, Double> m_col2pos_d; // only for mpop files
	protected ArrayList<Double> m_freqs_est;
	protected ArrayList<Integer> m_freqs_abs;
	protected HashMap<Integer, Integer> m_freq_spect, m_scaled_freq_spect; // bin these
	protected int m_maxPos, m_minPos;
	protected int m_numSites, m_numHaps;
	protected String m_chr;
	
	protected boolean m_is_mpop; // true if read mpop format
	protected boolean m_skip_MA_snps; // ignore multi allelic SNPs
	protected int m_lastWinStartCol;
	protected int m_lastWinEndCol;
	
	protected int window_size;
	
	/******************** Constructor ********************/
	public SnpMat(File vcf, boolean skipMAsnps, int win_size){
		m_snp_file = vcf;
		m_skip_MA_snps = skipMAsnps;
		window_size = win_size;
		
		m_numHaps = 0; m_numSites = 0;
		m_maxPos = -1; m_minPos = -1;
		m_lastWinStartCol = 0; m_lastWinEndCol = 0;
		
		m_freqs_est = new ArrayList<Double>();
		m_freqs_abs = new ArrayList<Integer>();
		m_freq_spect = new HashMap<Integer, Integer>();
		m_scaled_freq_spect = new HashMap<Integer, Integer>();
		m_is_mpop = false;
	}
	
	/**************** first_col_in_window ****************/
	public int first_col_in_window(int winStart) {
		// Returns first column with SNP position >= winStart, or -1 if none exists
		// Note: returned column (SNP) may also have position > winEnd, in which case winStart == winEnd
		int check = m_lastWinStartCol;
		while(check < m_numSites){
			if(m_col2pos.get(check) >= winStart){
				break;
			}
			check++;
		}
		m_lastWinStartCol = check;
		if(check == m_numSites || m_col2pos.get(check) < winStart)
			return -1;
		return check;  
	}

	/***************** last_col_in_window *****************/
	public int last_col_in_window(int winEnd) {
		// Returns column to the right of last SNP with location <= winEnd, or -1 if none exists
		int check = m_lastWinEndCol;
		while(check < m_numSites){
			if(m_col2pos.get(check) <= winEnd){
				check++;
				continue;
			}
			break;
		}
		m_lastWinEndCol = check;
		if(check == 0 || m_col2pos.get(check - 1) > winEnd)
			return -1;
		return check;
	}
	
	/*************** file2mat_n_est_freqs ****************/
	public int file2mat_n_est_freqs(PosFilter posFlt, boolean write_freqs, String freq_out){
		
		// create proper population reader and read input file
		Reader pr; 
		if(is_mpop(m_snp_file)){
			pr = new ReaderMpop(m_snp_file, window_size);
			m_is_mpop = true;
		}
		else
			pr = new ReaderVCF(m_snp_file, m_skip_MA_snps, posFlt);
		
		// get info from reader
		pr.read();
		m_mat = pr.get_mat();
		m_numSites = pr.get_num_sites();
		m_numHaps = pr.get_num_haps();
		m_col2pos = pr.get_col2pos();
		if(m_is_mpop)
			m_col2pos_d = pr.get_col2pos_d();
		
		m_chr = pr.get_chr();
		m_maxPos = pr.get_max_pos(); m_minPos = pr.get_min_pos();
		m_row2SampNames = pr.get_row2sampNames();
		
		// estimate frequencies
		est_f_fspectra_n_flt_nonvars(write_freqs, freq_out);
		
		// return number of haplotypes
		return m_numHaps;
	}

	/*********** est_f_fspectra_n_flt_nonvars ***********/
	private void est_f_fspectra_n_flt_nonvars(boolean output, String outfile){
		
		//******* estimate freq & mark non-variants (100% major allele!) to remove *******//
		int flt_count = 0, kept_count = 0;
		ArrayList<Integer> to_remove = new ArrayList<Integer>();
		HashMap<Integer, Integer> col2pos_new   = new HashMap<Integer, Integer>();
		HashMap<Integer, Double>  col2pos_d_new = new HashMap<Integer, Double>(); // only for mpop files
		
		for (int site=0; site < m_numSites; site++){
			
			double tot_ma = 0.0, tot_observed = 0.0; // accounts for missing calls 
			
			for (int hap=0; hap < m_numHaps; hap++){
				if( m_mat.get(hap).get(site) == 1) tot_ma++;
				if( m_mat.get(hap).get(site) >= 0) tot_observed++;
			}
			
			// frequency estimate with #-observed-haplotypes-at-site as denominator
			double maf_est = 0.0;
			if(tot_observed > 0.0)
				maf_est = tot_ma / tot_observed;
			
			if(maf_est == 0.0){
				to_remove.add(site); // remove nonvariant (100% major allele) site
				flt_count++;
			}
			else{
				// keep site - variant
				m_freqs_est.add(kept_count, maf_est);
				m_freqs_abs.add(kept_count, (int)tot_ma); // #-occurrences in population (absolute frequency)
				col2pos_new.put(kept_count, m_col2pos.get(site));
				if(m_is_mpop)
					col2pos_d_new.put(kept_count, m_col2pos_d.get(site));
				kept_count++;
			}
		}
		
		// update state & write frequencies
		m_numSites = kept_count;
		m_col2pos = col2pos_new;
		m_col2pos_d = col2pos_d_new;
		
		//****** remove (last to first) non variant sites from matrix ******//
		Collections.reverse(to_remove);
		for(int site : to_remove){
			for (int hap=0; hap < m_numHaps; hap++){
				m_mat.get(hap).remove(site);
			}
		}
		
		//******************* generate frequency spectrum *******************//
		// calculate exact site frequency spectra (SSF) with one bin per haplotype
		// unlike the frequency estimates, assumes all haplotypes were observed at each site
		for(int i=1; i <= m_numHaps; i++)
			m_freq_spect.put(i, 0);
		
		for (int site=0; site < m_numSites; site++){
			int f = m_freqs_abs.get(site); // absolute frequency
			m_freq_spect.put(f, m_freq_spect.get(f) + 1 );
		}
		
		//**************** generate scaled frequency spectrum ****************//
		// exact scaled site frequency spectra (sSFS) with one bin per haplotype
		// unlike the frequency estimates, assumes all haplotypes were observed at each site
		for(int i=1; i <= m_numHaps; i++)
			m_scaled_freq_spect.put(i, 0);
		
		for (int site=0; site < m_numSites; site++){
			int f = m_freqs_abs.get(site); // absolute frequency
			m_scaled_freq_spect.put(f, m_scaled_freq_spect.get(f) + f );
		}
		
		// output frequencies
		if(output) write_freq(outfile);
		System.out.println("[SnpMat]::Estimated Freq & SFS (removed " + flt_count + " non-variant sites)");
	}
	
	/********************* write_freq *********************/
	private void write_freq(String outFname) {
		
		// make file name
		File parentDir = new File(outFname).getAbsoluteFile().getParentFile();
		String name_freq = m_snp_file.getName();
		if(name_freq.endsWith(".vcf") || name_freq.endsWith(".VCF")){ 
			name_freq = name_freq.substring(0, name_freq.length() - 4);
		}
		name_freq += ".freq";
		String freq_path = parentDir.toString() + File.separator + name_freq;
		
		// write frequencies
		try{	
			PrintWriter out_w = new PrintWriter(new BufferedWriter(new FileWriter(freq_path)));
			out_w.print("#nHaps\t" + m_numHaps + "\n");
			for (int site_col=0; site_col < m_numSites; site_col++){
				String to_write = m_col2pos.get(site_col) + "\t" + m_freqs_est.get(site_col);
				if(m_is_mpop)
					to_write += "\t" + m_col2pos_d.get(site_col);
				to_write += "\n";
				out_w.print(to_write);
			}
			out_w.close();
		}
		catch(IOException e){ e.printStackTrace(); System.exit(1); }		
	}

	/********************* is_mpop **********************/
	private boolean is_mpop(File file){
		boolean ans = false;
		// read first line to determine if MPOP format
		try {
			 BufferedReader mpopFile = new BufferedReader(new FileReader(file));
			 String line = mpopFile.readLine();
			 if(line.startsWith("mpop")) ans = true;
	         mpopFile.close();
		}
		catch (IOException e){
			e.printStackTrace();
			Utils.exit_print("[SnpMat]::problem reading input file: " + file.toString());
		}
		
		return ans;
	}
	
}
