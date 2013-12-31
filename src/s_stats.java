import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class s_stats {
	
	protected boolean m_skipMultiAsnps = true; // ignore multi allelic SNPs
	protected boolean m_writeFreqs = true; // write estimated SNP-frequencies 
	protected int m_windowSize = 50000, m_stepSize = 2000;
	protected double m_pseudoCount = 0.1; // pseudo-count for empty windows
	protected int m_windows;
	
	protected SnpMat m_matCase, m_matCont, m_matOutg;
	protected Tests tests;
	protected int m_caseSize, m_contSize, m_outgSize;
	protected PosFilter m_posFlt; // positional filter (e.g. exome)
	protected String m_outFile;
	protected Long m_startTime;
	
	/*********************** MAIN ************************/
	public static void main(String[] args) {	
		new s_stats(args).go();
	}
	
	/************************ GO *************************/
	private void go(){
		System.out.println();
		
		// read files, create SNP matrices
		m_posFlt.build_filter();
		m_caseSize = m_matCase.file2mat_n_est_freqs(m_posFlt, m_writeFreqs, m_outFile);
		m_contSize = m_matCont.file2mat_n_est_freqs(m_posFlt, m_writeFreqs, m_outFile);
		if(m_matOutg != null) m_outgSize = m_matOutg.file2mat_n_est_freqs(m_posFlt, m_writeFreqs, m_outFile);
				
		// verify all data from single chromosome
		if(
			 !m_matCase.m_chr.equals(m_matCont.m_chr) || 
			( m_matOutg != null && !m_matCase.m_chr.equals(m_matOutg.m_chr) ) ||
		    ( m_posFlt.m_chr != null && !m_matCase.m_chr.equals(m_posFlt.m_chr) ) 
		  )	
			Utils.exit_print("[Selection]::Error::Input should be from single chromosome");
		
		// run sliding window
		tests = new Tests(m_caseSize, m_contSize, m_outgSize, m_matCase, m_matCont, m_matOutg);
		try{
			// output file
			PrintWriter outF = new PrintWriter(new BufferedWriter(new FileWriter(m_outFile)));
			outF.print("#Chr\tStart\tEnd\tcase-sites\tcase-sitesNF\tcont-sites\tcont-sitesNF\toutG-sites\toutG-sitesNF\t" +
						"caseWT\tcontWT\tcaseTT\tcontTT\tcase-cont-logTT\tcaseTD\tcontTD\tcaseSF\tcontSF\tcase-cont-logSF\tcase-contFst\tcont-outgFst\tcase-outgFst\tPBS\tcase-cont-H\n");
			
			// perform tests on window
			System.out.println("\n[Selection]::Performing tests of selection");
			int min_snp = Math.min(m_matCase.m_minPos, m_matCont.m_minPos);
			int max_snp = Math.max(m_matCase.m_maxPos, m_matCont.m_maxPos);
			
			if(max_snp - min_snp < m_windowSize){
				do_window(min_snp, outF); // one window
			}
			else{
				// sliding window
				for(int win_start = min_snp; win_start + m_windowSize <= max_snp; win_start += m_stepSize){
					do_window(win_start, outF);
				}
			}
			outF.close();
		}
		catch(IOException e){ e.printStackTrace(); System.exit(1); }
		
		finish_up();
	}	
	
	/********************* do_window *********************/
	private void do_window(int win_start, PrintWriter out){

		// get start & end columns for window (case/control and outgroup)
		int caseStartCol = m_matCase.first_col_in_window(win_start);
		int caseEndCol = m_matCase.last_col_in_window(win_start + m_windowSize);
		int contStartCol = m_matCont.first_col_in_window(win_start);
		int contEndCol = m_matCont.last_col_in_window(win_start + m_windowSize);
		int outgStartCol = -1, outgEndCol = -1;
		if(m_matOutg != null){
			outgStartCol = m_matOutg.first_col_in_window(win_start);
			outgEndCol = m_matOutg.last_col_in_window(win_start + m_windowSize);;
		}
		
		// tests
		String result = tests.testWindow(caseStartCol, caseEndCol, contStartCol, contEndCol, outgStartCol, outgEndCol, win_start, win_start + m_windowSize);				
		out.print(result);
		m_windows++;
	}

	/******************** Constructor ********************/
	public s_stats(String[] args){
		
        try {
        	m_startTime = System.currentTimeMillis(); 
        	m_windows = 0;
        	String casef, contf, outgroupf, bed, fltType;
        	
        	// parse command line
        	Options options = new Options();
            options.addOption("case", true, "case population (VCF)");
            options.addOption("cont", true, "control population (VCF)");
            options.addOption("outg", true, "out group (VCF)");
            options.addOption("o", true, "out file name");
            options.addOption("b", true, "regions-of-interest (BED), discards variants outside these regions");
            options.addOption("f", true, "filter: \"full\" (for multiple small regions, i.e. exons) or \"region\" (for fewer large regions, i.e. genes)");
            options.addOption("h", false, "show this helpful help message");
        	
        	CommandLine cmd = new BasicParser().parse(options, args);
        	
        	// help
        	if(cmd.hasOption("h"))
				Utils.exit_usage();
        	
        	// output file
    		if(cmd.hasOption("o"))
    			m_outFile = cmd.getOptionValue("o");
    		else
    			m_outFile = "tests.out";
    		
    		// regions file & filter type
    		if(cmd.hasOption("b")){
    			bed = cmd.getOptionValue("b");
    			File aoiBed = new File(bed);
    			if(!aoiBed.canRead()) Utils.exit_print("Can't read file " + bed + ". Quitting...");
    			
    			if(!cmd.hasOption("f"))
    				Utils.exit_print("filter-type must be specified when using regions file, see -h for help. Quitting...");
    			
    			fltType = cmd.getOptionValue("f");
    			m_posFlt = new PosFilter(aoiBed, fltType);	
    		}
    		else
    			m_posFlt = new PosFilter(null, ""); // non-filter
    		
        	// case VCF file & matrix
        	if(cmd.hasOption("case")){
        		casef = cmd.getOptionValue("case");
        		File caseFile = new File(casef);
    			if(!caseFile.canRead()) Utils.exit_print("Can't read file " + casef + ". Quitting...");
    			m_matCase = new SnpMat(caseFile, m_skipMultiAsnps, m_windowSize);
        		m_caseSize = 0;
        	}
        	else
        		Utils.exit_usage();
        	
        	// control VCF file & matrix
        	if(cmd.hasOption("cont")){
        		contf = cmd.getOptionValue("cont");
	        	File contFile = new File(contf);
	    		if(!contFile.canRead()) Utils.exit_print("Can't read file " + contf + ". Quitting...");
	    		m_matCont = new SnpMat(contFile, m_skipMultiAsnps,  m_windowSize);
	    		m_contSize = 0;
	    	}
        	else
        		Utils.exit_usage();
        	
        	// outgroup VCF file & matrix
    		if(cmd.hasOption("outg")){
    			outgroupf = cmd.getOptionValue("outg");
    			File outgroupFile = new File(outgroupf);
    			if(!outgroupFile.canRead()) Utils.exit_print("Can't read file " + outgroupf + ". Quitting...");
    			m_matOutg = new SnpMat(outgroupFile, m_skipMultiAsnps, m_windowSize);
    			m_outgSize = 0;
    		}
    		else
    			m_matOutg = null;
        }
        catch (ParseException e) {
			e.printStackTrace();
		}
	}

	/******************** finish_up **********************/
    private void finish_up(){
    	// get elapsed time
        long elapsedTimeMillis = System.currentTimeMillis() - m_startTime;
        float elapsedTimeMin = elapsedTimeMillis / (60 * 1000F);
    	
    	// report
    	System.out.println("[Selection]::Done!");
    	System.out.print("[Selection]::Tested " + m_windows + " windows (" + m_windowSize + " bp)");
		if(m_posFlt.m_chr != null) System.out.print(" on chr" + m_posFlt.m_chr);
		System.out.printf(", ran for %.1f minutes\n\n", elapsedTimeMin);
    }

}
