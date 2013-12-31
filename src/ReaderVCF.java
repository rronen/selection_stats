import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReaderVCF extends Reader {
	
	/****************** class variables ******************/
	private File vcf;
	private PosFilter posFlt;
	private boolean skip_MA_snps;
	protected int col_1st_gen;	
	protected int skippedMA, skippedAOI; // statistics
	
	/*********************** C'tor ***********************/
	public ReaderVCF(File VCF, boolean skip_MA, PosFilter pFlt){
		super();
		posFlt = pFlt; 
		vcf = VCF;
		skip_MA_snps = skip_MA;
		col_1st_gen = 0; 
		skippedMA = 0; skippedAOI = 0;
	}
	
	/*********************** read ***********************/
	public void read(){
		System.out.print("\n[VCFReader]::Reading "+ vcf.getName() + " into SNP matrix, ");
		String line;
		int seenSnps = 0;
		try {
			 BufferedReader vcfFile = new BufferedReader(new FileReader(vcf));
	         
	         // info lines
	         while((line = vcfFile.readLine()) != null)
	        	 if(line.startsWith("##")) 
	        		 continue;
	        	 else
	        		 break;
	         
	         // get sample info from VCF header
	         sample_info_from_header(vcfFile, line);	         
	         System.out.println("total " + numHaps / 2 + " samples");
	         
	         // read SNP data	         
	         VCF:
	         while ((line = vcfFile.readLine()) != null) {
	        	 seenSnps++;
	        	 String[] vcfLine = line.split("\t");

	        	 // verify VCF sorted
	        	 verify_sort(vcfLine[0], vcfLine[1]);
	        	 
	        	// check SNP in area-of-interest
	        	 if(posFlt.in_area_of_interest(currPos)){
	        		 
	        		 String[] altGenoypes = vcfLine[4].split(",");
	        		 if(skip_MA_snps && altGenoypes.length > 1){
		        		 // verify not multi-allelic
	        			 skippedMA++;
	        			 continue VCF;
	        		 }
	        		 
	        		 /* possibly filter on [5]Quality [6]Filter [7]INFO [8]FORMAT */
	        		 
	        		 // read & store genotypes
	        		 for (int i = col_1st_gen; i < vcfLine.length; i++){
		        		 int hapRow1 = (i - col_1st_gen) * 2;
		        		 int hapRow2 = (i - col_1st_gen) * 2 + 1;
	        			 String[] gen_field = vcfLine[i].split(":");
	        			 String[] gen = gen_field[0].split("/|\\|"); // 1 or 2 alleles from {0,1,.}, possibly from {0,1,2,3,.}

	        			 int geno0 = -1, geno1 = -1;
	        			 
	        			 if(gen.length == 1){
	        				 // haploid genotype (chrY/X for males), count as homozygous (for frequency calc.)  
	        				 if(gen[0].equals(".")){
		        				 // missing genotype, leave -1's 
	        				 }
	        				 else{
	        					 geno0 = Integer.parseInt(gen[0]);
		        				 geno1 = Integer.parseInt(gen[0]);	 
	        				 }
	        			 }
	        			 else if(gen.length == 2){
	        				 // diploid genotype
	        				 if(gen[0].equals(".") && gen[1].equals(".")){
		        				 // missing genotype, leave -1's
		        			 }
		        			 else if(gen[0].equals(".") || gen[1].equals(".")){
		        				 Utils.exit_print("[VCFReader]::Unexpected genotype field: " + gen_field[0]);
		        			 }
		        			 else{
		        				 geno0 = Integer.parseInt(gen[0]);
		        				 geno1 = Integer.parseInt(gen[1]);		        				 
		        			 }
	        			 }
	        			 else{
	        				 // sanity check
	        				 Utils.exit_print("\n[VCFReader]::Genotype field error: " + gen_field[0]);
	        			 }

	        			 mat.get(hapRow1).add(numSites, geno0);
	        			 mat.get(hapRow2).add(numSites, geno1);
		        	 }
	        		 
	        		 col2pos.put(numSites, currPos);
	        		 if(numSites == 0)
	        			 min_pos = currPos;
	        		 
	        		 numSites++;
	        	 }
	        	 else
	        		 skippedAOI++; // not in area of interest, skip
	        
	         } // end while (VCF_LINE)
	         
	         max_pos = currPos;
	         System.out.println("[VCFReader]::Read " + seenSnps + " SNPs, kept " + numSites + " (bi-allelic SNPs), min-pos " + min_pos + ", max-pos " + max_pos);
		 }
		catch (IOException e) {
			e.printStackTrace();
			Utils.exit_print("[VCFReader]::problem reading VCF file.");
		}

	}
	
	/*************** sample_info_from_header **************/
	private void sample_info_from_header(BufferedReader vcfFile, String line) throws IOException {
		// verify header & get samples-info
		if (!line.startsWith("#")) {
			vcfFile.close();
			Utils.exit_print("[SnpMat] header missing from VCF file. Quitting...");
		}
		 
		String[] header = line.split("\t");
		col_1st_gen = get_first_geno_col(header);
		if(col_1st_gen == 0){
			vcfFile.close();
			Utils.exit_print("[SnpMat] Error in VCF header, FORMAT missing. Quitting...");
		}
		if(header.length <= col_1st_gen){
			vcfFile.close();
			Utils.exit_print("[SnpMat] VCF header indicates 0 samples. Quitting...");
		}
		 
		for (int i = col_1st_gen; i < header.length; i++){
			// remember name & create list for hap1
			row2sampNames.put(numHaps, header[i]);
			mat.add(numHaps, new ArrayList<Integer>());
			numHaps++;
			
			// remember name & create list for hap2	        	 
			row2sampNames.put(numHaps, header[i]);
			mat.add(numHaps, new ArrayList<Integer>());
			numHaps++;
		}
	}
	
	/***************** get_first_geno_col ****************/
	private int get_first_geno_col(String[] headerLine) {
		for (int i=0; i < headerLine.length; i++)
			if(headerLine[i].equals("FORMAT")) 
				return i+1;
		return 0;
	}

	/********************* check_sort ********************/
	private void verify_sort(String chr, String pos){
	   	 /* verifies VCF is position-sorted */
		
		if(prevChrom.equals("-1")) prevChrom = chr; // first data line
		
		// verify same chromosome		
	   	 if(!prevChrom.equals(chr))
	   		Utils.exit_print("\nError, expecting data from one chromosome only!");
	   	 
	   	 // read chromosome and position
	   	 prevChrom = chr;
	   	 try{
	   		currPos = Integer.parseInt(pos);
	   	 }
	   	 catch(NumberFormatException e){
	   		e.printStackTrace(); 
	   	 }
	   	 
	   	 // !!! can be STRICTLY less than, to handle double called SNPs !!!
	   	 
	   	 if(currPos < prevPos)
	   		Utils.exit_print("\n[ReaderVCF]::Error: input VCF unsorted! (expecting sorted within chromosome)");
	   	 
	   	 prevPos = currPos;
	}
		
}
