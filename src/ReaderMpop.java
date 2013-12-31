import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReaderMpop extends Reader {
	
	/****************** class variables ******************/
	private File mpop;
	private int win_size;
	
	/*********************** C'tor ***********************/
	public ReaderMpop(File mpop_file, int w_size){
		super();
		mpop = mpop_file;
		win_size = w_size;
	}
	
	/*********************** read ************************/
	public void read(){
		String line;
		boolean have_pos = false;
		try {
			 BufferedReader mpopFile = new BufferedReader(new FileReader(mpop));
	         
			 MPOP:
			 while((line = mpopFile.readLine()) != null){

				 // get info
	        	 if(!have_pos){
	        		 have_pos = get_info(line);
	        		 continue MPOP;
	        	 }
	        	 
	        	 // read haplotype
	        	 mat.add(new ArrayList<Integer>());
	        	 
	        	 for(int i = 0; i < line.length(); i++){
	        		 int gen = Integer.parseInt( Character.toString(line.charAt(i)) );
	        		 mat.get(numHaps).add(gen);
	        	 }
	        	 numHaps++;
	         }
			 
	         mpopFile.close();
		}
		catch (IOException e) {
			e.printStackTrace();
			Utils.exit_print("[MpopReader]::problem reading VCF file.");
		}
		
		// update haplotype count
		min_pos = col2pos.get(0);
		max_pos = col2pos.get(numSites-1);
		
		// report
		System.out.println("\n[ReaderMpop]::Read data for " + numHaps + " haplotypes, over " + numSites + " sites (min-pos: " + min_pos + " max-pos: " + max_pos +")");
	}

	/********************* get_info **********************/
	private boolean get_info(String line){
		
		if(line.startsWith("positions")){
			// get SNP positions
			String[] sp = line.split("\\s+");
			for(int i=1; i < sp.length; i++){
				double pos = Double.parseDouble(sp[i]);
				col2pos.put(i-1, (int) Math.floor(pos * win_size) ); // some keys (SNP index) may end up having identical values (SNP position)
																	 // this is since mpop & ms assume infinite sites, and use [0,1] positions
				
				col2pos_d.put(i-1, pos); // actual [0,1] position
			}
			return true;
		}
		else if(line.startsWith("derived_sites")){
			// get number of derived sites
			String[] sp = line.split("\\s+");
			numSites = Integer.parseInt(sp[1]);
		}
		return false;
	}

}
