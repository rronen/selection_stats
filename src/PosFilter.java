import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class PosFilter {
	protected File m_bed;
	protected boolean m_useAll;
	protected HashMap<Integer, String> m_aoi;
	protected HashMap<Integer, Integer> m_regions;
	protected String m_chr;
	protected boolean m_fullFlt;
	
	/******************** CONSTRUCTOR ********************/	
	public PosFilter(File bed, String fltType){
		
		// determine if filtering needed
		m_bed = bed;
		if(m_bed == null) 
			m_useAll = true;
		else 
			m_useAll = false;
		
		// determine filter type
		if(fltType.equalsIgnoreCase("full")) 
			m_fullFlt = true;
		else if(fltType.equalsIgnoreCase("region")) 
			m_fullFlt = false;
		else if(!m_useAll) 
			Utils.exit_print("[PosFilter]::Error::Filter type should be \"full\" or \"region\"");
		
		m_chr = null;
		m_aoi = new HashMap<Integer, String>();
		m_regions = new HashMap<Integer, Integer>();
	}
	
	/******************** build_filter *******************/
	public void build_filter(){
		if(m_useAll){
			System.out.println("[PosFilter]::No positional filter");
			return;
		}
		else{
			// build filter
			System.out.print("\nBuilding positional filter...");
			int approxTotal = 0;
			try {
				 BufferedReader bedFile = new BufferedReader(new FileReader(m_bed));
		         String line;
		         // read BED file
		         while((line = bedFile.readLine()) != null){
		        	 
		        	 if(line.startsWith("#")) continue; // ignore comment lines
		        	 String[] split_line = line.split("\t");
		        	 String chr = split_line[0];
		        	 
		        	 // sanity check
		        	 if(m_chr == null) m_chr = chr;
		        	 else if(!m_chr.equals(chr)){
		        		 System.out.println("[PosFilter] Expecting one chromosome. Quitting...");
		        		 System.exit(1);
		        	 }
		        	 
		        	 // get gene name in annotated BED
		        	 String gene = "";
		        	 if(split_line.length > 3) gene = split_line[3];
		        	 
		        	 // get BED interval
		        	 try{
		        		 int start = Integer.valueOf(split_line[1]);
		        		 int end = Integer.valueOf(split_line[2]);
			        	 // update AOI table
			        	 if(m_fullFlt)
			        		 for(int i=start; i <= end; i++) m_aoi.put(i, gene);
			        	 else
			        		 m_regions.put(start, end);
			        	 
			        	 approxTotal += end - start + 1;
		        	 }
		        	 catch(NumberFormatException e){
		        		 e.printStackTrace();
		        		 System.exit(1);
		        	 }
		         }
		         bedFile.close();
			 }
			 catch (IOException e) {
				 e.printStackTrace();
		         System.exit(1);
			 }
			 
			 // write summary information
			 System.out.println("done.");
			 if(m_fullFlt)
			 	System.out.println("\t[PosFilter]::Will consider SNPs from " + m_aoi.size() + " chromosomal positions");
			 else
				 System.out.println("\t[PosFilter]::Will consider SNPs from ~" + approxTotal + " chromosomal positions");
		}
	}
	
	/**************** in_area_of_interest ****************/
	public boolean in_area_of_interest(int pos){
		if(m_useAll){
			// no filter
			return true;
		}
		else{
			if(m_fullFlt){
				// full filter
				return m_aoi.containsKey(pos);
			}
			else{
				// region filter
				for(Integer st : m_regions.keySet()){
					if(pos >= st && pos <= m_regions.get(st))
						return true;
				}
				return false;
			}
		}
	}
	
} // class
