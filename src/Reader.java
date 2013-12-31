import java.util.ArrayList;
import java.util.HashMap;

/**
 * General population-sample reader class.
 * Contains the logic and functionality that is shared across VCF & MPOP readers. 
 */
public class Reader {
	
	/****************** class variables ******************/
	protected ArrayList<ArrayList<Integer>> mat;
	protected HashMap<Integer, String> row2sampNames;
	protected HashMap<Integer, Integer> col2pos;
	protected HashMap<Integer, Double> col2pos_d;
	protected int min_pos, max_pos, currPos, prevPos;
	protected int numSites, numHaps;
	protected String prevChrom;
	
	/*********************** C'tor ***********************/
	public Reader(){
		
		prevChrom = "-1";
		numHaps = 0; numSites = 0;
		currPos = 0; prevPos = 0;
		max_pos = -1; min_pos = -1;

		mat = new ArrayList<ArrayList<Integer>>();
		row2sampNames = new HashMap<Integer, String>();
		col2pos = new HashMap<Integer, Integer>();
		col2pos_d = new HashMap<Integer, Double>();
	}
	
	/********************** read **********************/
	public void read(){ /* overridden in all extending classes */ }
	
	/********************** get_mat **********************/
	public ArrayList<ArrayList<Integer>> get_mat(){ return mat; }
	
	/***************** get_row2sampNames *****************/
	public HashMap<Integer, String> get_row2sampNames(){ return row2sampNames; }
	
	/******************** get_col2pos ********************/
	public HashMap<Integer, Integer> get_col2pos() { return col2pos; }
	
	/******************* get_col2pos_d *******************/
	public HashMap<Integer, Double> get_col2pos_d() { return col2pos_d; }
	
	/******************** get_min_pos ********************/
	public int get_min_pos(){ return min_pos; }
	
	/******************** get_max_pos ********************/
	public int get_max_pos(){ return max_pos; }

	/******************** get_num_haps *******************/
	public int get_num_haps(){ return numHaps; }
	
	/******************* get_num_sites ******************/
	public int get_num_sites(){ return numSites; }

	/********************** get_chr *********************/
	public String get_chr(){ return prevChrom; }	
	
}
