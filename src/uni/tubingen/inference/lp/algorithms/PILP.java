package uni.tubingen.inference.lp.algorithms;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;

import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.GlpkException;
import org.gnu.glpk.SWIGTYPE_p_double;
import org.gnu.glpk.SWIGTYPE_p_int;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataRow;
import org.knime.core.data.RowKey;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;

import weka.core.Utils;


public class PILP extends BasicAlg{
	
	/** the array of uni.tubingen.protein group probabilities*/
	private double[] proteingroupProb;
	
	/** the number of non-zero variables in the matrix P */
	private int nonZeroEntriesM2 = 0;
	
	/** the set of peptides corresponds to each uni.tubingen.protein group */
	private ArrayList<String>[] peptideListInEachProteinGroup;
	
	/** the set of peptide IDs in each uni.tubingen.protein group */
	private Hashtable<Integer, Integer>[] candidateList;
	
	
	/* Load the peptide identification data and make it suitable for proteinLP processing*/
	public void initialization(String resultFile){
		
		loadPeptideFile(resultFile);
		getcoef("max");
		
		peptideListInEachProteinGroup = new ArrayList[TotalProteinsGroups];
		for(Iterator h = proteingroup.keySet().iterator(); h.hasNext(); ){
			 String key = (String) h.next();
			 ArrayList<String> al = (ArrayList<String>)proteingroup.get(key);
			 int pos = ((Integer)proteingroup_pos.get(key)).intValue();
			 peptideListInEachProteinGroup[pos] = new ArrayList<String>();
			 for(int i=0;i<al.size();i++){
				 peptideListInEachProteinGroup[pos].add(al.get(i));
			 }
		}
		
		candidateList = new Hashtable[TotalProteinsGroups];
		for(int i=0;i<TotalProteinsGroups;i++){
			candidateList[i] = new Hashtable<Integer, Integer>();
			for(int j=0;j<peptideListInEachProteinGroup[i].size();j++){
				String peptideseq = peptideListInEachProteinGroup[i].get(j);
				int num = ((Integer)distinct_peptide.get(peptideseq));
				candidateList[i].put(num, num);
			}
		}
		
		//count the number of non-zero variables in the matrix P
		int counter=0;
		for(int i=0;i<TotalProteinsGroups;i++){
	        for(int j=0;j<TotalPeptides;j++){
	        	if(candidateList[i].get(new Integer(j))!=null){
	        		counter++;
	        	}
	        }
		}
		nonZeroEntriesM2=counter;
	}
	
	
	public void ProteinLP(double t){
		System.out.println("TotalPeptides:"+TotalPeptides);
		System.out.println("TotalProteins:"+TotalProteins);
		
		glp_prob solver;
		glp_smcp parm;
		SWIGTYPE_p_int row, col;
		SWIGTYPE_p_double coef;
		int ret;
		
		try {	
		// Create problem
		solver = GLPK.glp_create_prob();
		System.out.println("Problem created");
		GLPK.glp_set_prob_name(solver, "LP_LP");
		GLPK.glp_set_obj_name(solver, "z");
		GLPK.glp_set_obj_dir (solver, GLPKConstants.GLP_MAX); 
     
		        
		/* Specify the number of constraints and their bounds */
		//GLPK.glp_add_rows(solver, TotalPeptides*2+TotalProteins*2);
		/* Specify the number of constraints and their bounds */
		//GLPK.glp_add_rows(solver, TotalPeptides*2+TotalProteins*2);
		GLPK.glp_add_rows(solver, TotalPeptides*2+nonZeroEntriesM2);
		for(int i=0;i<TotalPeptides;i++){
			double para=t;
			double temp=0.0;
			temp = Math.log(1-para-probability[i]);
			
			GLPK.glp_set_row_name(solver, i+1, "c"+(i+1));
			GLPK.glp_set_row_bnds(solver, i+1, GLPKConstants.GLP_LO,temp,TotalProteinsGroups);
		}
		for(int i=0;i<TotalPeptides;i++){
			double para=t;
			double temp=0.0;
			temp = Math.log(1+para-probability[i]);
			
			GLPK.glp_set_row_name(solver, TotalPeptides+i+1, "c"+(TotalPeptides+i+1));
			GLPK.glp_set_row_bnds(solver, TotalPeptides+i+1, GLPKConstants.GLP_UP,-TotalProteinsGroups,temp);
		}
		int num=0;
		for(int j=0;j<TotalProteinsGroups;j++){
			for(int i=0;i<TotalPeptides;i++){
				if(candidateList[j].get(new Integer(i))!=null){
		        		GLPK.glp_set_row_name(solver, TotalPeptides*2+num+1, "c"+(TotalPeptides*2+num+1));
		    			GLPK.glp_set_row_bnds(solver, TotalPeptides*2+num+1, GLPKConstants.GLP_LO,0.0d,Double.POSITIVE_INFINITY);
		    			num++;
		        }
		 	}		
		}
		num=0;
		/* Specify the number of variables and precise their types, bounds
	 	* and coefficients in the objective function */
		GLPK.glp_add_cols(solver, TotalProteinsGroups+nonZeroEntriesM2);
		for(int j=0;j<TotalProteinsGroups;j++){
			for(int i=0;i<TotalPeptides;i++){
				if(candidateList[j].get(new Integer(i))!=null){
		        		num++;
		        		GLPK.glp_set_col_name(solver,num, "x"+(num));
		    			GLPK.glp_set_col_kind(solver,num, GLPKConstants.GLP_CV);
		    		    GLPK.glp_set_col_bnds(solver,num, GLPKConstants.GLP_DB, Double.NEGATIVE_INFINITY, 0.0d);
		    			GLPK.glp_set_obj_coef(solver,num, 0);
		        }
		 	}		
		}
		for(int j=0;j<TotalProteinsGroups;j++){
			   GLPK.glp_set_col_name(solver,nonZeroEntriesM2+j+1, "x"+(nonZeroEntriesM2+j+1));
			   GLPK.glp_set_col_kind(solver,nonZeroEntriesM2+j+1, GLPKConstants.GLP_CV);
			   GLPK.glp_set_col_bnds(solver,nonZeroEntriesM2+j+1, GLPKConstants.GLP_DB, Double.NEGATIVE_INFINITY, 0.0d);
			   GLPK.glp_set_obj_coef(solver,nonZeroEntriesM2+j+1, 1);
		}
		
		        
		/*Specify the constraints in a matrix */
		int nonZeroEntries = nonZeroEntriesM2*4;
		System.out.println("nonZeroEntries:"+nonZeroEntries);
		row = GLPK.new_intArray(nonZeroEntries+1);
		col = GLPK.new_intArray(nonZeroEntries+1);
		coef = GLPK.new_doubleArray(nonZeroEntries+1);
		        
		int counter=1;
		num=0;
		for(int j=0;j<TotalProteinsGroups;j++){
			for(int i=0;i<TotalPeptides;i++){
				if(candidateList[j].get(new Integer(i))!=null){
					GLPK.intArray_setitem(row, counter, i+1);
					GLPK.intArray_setitem(col, counter, num+1);
					GLPK.doubleArray_setitem(coef, counter,1);
		        	counter++;
		        	num++;
		        }
		    }
		}
		
		num=0;
		for(int j=0;j<TotalProteinsGroups;j++){
			for(int i=0;i<TotalPeptides;i++){
				if(candidateList[j].get(new Integer(i))!=null){
					GLPK.intArray_setitem(row, counter, TotalPeptides+i+1);
					GLPK.intArray_setitem(col, counter, num+1);
					GLPK.doubleArray_setitem(coef, counter,1);
		        	counter++;
		        	num++;
		        }
		    }
		}
		
		num=0;
		for(int j=0;j<TotalProteinsGroups;j++){
			for(int i=0;i<TotalPeptides;i++){
				if(candidateList[j].get(new Integer(i))!=null){
					GLPK.intArray_setitem(row, counter, TotalPeptides*2+num+1);
					GLPK.intArray_setitem(col, counter, num+1);
					GLPK.doubleArray_setitem(coef, counter, 1);
		        	counter++;
		        	GLPK.intArray_setitem(row, counter, TotalPeptides*2+num+1);
					GLPK.intArray_setitem(col, counter, nonZeroEntriesM2+j+1);
					GLPK.doubleArray_setitem(coef, counter,-1);
		        	counter++;
		        	num++;
		        }
		 	}		
		}

		GLPK.glp_load_matrix(solver,nonZeroEntries, row, col, coef);    
		parm = new glp_smcp();
		GLPK.glp_init_smcp(parm);
		ret = GLPK.glp_simplex(solver, parm);
		
		if (ret == 0) {
			write_lp_solution(solver);
		} else {
			System.out.println("The problem could not be solved");
		}
		// Free memory
		GLPK.glp_delete_prob(solver);       
		} catch (GlpkException ex) {
			ex.printStackTrace();
		}
	}
	
	void write_lp_solution(glp_prob lp) {
		int i,j;
		int n;
		int num=0;
		String name;
		double val;
		name = GLPK.glp_get_obj_name(lp);
		val = GLPK.glp_get_obj_val(lp);
		System.out.print(name);
		System.out.print(" = ");
		System.out.println(val);
		n = GLPK.glp_get_num_cols(lp);
		proteingroupProb = new double[TotalProteinsGroups];
		double[][] weigh = new double[TotalPeptides][TotalProteinsGroups];
		for (j = 0; j <TotalProteinsGroups; j++) {
			for (i = 0; i <TotalPeptides; i++) {
				weigh[i][j]=0.0;
				if(candidateList[j].get(new Integer(i))!=null){
					name = GLPK.glp_get_col_name(lp, num+1);
					val = GLPK.glp_get_col_prim(lp, num+1);
					weigh[i][j]=val;
					//System.out.print(name);
					//System.out.print(" = ");
					//System.out.println(val);
					num++;
				}
			}
		}
		for (j = 0; j <TotalProteinsGroups; j++) {
			proteingroupProb[j]=1.0;
			for (i = 0; i <TotalPeptides; i++) {
				if(candidateList[j].get(new Integer(i))!=null){
					proteingroupProb[j]=proteingroupProb[j]*(Math.pow(Math.E,weigh[i][j]));
				}
			}
			proteingroupProb[j]=1-proteingroupProb[j];
		}
	}
	
	
	public void print(BufferedDataContainer container){
		try{
			int[] pos = Utils.sort(proteingroupProb);
	    	
	    	for(int i=0;i<TotalProteinsGroups;i++){
				int bestIndex =TotalProteinsGroups-i-1;
				ArrayList pt = (ArrayList)proteinListInEachGroup.get(proteingroupnames[pos[bestIndex]]);
				
				StringBuilder sb = new StringBuilder();
	    		for(int j=0; j < pt.size(); j++) {
	    			int num = ((Integer)pt.get(j)).intValue();
	    					
	    			if (sb.length() > 0) {
	    				sb.append(';');
	    			}
	    			sb.append(proteinNames[num]);
	    		}
	    		
    			RowKey key = new RowKey(sb.toString());
                DataCell[] cells = new DataCell[4];	    		
        		cells[0] = new StringCell(sb.toString());
        		cells[1] = new DoubleCell(proteingroupProb[pos[bestIndex]]);
        		
        		Set<String> sequenceWithoutMods = new HashSet<String>();
        		for (String modSeq : peptideListInEachProteinGroup[pos[bestIndex]]) {
        			sequenceWithoutMods.add(modSeq.replaceAll("\\([^\\)]+\\)", ""));
        		}
        		
        		cells[2] = new IntCell(peptideListInEachProteinGroup[pos[bestIndex]].size());
        		cells[3] = new IntCell(sequenceWithoutMods.size());
        		
        		DataRow row = new DefaultRow(key, cells);
        		container.addRowToTable(row);
			}	

	    }
		catch (Exception x){
			x.printStackTrace(System.err);
	    }
	    finally{
	    	if (null != container)
				container.close();
	    }
	}
		
}