package uni.tubingen.algorithms;


/*
 * This is the most base class for all consequent uni.tubingen.algorithms
 * 
 * It contains those commonly used functions
 * 
 */

import java.util.StringTokenizer;
import java.util.Vector;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Hashtable;

//import jxl.Workbook;
//import jxl.format.UnderlineStyle;
//import jxl.write.Label;
//import jxl.write.Number;
//import jxl.write.NumberFormat;
//import jxl.write.WritableCellFormat;
//import jxl.write.WritableFont;
//import jxl.write.WritableSheet;
//import jxl.write.WritableWorkbook;

import uni.tubingen.protein.FastaLoader;
import uni.tubingen.protein.PeptideGenerator;
import uni.tubingen.protein.Protein;
import weka.core.Utils;



public class BasicAlg {
	
	public int TotalProteins,TotalPeptides; 	// the number of candidate proteins and identified peptides in the bipartite graph
	public int TotalProteinsGroups; 			// the number of candidate uni.tubingen.protein groups in the bipartite graph after clustering the proteins with the same set of identified peptides into a group
	public double[][] coef; 					//the relationship matrix between peptides and proteins
	public double[] probability;				//the array of peptide probabilities
	public double[] coffie;						//the array of ln(1-x) transformed peptide probabilities
	public String[] proteinNames;				//the array of uni.tubingen.protein names
	public String[] peptideSeq;					//the array of peptide sequences
	public String[] proteingroupnames;			//the array of uni.tubingen.protein group names

	Hashtable distinct_peptide = new Hashtable();		// the hashtable of peptide sequences
	Hashtable peptide_prob = new Hashtable(); 			// the hashtable of peptide probabilities
	Hashtable distinct_protein = new Hashtable();		// the hashtable which saves the relationship between candidate proteins and identified peptides
	Hashtable distinct_protein_pos = new Hashtable();	// the hashtable of uni.tubingen.protein names

	Hashtable proteinListInEachGroup = new Hashtable();	// the set of proteins in each uni.tubingen.protein group
	Hashtable proteingroup = new Hashtable();			// the hashtable which saves the relationship between each uni.tubingen.protein groups and identified peptides
	Hashtable proteingroup_pos = new Hashtable();		// the hashtable of uni.tubingen.protein group names
	
	//load the information in the peptide identification file
	public void loadPeptideFile(String resultFile){
	    try{
	    	
	    		 //load the peptide sequences;
	    		 File sfile = new File(resultFile);
			     BufferedReader s_in = new BufferedReader(new FileReader(sfile));
			     String tuple;
				 StringTokenizer token = null;
				 int counter=0,seq_id=0;
				    while((tuple=s_in.readLine())!=null){
				    		token = new StringTokenizer(tuple, "\t");
				    		String peptideseq = token.nextToken();//peptide sequences
				    		Integer po = (Integer)distinct_peptide.get(peptideseq);//the hashtable stores all the peptide sequences.
					    	if(po==null){
					    		int pos = distinct_peptide.size();
					    		distinct_peptide.put(peptideseq, new Integer(pos));
					    	}
				    }
				    s_in.close();
				    TotalPeptides = distinct_peptide.size();
				   
				    
				    
				  //load the uni.tubingen.protein information
				    sfile = new File(resultFile);
				    s_in = new BufferedReader(new FileReader(sfile));
					token = null;
				    counter=0;
				    seq_id=0;
				    peptide_prob = new Hashtable(TotalPeptides);
				    while((tuple=s_in.readLine())!=null){
				    		token = new StringTokenizer(tuple, "\t");
				    		String peptideseq = token.nextToken();//peptide sequences
				    		String proteinName = token.nextToken();//uni.tubingen.protein names
				    		String peptideprob = token.nextToken();//peptide probabilities
				    		
				    		seq_id=(Integer)distinct_peptide.get(peptideseq);
				    		ArrayList prob=(ArrayList)peptide_prob.get(peptideseq);//peptide_prob hashtable stores all the probabilities corresponding to a peptide.
					    	if(prob==null){
					    		prob= new ArrayList();
					    		prob.add(peptideprob);
					    		peptide_prob.put(peptideseq,prob);
					    	}
					    	else{
					    		prob.add(peptideprob);
					    	}
					    	
					    	//deal with uni.tubingen.protein names in the peptide identification file so that these uni.tubingen.protein names can match with those in the peptide detectability file.
					    /*	if(proteinName.contains("|")){		
						    	int firstLine = proteinName.indexOf("|");
						    	int lastLine = proteinName.indexOf("|",firstLine+1);
						    	if(lastLine==-1){
						    		if(!proteinName.contains("gi|")){
						    		proteinName = proteinName.substring(0,firstLine);
						    		}
						    	}
						    	else proteinName = proteinName.substring(firstLine+1,lastLine);
					    	}*/
					    	
						    ArrayList pt = (ArrayList)distinct_protein.get(proteinName);//the hashtable stores the relationship between candidate proteins and identified peptides.
						    if(pt==null){
						    	pt = new ArrayList();
						    	pt.add(peptideseq);
						    	distinct_protein.put(proteinName, pt);
						    	int pt_pos = distinct_protein_pos.size();//the hashtable stores all the proteins.
						    	distinct_protein_pos.put(proteinName, new Integer(pt_pos));		
						    }else{
						    	if(!pt.contains(peptideseq)){
						    	pt.add(peptideseq);	
						    	}
						    }
				    }
				    
				    TotalProteins = distinct_protein.size();
				    s_in.close();
				    
				    
				    peptideSeq = new String[TotalPeptides];
				    for(Iterator h = distinct_peptide.keySet().iterator();h.hasNext(); ){
						 String key = (String) h.next();
						 int pos = ((Integer)distinct_peptide.get(key)).intValue();
						 peptideSeq[pos] = key;
				    }
					
				    //Clustering the proteins with the same set of identified peptides into a group
				    proteinNames = new String[TotalProteins];
					for(Iterator h = distinct_protein.keySet().iterator();h.hasNext(); ){
						 StringBuffer sb2 = new StringBuffer();
						 String key = (String) h.next();
						 int pos = ((Integer)distinct_protein_pos.get(key)).intValue();
						 proteinNames[pos] = key;
						 
						 ArrayList al = (ArrayList)distinct_protein.get(key);
						 for(int i=0;i<al.size();i++){
							 String peptideseq = (String)al.get(i);
							 int num = ((Integer)distinct_peptide.get(peptideseq)).intValue();
							 sb2.append(num);
						 }
						 String s2 = sb2.toString();
				
						 ArrayList pt = (ArrayList)proteinListInEachGroup.get(s2);
						 if(pt==null){
						    pt = new ArrayList();
						    pt.add(pos);
						    int pt_pos = proteinListInEachGroup.size();
						    proteinListInEachGroup.put(s2, pt);	
						    proteingroup.put(s2, al);
						    proteingroup_pos.put(s2,new Integer(pt_pos));	
						 }else{
						    if(!pt.contains(pos)){
						    	pt.add(pos);	
						    }
						 }
					}
					TotalProteinsGroups = proteingroup.size();

					
					proteingroupnames = new String[TotalProteinsGroups];
					for(Iterator h = proteingroup.keySet().iterator();h.hasNext(); ){
						 String key = (String) h.next();
						 ArrayList al = (ArrayList)proteingroup.get(key);
						 int pos = ((Integer)proteingroup_pos.get(key)).intValue();
						 proteingroupnames[pos] = key;
			
					 }
	    }catch (Exception x){
	    	x.printStackTrace(System.err);
	    }
	}
	
	//According to hashtable "peptide_prob", we initialize the array of peptide probabilities and its transformed value: ln(1-x).
	public void getcoef(String tag){	
		try{
			double peptideprob;
			probability=new double[TotalPeptides];
			coffie=new double[TotalPeptides];
			for(int j=0;j<TotalPeptides;j++){	
				ArrayList arl = (ArrayList)peptide_prob.get(peptideSeq[j]);
				if(tag.equals("average")){
					peptideprob=0;
					for(int i=0;i<arl.size();i++){
						peptideprob=peptideprob+Double.parseDouble((String)arl.get(i));
					}
					probability[j]=peptideprob/arl.size();
				}
				else if(tag.equals("max"))
				{
					peptideprob=Double.parseDouble((String)arl.get(0));
					for(int i=0;i<arl.size();i++){
						if(peptideprob<Double.parseDouble((String)arl.get(i))){
							peptideprob=Double.parseDouble((String)arl.get(i));
						}
					}
					probability[j]=peptideprob;
				}
				else
					probability[j]=1;
				
				if(probability[j]==1.0){
					probability[j]=0.99999;
					double t = 1-probability[j];
					coffie[j]=Math.log(t);
				}
				else 
				{
					double t = 1-probability[j];
					coffie[j]=Math.log(t);
				}

			}
	    }
		catch (Exception x){
		    x.printStackTrace(System.err);
		}
	}
}
