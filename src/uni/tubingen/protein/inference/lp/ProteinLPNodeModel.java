package uni.tubingen.protein.inference.lp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.lang3.StringUtils;
import org.gnu.glpk.GLPK;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowIterator;
import org.knime.core.data.StringValue;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import uni.tubingen.algorithms.PILP;


/**
 * This is the model implementation of ProteinLP.
 * This node is a implementation of ProteinLP algorihtm for protein inference analysis.
 *
 *  A linear programming model for protein inference problem in shotgun proteomics. Bioinformatics, 2012
 *  Ting Huang and Zengyou He
 *  
 *  
 * @author enrique
 */
public class ProteinLPNodeModel extends NodeModel {
	
	
private static final NodeLogger logger = NodeLogger.getLogger("ProteinLP probabilities");
    
    static String CFGKEY_PEPTIDES = "peptides";
	static String CFGKEY_PROTEIN = "protein";
	static String CFGKEY_PROBABILITIES = "probabilities";
	//static String CFGKEY_LAMMDA_PARAMETER = "lammda";
	
	//fields to link execute variable with input variable...
	private final SettingsModelString m_peptide_column = new SettingsModelString(CFGKEY_PEPTIDES, "Peptides");
	private final SettingsModelString m_protein_column   = new SettingsModelString(CFGKEY_PROTEIN, "Protein");
	private final SettingsModelString m_probability_column   = new SettingsModelString(CFGKEY_PROBABILITIES, "Probabilities");
	//private final SettingsModelDoubleBounded lammda_parameter = new SettingsModelDoubleBounded(CFGKEY_LAMMDA_PARAMETER, 0.9, 0.0, 1.0);
	
	//fields to manage the input table (column index)...
		static int pep_idx    = 0;
		static int accsn_idx  = 0;
		static int proba_idx  = 0;
		
		static File tmp_identification_file = null;
		
		static BufferedDataContainer container = null;
		
    
    /**
     * Constructor for the node model.
     */
    protected ProteinLPNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
          //loading native library (version-specific) for solve optimization problems...
            try {
            	System.load("/media/Datos/TRABAJO/KNIME workflow/workspace/ProteinInferenceLP/Native-Library/libglpk_java.so");
            } catch (UnsatisfiedLinkError e) {
              System.err.println("Native code library failed to load.\n" + e);
              System.exit(1);
            }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	  this.checkTableConfiguraion(inData);
          
          DataTableSpec new_spec_table = new DataTableSpec(make_output_spec());  	
          container = exec.createDataContainer(new_spec_table);
          
          System.out.println(GLPK.glp_version( )) ;
  		  
          buildTemporalIdentificationFile (inData[0]); //building temporal_identification_file FILE...

  		double para=0.0;
  		double startTime = System.currentTimeMillis();
  		
  		PILP alg = new PILP(); 
  		
  		if (tmp_identification_file.exists()){
  			alg.initialization(tmp_identification_file.getAbsolutePath());
  		}else{
  			throw new IOException("File identification not create...");
  		}
  		
  		alg.ProteinLP(para);
  		System.out.println("------------------------------------------");
  		System.out.println("The value of parameter sigma="+para);
  		 
  		alg.print(container);
  		
  		double endTime = System.currentTimeMillis();
  		double running_time = (endTime-startTime)/(double)1000;
  		System.out.println("Running Time:"+running_time);
          
         tmp_identification_file.delete();
         
         return new BufferedDataTable[]{ container.getTable() };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        // TODO: generated method stub
    	return new DataTableSpec[]{new DataTableSpec(this.make_output_spec())};
    }
    
    
  //method for checking the table configuration coming...
    private void checkTableConfiguraion(BufferedDataTable[] inData) throws Exception{
    	
    	pep_idx  = inData[0].getDataTableSpec().findColumnIndex(m_peptide_column.getStringValue());
    	accsn_idx= inData[0].getDataTableSpec().findColumnIndex(m_protein_column.getStringValue());
    	proba_idx= inData[0].getDataTableSpec().findColumnIndex(m_probability_column.getStringValue());
    	
    	if (pep_idx < 0 || accsn_idx < 0 || proba_idx < 0  || pep_idx == accsn_idx ) {
    		throw new Exception("Illegal columns: "+m_peptide_column+" "+m_protein_column+" "+m_probability_column+" , re-configure the node!");
    	}
    }
    
 
     private DataColumnSpec[]  make_output_spec() {
    	
    	DataColumnSpec cols[] = new DataColumnSpec[2];
    	cols[0] = new DataColumnSpecCreator("Protein ID", StringCell.TYPE).createSpec();
    	cols[1] = new DataColumnSpecCreator("ProteinLP Probability", StringCell.TYPE).createSpec();
	
      return cols;
	}
    
     
     //building Identification file. Mandatory input for ProteinLP algorithm
     private static void buildTemporalIdentificationFile (BufferedDataTable inData) throws IOException{
    	  
    	  tmp_identification_file = File.createTempFile("identification_file", ".txt");
      	  PrintWriter pw = new PrintWriter(new FileWriter(tmp_identification_file));
  		
    	   try {			
	  
    	    RowIterator row_it = inData.iterator();
	    	while (row_it.hasNext()) {
	    		
	     		//getting information from current row...
	    		DataRow r = row_it.next();
	    		DataCell pep_cell        =  r.getCell(ProteinLPNodeModel.pep_idx);
	    		DataCell accsn_cell     =  r.getCell(ProteinLPNodeModel.accsn_idx);
	    		DataCell proba_cell    =  r.getCell(ProteinLPNodeModel.proba_idx);
	
	    		//rows with missing cells cannot be processed (no missing values in PSM graph...)
	    		if (pep_cell.isMissing() || accsn_cell.isMissing() || proba_cell.isMissing()) {
	    			continue;
	    		}
	    		
	    		//getting value from cells
	    		String peptide_entry   =   ((StringValue) pep_cell).getStringValue();
	    		String protein_accsn   =   ((StringValue) accsn_cell).getStringValue();
	    		String proba_entry    =   ((StringValue) proba_cell).getStringValue();
	    		 
	    		    String [] protein_group = protein_accsn.split(";"); //compatible with idXML file in OpenMS workflow...
	    		 
	    	    	for(int i = 0; i < protein_group.length; i++){				
	                     pw.println(peptide_entry + "\t" + formattedProteinAccesion (protein_group[i]) + "\t" + proba_entry);
	    	        }			 
				}	
	    	
	    	}
	    	catch (Exception x){
	    		x.printStackTrace(System.err);
		    }
    	   pw.close();
      }
      
      
    //this function is for formatting long protein name (getting universal ID...)
  	static private String formattedProteinAccesion (String protein_accn){
  		      if(protein_accn.contains("|")){
  		    	  return StringUtils.substringBeforeLast(protein_accn, "|");
  		      }
  		      else if (protein_accn.contains(" ")){
  		    	  return protein_accn.trim(); //avoiding protein identifier with space...
  		      }	  
  		      else{
  		    	  return protein_accn;
  		      }
  	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_peptide_column.saveSettingsTo(settings);
        m_protein_column.saveSettingsTo(settings);
        m_probability_column.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_peptide_column.loadSettingsFrom(settings);
        m_protein_column.loadSettingsFrom(settings);
        m_probability_column.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_peptide_column.validateSettings(settings);
        m_protein_column.validateSettings(settings);
        m_probability_column.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}

