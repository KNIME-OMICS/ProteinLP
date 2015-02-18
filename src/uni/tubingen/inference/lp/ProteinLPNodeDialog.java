package uni.tubingen.inference.lp;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DoubleValue;
import org.knime.core.data.StringValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.ColumnFilter;


/**
 * <code>NodeDialog</code> for the "ProteinLP" Node.
 * This node is a implementation of ProteinLP algorihtm for protein inference analysis.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author enrique
 */
public class ProteinLPNodeDialog extends DefaultNodeSettingsPane {
	/**
	 * New pane for configuring the ProteinLP node.
	 */
	protected ProteinLPNodeDialog() {
		super ();
		
		//fields to match with coming table...
		final SettingsModelString matches_peptides = new SettingsModelString(ProteinLPNodeModel.CFGKEY_PEPTIDES, "Peptides");
		final SettingsModelString accsn_protein    = new SettingsModelString(ProteinLPNodeModel.CFGKEY_PROTEIN, "Protein");
		final SettingsModelString probabilities    = new SettingsModelString(ProteinLPNodeModel.CFGKEY_PROBABILITIES, "Probabilities");
		
		addDialogComponent(new DialogComponentColumnNameSelection(accsn_protein, "Proteins Column", 0, true, StringValue.class));
		
		addDialogComponent(new DialogComponentColumnNameSelection(matches_peptides, "Peptides Column", 0, true, 
				
				new ColumnFilter() {
			
			@Override
			public boolean includeColumn(DataColumnSpec colSpec) {
				if (colSpec.getType().isCollectionType() && colSpec.getType().getCollectionElementType().isCompatible(StringValue.class))
					return true;
				
				if (colSpec.getType().isCompatible(StringValue.class)) 
					return true;
				
				return false;
			}
			
			@Override
			public String allFilteredMsg() {
				return "No suitable columns (string or List/Set column) to select!";
			}
		}));
		
		addDialogComponent(new DialogComponentColumnNameSelection(probabilities, "Probabilities", 0, true, DoubleValue.class));
	}
}

