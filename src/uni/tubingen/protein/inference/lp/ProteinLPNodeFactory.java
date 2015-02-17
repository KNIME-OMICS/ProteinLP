package uni.tubingen.protein.inference.lp;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ProteinLP" Node.
 * This node is a implementation of ProteinLP algorihtm for uni.tubingen.protein inference analysis.
 *
 * @author enrique
 */
public class ProteinLPNodeFactory 
        extends NodeFactory<ProteinLPNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ProteinLPNodeModel createNodeModel() {
        return new ProteinLPNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<ProteinLPNodeModel> createNodeView(final int viewIndex,
            final ProteinLPNodeModel nodeModel) {
        return new ProteinLPNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new ProteinLPNodeDialog();
    }

}

