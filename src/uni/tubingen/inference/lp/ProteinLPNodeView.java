package uni.tubingen.inference.lp;

import org.knime.core.node.NodeView;



/**
 * <code>NodeView</code> for the "ProteinLP" Node.
 * This node is a implementation of ProteinLP algorihtm for uni.tubingen.protein inference analysis.
 *
 * @author enrique
 */
public class ProteinLPNodeView extends NodeView<ProteinLPNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link ProteinLPNodeModel})
     */
    protected ProteinLPNodeView(final ProteinLPNodeModel nodeModel) {
        super(nodeModel);
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    	ProteinLPNodeModel nodeModel = 
                getNodeModel();
            assert nodeModel != null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
        // TODO: generated method stub
    }

}

