Knime-ProteinLP
===============

# About KNIME-ProteinLPInference Node

ProteinLPInference is a KNIME node implementation to perform protein inference analysis on Mass Spectrometry data. The ProteinLP model introduces the joint probability that both a protein and its constituent peptide are present in the sample. It is based on a linear model, where a mathematical transformation of this joint probability is used. The marginal probability of a peptide being present is expressed as a formula in terms of the linear combination of these variables. The model assumes that the marginal probability of each identified peptide being present is known, then the protein inference problem is formulated as the following optimization problem: ‘minimize the number of proteins of non-zero probabilities while the calculated peptide probability should be as close to its known value as possible’.

# License

knime-protein_lp is a JOmics licensed under [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0.txt).

# How to cite it:

Huang T, He Z. A linear programming model for protein inference problem in shotgun proteomics. Bioinformatics. 2012;28:2956-62.

(http://www.ncbi.nlm.nih.gov/pubmed/22954624)

# Main Features

* Fast and efficient processing of large dataset (tabular format like input).

* Simple structure (one port input/one port output).

* Requires peptide's Posterior Error Probabilities (PEPs)like input parameters.

* Compute posterior error probabilities for proteins.

* Possibility of easy integration with OpenMS workflow.


# Node Input/Output

Basically, the ProteinLPInference node handle the data in tabular format. Input: Peptide list with proteins group and peptide probabilities. Output: report Protein list with probability associated.


**Note**: The node is still evolving, we are committed to expand the node and add more features such as additional settings in the configuration dialog to controlling the optimization process. Example: possibility of setup the epsilon parameter.

# Getting ProteinLPInference node

## Installation Requirements

* Java: KNIME SDK 9.2 (or above) and Java JRE 1.6 (or above), which you can download for free here. (Note: most computers should have Java installed already).

* Operating System: The current version has been tested on Linux and Max OS X, it may require additional adjustment for other platform. If you come across any problems on your platform, please contact us (enriquea@cim.sld.cu).

* Memory: MS dataset can be very large sometimes, in order to get good performance from this KNIME node, we recommend the following settings for VM arguments: -ea -Xmx1G -XX:MaxPermSize=512M. For additional information see https://tech.knime.org/test-your-node.

## Launch via Eclipse project

Click [here]
(https://github.com/KNIME-OMICS/ProteinLP) to launch directly the latest ProteinLPInference (KNIME node) implementation.

*This node is ready to use. However, taking into account you current platform, it may requiere additional adjustment. You can build the standalone node in few step (see https://tech.knime.org/developer/documentation/export).*

## Download

You can get the latest ProteinLPInference node implementation from our Download Section. Unzipping the file and load via Eclipse project.

##Maven Dependency


# Getting Help

If you have questions or need additional help, please contact us via Enrique Audain Martinez (enriquea@cim.sld.cu) and Yasset Perez-Riverol (yperez@ebi.ac.uk).

Please send us your feedback, including error reports, improvement suggestions, new feature requests and any other things you might want to suggest.

# Screenshots

Peptide/Protein Identification workflow (OpenMS-KNIME). (1) Mass spectrometry data pre-processing. (2) Database searching/Peptide identification. (3) idXML file processing. (4) Peptide detectability prediction based on SVM approach. (5) Protein inference analysis.

 




