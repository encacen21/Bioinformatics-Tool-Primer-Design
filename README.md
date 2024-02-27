# TFM
## Genome Marker Analysis Tool

### Overview

This Python script is designed for the analysis and optimization of DNA sequences for primer design, consensus sequence generation, and genome annotation visualization. It allows users to process multiple DNA sequences, identify optimal primers, generate consensus sequences, and visualize annotations on genomic sequences. This tool integrates functionalities for reading and processing genome records, optimizing primer sequences based on melting temperatures, extending primers to achieve desired Tm values, batch processing of primers, and exporting results.

### Features

- **Primer Optimization**: Select and optimize primers based on specific criteria such as melting temperature (Tm) and Tm differences between forward and reverse primers.
- **Consensus Sequence Generation**: Generate a consensus sequence from multiple DNA sequences to identify common nucleotides or patterns.
- **Genome Annotation Visualization**: Visualize annotations on genomic sequences, providing insights into the distribution and characteristics of genes or features.
- **Batch Primer Processing**: Facilitate the processing of primers in batches, applying different optimization criteria or steps for each batch.
- **Exporting Results**: Export optimized primers, their degeneracy scores, and genome annotations to text and Excel files for further analysis.

### Installation

Before running this script, ensure you have Python 3.x installed along with the following packages:
- Biopython
- gffutils
- matplotlib
- numpy
- pandas

You can install the required packages using pip.

### Usage

1. Prepare your DNA sequences in FASTA format and, if available.
2. Update the script's parameters to point to your input files and specify desired output locations.
3. Run the script from the command line.

### Key Classes and Functions

#### Class: `GenomeMarker`

##### `__init__(self, sequences=None, gtf_file=None)`
Initializes the GenomeMarker object with sequences and optionally a genome annotation file.

##### `generate_consensus_or_return_sequence(self)`
Generates a consensus sequence if multiple genome records are provided; otherwise, returns the single sequence.

##### `generate_consensus(self)`
Calculates a consensus sequence by analyzing the nucleotide composition at each position across all provided sequences.

##### `save_sequences_with_consensus_to_fasta(self, output_filename)`
Saves the consensus sequence and, if applicable, the original sequences to a FASTA file.

### Examples

An example usage scenario might involve generating a consensus sequence from several aligned DNA sequences and then visualizing the distribution of specific primer sequences along the consensus.

### Contributing

Contributions to improve the script or add new features are welcome. Please fork the repository and submit pull requests with your proposed changes. For reporting bugs or requesting features, please open an issue through the GitHub issue tracker.

### Citing and References

If you use this tool in your research, please cite it as follows:

> Enrique Caban Centeno. Genome Marker Analysis Tool. 2024. [URL to repository].

This tool builds upon the work and methodologies from [Reference 1], [Reference 2], etc.

### Contact

For further questions or collaborations, feel free to contact me at centenoenrique1963@gmail.com.

### License

This script is provided under the MIT License. See the LICENSE file for more details.

