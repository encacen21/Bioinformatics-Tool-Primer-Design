# GenomeMarker: Primer Design and Consensus Sequence Workflow

GenomeMarker is a Python-based bioinformatics workflow for generating consensus sequences from aligned genomic FASTA files, designing overlapping amplicon primer schemes, estimating simple primer melting temperature values, and exporting analysis-ready outputs.

This repository was originally developed as part of a Master's project in Bioinformatics. It is presented here as a portfolio example of sequence processing, primer-design logic, reproducible scientific scripting, and bioinformatics workflow organization.

## Overview

The workflow supports primer-design tasks from aligned FASTA sequences. It can generate a consensus sequence from one or more aligned genomic records, handle IUPAC ambiguity codes, place primer pairs across a target sequence, evaluate primer properties, support optimization logic, estimate primer degeneracy, and prepare outputs that can be reused in downstream analysis or reporting.

This project is intended as an educational and research-oriented workflow example. Primer candidates generated computationally should always be reviewed with dedicated primer-design software and validated experimentally before laboratory use.

## Main Features

- Consensus sequence generation from aligned FASTA records
- IUPAC ambiguity handling for variable genomic positions
- Overlapping amplicon primer-design logic
- Forward and reverse primer placement
- Basic melting temperature estimation
- Primer quality checks and optimization logic
- Primer degeneracy scoring
- Sequence visualization and reporting support
- Export of analysis-ready outputs

## Technologies Used

- Python
- Biopython
- pandas
- Matplotlib
- ReportLab
- openpyxl

## Repository Structure

```text
Bioinformatics-Tool-Primer-Design/
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
├── src/
│   └── genome_marker.py
├── examples/
│   ├── README.md
│   ├── example_alignment.fasta
│   └── expected_output/
└── docs/
```

## Installation

Clone the repository:

```bash
git clone https://github.com/encacen21/Bioinformatics-Tool-Primer-Design.git
cd Bioinformatics-Tool-Primer-Design
```

Create and activate a Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Example Data

A small synthetic aligned FASTA file is provided in:

```text
examples/example_alignment.fasta
```

This mock file is included only to demonstrate the expected input format. It is not derived from unpublished research data.

## Example Usage

The current script is a research prototype and may require manual input depending on the workflow configuration. A cleaned command-line interface is planned as part of future improvements.

The main workflow code is located in:

```text
src/genome_marker.py
```

## Expected Outputs

Depending on the selected workflow configuration, the script can generate outputs such as:

- consensus FASTA sequences
- primer tables
- primer batch summaries
- primer degeneracy scores
- Excel summary files
- sequence or primer visualization outputs
- report-ready files for downstream review

## Limitations

This repository is a portfolio and research workflow example, not a replacement for experimental primer validation.

Current limitations include:

- simple rule-based melting temperature estimation
- no full thermodynamic primer-specificity model
- no complete off-target search workflow
- no wet-lab validation
- assumes aligned input sequences
- may require manual review of ambiguous or gapped regions
- currently structured as a research prototype rather than a production package

## Future Improvements

Planned improvements include:

- expanding the mock FASTA example with expected output files
- creating a command-line interface with argparse
- adding lightweight tests for core functions
- improving modular structure
- adding workflow diagrams and usage examples

## Author

Enrique A. Caban Centeno  
PhD researcher and bioinformatics analyst

Focus areas: RNA/miRNA analysis, sequence analysis, biological annotation, qPCR data processing, reproducible Python/R workflows, and omics-style data analysis.

## License

This project is released under the MIT License. See the LICENSE file for details.
