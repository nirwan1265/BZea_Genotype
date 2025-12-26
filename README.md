# BZea Genotyping with ANGSD

This repository provides a comprehensive guide for genotyping maize (Zea mays, specifically B73 and related lines) using ANGSD (Analysis of Next Generation Sequencing Data).

## Overview

ANGSD is a powerful software for analyzing Next Generation Sequencing (NGS) data that can handle genotype uncertainty. Unlike traditional genotype callers, ANGSD uses genotype likelihoods rather than called genotypes, making it particularly useful for:

- Low to medium coverage sequencing data
- Ancient DNA or degraded samples
- Population genetic analyses
- Estimating allele frequencies and genetic diversity

This repository contains tutorials, example scripts, and workflows for performing genotyping analysis on maize samples using ANGSD.

## Contents

- **[Installation Guide](docs/installation.md)** - How to install ANGSD and dependencies
- **[Data Preparation](docs/data_preparation.md)** - Preparing your sequencing data for analysis
- **[Genotyping Workflow](docs/genotyping_workflow.md)** - Step-by-step genotyping pipeline
- **[Quality Control](docs/quality_control.md)** - QC procedures and best practices
- **[Analysis & Interpretation](docs/analysis.md)** - Interpreting results and downstream analyses
- **[Example Scripts](scripts/)** - Ready-to-use bash scripts for common tasks

## Quick Start

1. Install ANGSD following the [installation guide](docs/installation.md)
2. Prepare your BAM files as described in [data preparation](docs/data_preparation.md)
3. Run the basic genotyping workflow using scripts in the `scripts/` directory
4. Perform quality control and analyze results

## Prerequisites

- Basic knowledge of command line/bash
- Understanding of NGS data and BAM file format
- Familiarity with population genetics concepts (helpful but not required)
- Access to aligned sequencing data (BAM files)

## Repository Structure

```
.
├── docs/                    # Detailed documentation
├── scripts/                 # Example scripts and workflows
├── examples/                # Example data and outputs
└── README.md               # This file
```

## Citation

If you use ANGSD in your research, please cite:

> Korneliussen, T. S., Albrechtsen, A., & Nielsen, R. (2014). ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics, 15, 356.

## Resources

- [ANGSD Official Website](http://www.popgen.dk/angsd/index.php/ANGSD)
- [ANGSD GitHub Repository](https://github.com/ANGSD/angsd)
- [Maize Genetics and Genomics Database](https://www.maizegdb.org/)

## Contributing

Contributions to improve this tutorial are welcome! Please feel free to submit issues or pull requests.

## License

See [LICENSE](LICENSE) file for details.
