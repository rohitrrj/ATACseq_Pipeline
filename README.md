# ATAC-seq Analysis Pipeline
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-2.7+-blue.svg)
![Bioinformatics](https://img.shields.io/badge/bioinformatics-ATAC--seq-brightgreen.svg)

A comprehensive pipeline for processing and analyzing ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) data. This pipeline automates the entire workflow from raw FastQ files to peak calling, including quality control and visualization steps.

## Table of Contents
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Structure](#output-structure)
- [Quality Control](#quality-control)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)

## Features
- Automated processing of paired-end ATAC-seq data
- Comprehensive quality control metrics
- Adapter trimming and read preprocessing
- Alignment to reference genome
- Peak calling with MACS2
- Generate genome browser visualization files
- TSS enrichment analysis
- Insert size distribution analysis
- Library complexity estimation

## Prerequisites
- Python 2.7+
- FastQC (v0.11.2+)
- Bowtie2 (v2.2.4+)
- MACS2 (v2.1.1+)
- Samtools (v1.3+)
- Picard Tools (v1.92+)
- BedTools (v2.25.0+)
- R (v3.2.2+)
- IGVtools (v2.3.3+)
- Preseq (v1.0.2+)
- NGSplot (v2.47+)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/ATACseq_Pipeline.git
cd ATACseq_Pipeline
```

2. Ensure all required modules are available:
```bash
module load python fastqc/0.11.2 bowtie/2.2.4 MACS2/2.1.1 java/latest
module load samtools/1.3 picard-tools/1.92 preseq/1.0.2 bedtools/2.25.0
module load r/3.2.2 igvtools/2.3.3 ngsplot/2.47
```

3. Configure your project:
```bash
cp conf.txt.example conf.txt
# Edit conf.txt with your project-specific paths
```

## Usage

1. Prepare your sample file (one sample per line):
```bash
echo "sample1" > samples.txt
echo "sample2" >> samples.txt
```

2. Structure your input data:
```
/path/to/data/
├── sample1_R1.fastq
├── sample1_R2.fastq
├── sample2_R1.fastq
└── sample2_R2.fastq
```

3. Run the pipeline:
```bash
./ATAC_seq_preprocess.sh
```

## Pipeline Steps

1. **Quality Control** (FastQC)
   - Raw read quality assessment
   - Adapter content analysis
   
2. **Read Preprocessing**
   - Adapter trimming
   - Quality filtering
   
3. **Alignment** (Bowtie2)
   - Mapping to reference genome
   - Remove unmapped/low-quality reads
   
4. **Post-alignment Processing**
   - Remove duplicates
   - Remove mitochondrial reads
   - Filter for chromosomes 1-22
   
5. **Peak Calling** (MACS2)
   - Identify accessible regions
   - Generate signal tracks

## Output Structure
```
Analysis/
├── trim/                  # Trimmed FastQ files
├── fastqc/                # Quality control reports
├── bowtie/                # Alignment files
├── QC/                    # QC metrics and plots
└── Peaks/                 # MACS2 peak calls
```

## Quality Control
The pipeline generates several QC metrics:
- TSS enrichment scores
- Insert size distributions
- Library complexity estimates
- Alignment statistics
- Peak quality metrics

## Configuration
Edit `conf.txt` to specify:
- Data directory
- Project directory
- Scripts directory
- Sample file location

Example configuration:
```bash
myDATADIR="/path/to/data"
myPROJDIR="/path/to/project"
mySCRIPTSDIR="/path/to/scripts"
mySampleFile="samples.txt"
```

## Troubleshooting

### Common Issues

1. **Module Load Errors**
   - Ensure all required modules are available on your system
   - Check module versions match requirements

2. **Memory Issues**
   - Adjust -l h_vmem in script header
   - Split processing into smaller batches

3. **Missing Reference Files**
   - Verify reference genome paths
   - Check file permissions

### Error Messages

- `Error: cannot find input file` - Check file paths in conf.txt
- `Segmentation fault` - Increase memory allocation
- `Command not found` - Ensure all modules are loaded

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this pipeline in your research, please cite:
```
Author et al. (Year). ATACseq Pipeline: A comprehensive processing pipeline for ATAC-seq data.
Repository: https://github.com/yourusername/ATACseq_Pipeline
```

## Contributing
Contributions are welcome! Please read the contributing guidelines before submitting pull requests.
