
# snps_mags Pipeline

**snps_mags** is a Nextflow pipeline designed for **calling point mutations in Metagenome-Assembled Genomes (MAGs)** using **InStrain**. The pipeline integrates state-of-the-art tools such as **BWA-MEM2** for read alignment, **SAMtools** for BAM file processing, **Qualimap** for quality assessment, and **MultiQC** for comprehensive reporting. It is optimized for high-throughput analysis and can be executed in various environments using configurable execution profiles.

---

## Key Features
- **Point Mutation Detection**: Utilizes **InStrain** to identify single nucleotide variants (SNVs) in MAGs.
- **End-to-End Workflow**: From read alignment to variant calling and quality control.
- **Containerized Execution**: Leverages **Singularity** containers for reproducibility and portability.
- **Scalable**: Supports execution on **SLURM** clusters for large-scale analyses.
- **Comprehensive Reporting**: Generates detailed reports using **MultiQC** for easy interpretation of results.

---

## Overview

The pipeline performs the following tasks:

1. **BWA Indexing**: Indexes the reference genome using `bwa-mem2`.
2. **SAMtools Indexing**: Indexes the reference genome using `samtools`.
3. **Alignment**: Aligns paired-end reads to the reference genome using `bwa-mem2`, sorts the output, and marks duplicates.
4. **Qualimap**: Runs Qualimap to assess the quality of the aligned BAM files.
5. **InStrain Variant Calling**: Performs variant calling using InStrain.
6. **MultiQC**: Generates a MultiQC report summarizing the results.

---

## Pipeline Processes

### 1. `bwa_index`
- **Description**: Indexes the reference genome using `bwa-mem2`.
- **Input**: Reference genome file.
- **Output**: BWA index files .
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

### 2. `samtools_index`
- **Description**: Indexes the reference genome using `samtools`.
- **Input**: Reference genome file.
- **Output**: SAMtools index file (`.fai`).
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

### 3. `aln_pipe`
- **Description**: Aligns paired-end reads to the reference genome using `bwa-mem2`, sorts the output, and marks duplicates.
- **Input**: 
  - Sample name and paired-end reads.
  - BWA index files.
- **Output**: 
  - Marked BAM file (`.marked.bam`).
  - BAM index file (`.marked.bam.bai`).
  - Log file (`.log.bwamem`).
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

### 4. `qualimap`
- **Description**: Runs Qualimap to assess the quality of the aligned BAM files.
- **Input**: 
  - Sample name and marked BAM file.
  - Reference genome file.
- **Output**: Qualimap results directory.
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

### 5. `instrain_variant_calling`
- **Description**: Performs variant calling using InStrain.
- **Input**: 
  - Sample name and marked BAM file.
  - Reference genome file.
  - SAMtools index file (`.fai`).
- **Output**: InStrain results directory.
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

### 6. `multiqc`
- **Description**: Generates a MultiQC report summarizing the results.
- **Input**: All output files from previous processes.
- **Output**: MultiQC report (`.html`).
- **Container**: `oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f`.

---

## Parameters

The pipeline is configured using the following parameters:

| Parameter       | Default Value          | Description                                                                 |
|-----------------|------------------------|-----------------------------------------------------------------------------|
| `reads`         | `"data/*_{1,2}.fastq.gz"` | Path to input paired-end reads.                                             |
| `reference`     | `"reference.fa"`       | Path to the reference genome file (set of MAGs).                                          |
| `outdir`        | `"results"`            | Directory to store pipeline outputs.                                        |
| `rlibrary`      | `"pair-ends"`          | Read group library name.                                                    |
| `rplat`         | `"ILL"`                | Read group platform.                                                        |
| `debug`         | `true`                 | Debug mode. If `true`, commands are printed but not executed.               |

---

## Resource Allocation

Resource allocation for specific processes can be configured as follows:

| Process                     | Memory  | CPUs |
|-----------------------------|---------|------|
| `aln_pipe`                  | 5.GB    | 2    |
| `qualimap`                  | 5.GB    | 4    |
| `instrain_variant_calling`  | 5.GB    | 4    |
| `multiqc`                   | 5.GB    | 1    |

---

## Execution Profiles

The pipeline supports multiple execution profiles for different environments:

### 1. `kutral`
- **Executor**: SLURM
- **Queue**: `ngen-ko`
- **Queue Size**: 10
- **Container**: Singularity
- **Bind Path**: `/mnt/beegfs/labs/`

### 2. `leftraru`
- **Executor**: SLURM
- **Queue**: `slims`
- **Queue Size**: 200
- **Container**: Singularity

---

## Execution Reports

The pipeline generates the following execution reports:
- **Timeline**: HTML file showing the execution timeline.
- **Report**: HTML file summarizing the execution.
- **Trace**: Text file containing execution trace information.
- **DAG**: HTML file showing the pipeline's Directed Acyclic Graph (DAG).

---


## Running the Pipeline

To run the pipeline, use the following command:

```bash
nextflow run pipeline.nf -profile <profile_name>
```

Replace <profile_name> with the desired execution profile (kutral or leftraru).

---

## Debug Mode
To run the pipeline in debug mode (commands are printed but not executed), set the debug parameter to true:

```
nextflow run pipeline.nf -profile <profile_name> --debug true
```
---

## Output Directory Structure

The pipeline outputs are organized as follows:

```
results/
├── bwa_index/
├── samtools_index/
├── alignment_results/
├── qualimap/
├── instrain/
├── multiqc/
└── pipeline_info/
    ├── execution_timeline_<timestamp>.html
    ├── execution_report_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    └── pipeline_dag_<timestamp>.html
```
---

## Notes
Ensure that the input files (reads and reference) are correctly specified in the params section.

Adjust resource allocation (memory and CPUs) based on your system's capabilities.

---

## Manifest

| Key            | Value                                      |
|----------------|--------------------------------------------|
| `defaultBranch`| `main`                                     |
| `homePage`     | `https://github.com/digenoma-lab/mag_snps` |
| `author`       | `Alex Di Genova`                           |
| `version`      | `0.0.1`                                    |

