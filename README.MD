# Quality Control Pipeline for Sequencing Data

## Overview
This script automates the quality control (QC) process for sequencing data, including:
- Raw data quality checks
- Cleaning and trimming reads using **Fastp**
- Post-cleaning QC
- Generating comprehensive summary reports

The updated version includes a dedicated **reports** folder for Fastp logs and the read summary file.

---

## Requirements

### Software Dependencies
Ensure the following tools are installed and accessible in the system PATH:
- **fastqc**
- **fastp**
- **multiqc**

### Conda Environment
The script expects the **`tools_qc`** Conda environment to be active. If not, it will attempt to activate it automatically.

You can create the environment using:

```bash
conda create -n tools_qc fastqc fastp multiqc -c bioconda -c conda-forge
```

---

## Sample Naming Convention
The script identifies sequencing files based on their naming patterns. Ensure your files follow one of the following conventions:

### Paired-end Files
- `<sample_name>_1.fastq.gz` and `<sample_name>_2.fastq.gz`
- `<sample_name>R1.fastq.gz` and `<sample_name>R2.fastq.gz`
- `<sample_name>_1_001.fastq.gz` and `<sample_name>_2_001.fastq.gz`

### Single-end Files
- `<sample_name>.fastq.gz`
- `<sample_name>.fastq`

### Handling Missing or Mismatched Pairs
- **Missing Pair:** If a paired-end file (e.g., `_2.fastq.gz`) is not found for a given `_1.fastq.gz` file, the script will log a warning and process the `_1.fastq.gz` file as single-end data.
- **Mismatched Pair:** If paired files have inconsistent prefixes or formats, they will be skipped, and the script will log an error. Ensure proper naming conventions to avoid such issues.

---

## Input and Output Directory Setup

### Input Files
Place all input files (FASTQ/FASTQ.GZ) in the **same directory** as the script or specify the path in the `INPUT_DIR` variable in the script.

### Results
The script creates a structured results directory:

```plaintext
results/
├── fastqc/
│   ├── before_cleaning/
│   └── after_cleaning/
├── fastp_cleaned/
├── multiqc/
├── reports/
│   ├── <sample_name>_fastp.log
│   └── read_summary.csv
```

---

## Running the Script

### 1. Preparation
- Place the script and all raw FASTQ files in the same directory, or modify the `INPUT_DIR` variable in the script to point to the directory containing the files.
- Activate the `tools_qc` environment:

```bash
conda activate tools_qc
```

### 2. Execute the Script
Run the script using:

```bash
bash fastp_pipeline.sh
```

### 3. Parameter Configuration
The script provides default parameters for Fastp, which are displayed at the start. You can:
- Accept the defaults by typing `yes`.
- Modify parameters interactively by typing `no`.

---

## Outputs

### 1. **FastQC Reports**
- Raw data reports are in `results/fastqc/before_cleaning/`.
- Cleaned data reports are in `results/fastqc/after_cleaning/`.

### 2. **Cleaned FASTQ Files**
- Saved in `results/fastp_cleaned/` with `_cleaned_1.fastq.gz` and `_cleaned_2.fastq.gz` suffixes for paired-end data.

### 3. **MultiQC Report**
- Combined quality report saved in `results/multiqc/`.

### 4. **Reports**
- Fastp logs for each sample are saved in `results/reports/`.
- `read_summary.csv` contains statistics for each sample:

```csv
Sample_Name,Initial_Reads,Processed_Reads,Discarded_Reads
```

---

## Error Handling
- **Missing Tools:** If required tools (e.g., fastqc, fastp) are not found in the system PATH, the script will log an error and terminate. Install the missing tools and re-run the script.
- **Conda Environment:** If the `tools_qc` environment is not active and cannot be activated automatically, the script will log an error and terminate. Activate the environment manually before running the script.

---

## Example of Directories Setup

### Directory Structure Before Running:

```plaintext
project_dir/
├── fastp_pipeline.sh
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
├── sample2_2.fastq.gz
```

### Directory Structure After Running:

```plaintext
project_dir/
├── fastp_pipeline.sh
├── results/
│   ├── fastqc/
│   │   ├── before_cleaning/
│   │   └── after_cleaning/
│   ├── fastp_cleaned/
│   │   ├── sample1_cleaned_1.fastq.gz
│   │   ├── sample1_cleaned_2.fastq.gz
│   │   ├── sample2_cleaned_1.fastq.gz
│   │   └── sample2_cleaned_2.fastq.gz
│   ├── multiqc/
│   │   └── multiqc_report.html
│   ├── reports/
│   │   ├── sample1_fastp.log
│   │   ├── sample2_fastp.log
│   │   └── read_summary.csv
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
├── sample2_2.fastq.gz
