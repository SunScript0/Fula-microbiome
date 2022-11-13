This repository contains the pipeline and analysis scripts used in the paper "Urbanization restricts gut microbiome diversity and delays its maturation in infants"

Intructions:
1. Install and activate qiime2. The version used in the paper is qiime2-2021.11.
2. Download fastq files and place them in 01_fastqs (optionally verify checksums)
3. Edit manifest.tsv to have the absolute paths of all the fastq files
4. Run "run_qiime2.sh", outputs will be in 02_pipeline_outputs
5. Run "prepro_fula.R" (Might require installing some R packages first), outputs will be in 03_processed_data
6. Run "analysis_fula.R" (Might require installing some R packages first), outputs will be in 04_analysis_output

Note: 03_processed_data is pre-populated in this repository, so it is possible to only run the analysis step ("analysis_fula.R")
