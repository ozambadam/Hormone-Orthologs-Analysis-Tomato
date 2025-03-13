# Hormone-Orthologs-Analysis-Tomato
This repository contains the code and data files used for the analysis of orthologs related to hormones in tomato. The process includes integrating data from different species, identifying orthologous genes involved in auxin and brassinosteroid pathways, and generating output files with information 
## Project Description

This project involves analyzing gene orthologs related to two plant hormones, **Brassinosteroids** and **Auxins**, in tomato. It involves two primary steps:

1. **Data Collection**: Gather orthologous genes across different species, especially focusing on tomato and its related species.
2. **Data Processing and Analysis**: Use R to join, clean, and manipulate orthologous data to identify common genes. The analysis includes:
   - Merging data with orthologs.
   - Identifying duplicates or genes found multiple times across different species.
   - Processing hormone-specific genes.
   - Generating output files for further analysis.

## Files in this Repository

- `Hormone-Ortholog-Analysis.R`: The main R script containing code for data processing, merging datasets, and analyzing hormone-specific orthologs.
- `BrassinoMipedit.txt`: Input file containing data on Brassinosteroid-related genes in tomato.
- `AuxinMipedit.txt`: Input file containing data on Auxin-related genes in tomato.
- `BrassinoMipOrthoSpecies.txt`: Output file with orthologous gene data related to Brassinosteroids.
- `AuxMiPOrthoSpecies.txt`: Output file with orthologous gene data related to Auxins.
- `mip_net.txt`: Output network file for Brassinosteroids.
- `AuxMip_net.txt`: Output network file for Auxins.
- `duplicated_genesSpeciesBrassinoMip.txt`: File with duplicated genes related to Brassinosteroids across species.
- `duplicated_genesSpeciesAuxMip.txt`: File with duplicated genes related to Auxins across species.

## Key Steps

1. **Reading the Data**:
   - Load the datasets that contain the information on hormone-specific genes and orthologs.

2. **Merging Data**:
   - Use `left_join` to merge the datasets with orthologs data, ensuring that all relevant genes are captured across species.

3. **Data Transformation**:
   - Clean and modify the data to replace "tomato" with the appropriate species in the `Species_b` column.
   - Normalize the hormone-specific data by setting all values in `Species_a` to "tomato".

4. **Network Analysis**:
   - Merge the gene data with the corresponding network data to build a complete representation of the orthologous relationships between species.

5. **Identifying Duplicates**:
   - Group the data by gene and species to identify duplicates—genes that appear in multiple species.

6. **Final Output**:
   - Save the processed data into output files for further analysis or reporting.

## Installation

To run the analysis, make sure you have R installed on your system. You will also need to install the required R packages. You can install them using the following commands:

```R
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("ggalluvial")

```

## Usage
After installing the necessary packages, you can run the script Hormone-Ortholog-Analysis.R. Make sure the input files (BrassinoMipedit.txt, AuxinMipedit.txt, etc.) are in the correct directory. The script will generate output files containing orthologous genes and associated network data.
