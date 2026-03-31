# Host oxidative stress primes mycobacteria for rapid antibiotic resistance evolution

## Abstract

The rapid emergence of multidrug-resistant *Mycobacterium tuberculosis* (Mtb) threatens global TB control, yet the mechanisms enabling rapid evolution of resistance in Mtb remain poorly understood. Here, we show that pre-existing mutations in oxidative stress response genes create permissive genomic backgrounds that accelerate high-level isoniazid resistance (INH<sup>R</sup>), challenging the paradigm that resistance mutations must precede compensatory adaptation. Using *Mycobacterium smegmatis* mc<sup>2</sup>155 (Msm) as a model, we demonstrate that brief exposure to sublethal isoniazid (INH) enriches for “low-level resistance and tolerance” (LLRT) mutants in a single step. LLRT mutants, particularly those with *ohrR* loss-of-function mutations, acquire high-level resistance (>500× IC<sub>50</sub>) at ~6-fold higher rates than wildtype, primarily through otherwise deleterious mycothiol biosynthesis mutations that become tolerable in an oxidative stress-buffered background.

Crucially, sublethal oxidative stress alone, mimicking host immune pressure, nearly tripled the rate of INH<sup>R</sup> evolution. Analysis of 1,578 clinical Mtb isolates revealed significant enrichment of oxidative stress-related loci among those associated with INH<sup>R</sup>. Reanalysis of genome-wide CRISPRi data further linked oxidative stress response pathways to survival under multiple antibiotics. Together, these findings suggest that host-imposed oxidative stress and suboptimal drug exposure may prime Mtb populations for rapid resistance evolution, highlighting oxidative stress defenses as potential targets to limit resistance emergence.

## Repository Overview

This repository contains all code and data required to reproduce the analysis and figures in the associated <a href="https://doi.org/10.1101/2025.11.19.689367" target="_blank">manuscript</a>.

### Directory structure

    Mycobacterium_PreR/
    ├── code/          
    ├── data/          
    ├── figures/       
    ├── supplementary-figures/
    ├── supplementary-tables/
    ├── LICENSE       
    └── README.md      

## How to reproduce the entire analysis performed in this manuscript

1.  Clone the repository:

    ``` bash
    git clone https://github.com/evanpepper/Mycobacterium_PreR.git
    cd Mycobacterium_PreR/code
    ```

2.  Install dependencies for each R Markdown file, found at the beginning of each script.

3.  The Fig. 6 notebook will require downloading two datasets:

    Dataset 1: <a href="10.5281/zenodo.19207254" target="_blank">Bayes Null Distributions from Zenodo</a>
    Dataset 2: <a href="10.1038/s41564-022-01130-y" target="_blank">Li *et al* 2022, Source Data Fig. 1</a>

    Unzip Dataset 1:

    ``` bash
    cd Mycobacterium_PreR/ 
    mkdir Mycobacterium_PreR/data/Fig6/LDF/bayes-null-distributions/INH-global-download 	
    tar -xvf INH-global.tar.gz -C Mycobacterium_PreR/data/Fig6/LDF/bayes-null-distributions/INH-global-download
    ``` 

    Copy Dataset 2 into the following directory:
    ``` bash
    cd Mycobacterium_PreR/ 
    mkdir Mycobacterium_PreR/data/published-data
    mv 41564_2022_1130_MOESM4_ESM.xlsx li-et-al-2022-table-SD-F1.xlsx 	
    cp li-et-al-2022-table-SD-F1.xlsx /data/published-data/li-et-al-2022-table-SD-F1.xlsx
    ``` 

4.  Run the R Markdown file, either cell by cell or all at once. Figures will be deposited into `figures/` and generated results will appear in `data/`

## Dependencies

The analyses in this repository were performed using the following software and R packages:

- **R** (version 4.1.1)  
- **scales** (version 1.4.0)  
- **ggpubr** (version 0.6.1)  
- **rstatix** (version 0.7.2)  
- **naniar** (version 1.1.0)  
- **factoextra** (version 1.0.7)  
- **knitr** (version 1.50)  
- **drc** (version 3.0-1)  
- **lubridate** (version 1.9.4)  
- **forcats** (version 1.0.0)  
- **stringr** (version 1.5.1)  
- **dplyr** (version 1.1.4)  
- **purrr** (version 1.1.0)  
- **readr** (version 2.1.5)  
- **tidyr** (version 1.3.1)  
- **tibble** (version 3.3.0)  
- **tidyverse** (version 2.0.0)  
- **growthcurver** (version 0.3.1)  
- **ggridges** (version 0.5.7)  
- **ggplot2** (version 3.5.2)  
- **ggrepel** (version 0.9.6)  
- **ggupset** (version 0.4.1)  
- **see** (version 0.11.0)  
- **viridis** (version 0.6.5)  
- **seqinr** (version 4.2-36)  
- **Biostrings** (version 2.72.1)  

## Citation

<a href="https://doi.org/10.1101/2025.11.19.689367" target="_blank"><strong>Link to bioRxiv preprint</strong></a>

## Contact

**Evan Pepper** -- epepper@systemsbiology.org or epepper@uw.edu
