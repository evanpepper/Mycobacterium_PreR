# Host oxidative stress primes mycobacteria for rapid antibiotic resistance evolution

## Abstract

The rapid emergence of multidrug-resistant *Mycobacterium tuberculosis*
(Mtb) threatens global TB control, yet the mechanisms enabling rapid
evolution of drug resistance in Mtb remain poorly understood. Here we
reveal that pre-existing mutations in oxidative stress response (OSR)
genes create permissive genomic backgrounds that accelerate high-level
isoniazid resistance (INHR) without fitness costs, challenging the
paradigm that resistance mutations always precede their fitness
compensatory adaptations. Using *M. smegmatis* mc²155 (Msm) as a model,
we show that brief exposure to sublethal INH (2× IC₅₀) enriches for
"low-level resistance and tolerance" (LLRT) mutants in a single step.
These LLRT mutants, particularly those with **ohrR** loss-of-function
mutations, acquire high-level resistance (\>500× IC₅₀) at 6-fold higher
rates than wildtype, primarily through otherwise-deleterious mycothiol
biosynthesis mutations that become tolerable in the
oxidative-stress-buffered background.

Crucially, we demonstrate that sublethal oxidative stress alone,
mimicking host immune pressure, nearly tripled the rate of INH
resistance evolution in Msm. Bayesian analysis of 1,578 clinical Mtb
isolates from Vietnam confirmed that mutations in oxidative stress
response genes were significantly associated with the emergence of INHR
strains (p-value = 1.09 × 10⁻⁷). Independently, reanalysis of
genome-wide CRISPRi screens revealed that the OSR network and high Bayes
probability genes are functionally associated with treatment escape and
survival under multiple antibiotics, including isoniazid, rifampicin,
ethambutol, bedaquiline, vancomycin, clarithromycin, linezolid, and
streptomycin.

Our findings that host-imposed oxidative stress and inadequate drug
penetration may synergistically prime Mtb populations for rapid
resistance evolution suggest that targeting pre-resistance mechanisms,
such as oxidative stress defenses, could help slow the emergence of
antibiotic resistance in tuberculosis.

## Repository Overview

This repository contains all computational and analytical components
supporting the manuscript:

**Host oxidative stress primes mycobacteria for rapid antibiotic
resistance evolution**

### Directory structure

    manuscript/
    ├── code/          
    ├── data/          
    ├── figures/       
    ├── supplemental-figures/
    ├── supplemental-tables/
    ├── LICENSE       
    └── README.md      

## How to reproduce the entire analysis performed in this manuscript

1.  Clone the repository:

    ``` bash
    git clone https://github.com/evanpepper/Mycobacterium_PreR.git
    cd Mycobacterium_PreR/code
    ```

2.  Install dependencies for each R Markdown file, found at the beginning of each script.

3.  Run the R Markdown file, either cell by cell or all at once. Figures will be deposited into `figures/` and generated results will appear in `data/`

## Citation

*(DOI)*

## Contact

**Evan Pepper** -- epepper@systemsbiology.org or epepper@uw.edu
