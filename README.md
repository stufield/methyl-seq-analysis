
# Overview

A quick directory structure and setup for the Methyl Seq Analysis at ARBL.

------------------

## Main files and directories

- `R/`: contains `R` code
  - `lib.R`: contains user defined functions to be used during analysis
  - `setup.R`: contains package installation instructions (must only be run once)
  - `run.R`: contains code that runs the code in `lib.R` once `setup.R`
    has been installed fully
- `data/`: all source data to be used during analysis
  - `*.bam`: source DNA sequence files stored as `*.bam` files
  - `*.csv`: any clinical data that needs to be merged (regarding samples?)
- `output/`: generated output
  - `plots/`: contains any plots that may be generated by `run.R`
  - `tables/`: contains any `*.csv`, `*.txt`, or table-like output from `run.R`
- `reports/`: 
  - contains analysis output as reports (e.g. results from `*.Rmd` files)


## Example Structure

````markdown
|--methyl-seq-analysis.Rproj    # RStudio file marks project root
|--R/
    |--setup.R                  # everything necessary to recreate analysis environment
    |--lib.R                    # analysis-specific function library
    |--run.R                    # calls `lib.R`; runs & generates all plots/tables
|--data/
    |--dna-seq-01.bam           # source DNA sequence data
    |--dna-seq-02.bam           # source DNA sequence data
    |--clinical_meta.csv        # source data
    |--training-data.rds        # derived from `*.bam`
|--output/
    |--plots/
        |--figure1_2018-11-09.pdf  # figs resulting from run.R
        |--figure2_2018-11-09.png
    |--tables/
        |--result_table1.csv       # tables resulting from run.R
        |--result_table2.csv
    |--models/
        |--output_model.rds        # any models used
|--reports/
    |--analysis-report.Rmd      # final report in Rmarkdown
                                # refer `source()` `lib.R` in setup chunk
    |--analysis-report.html     # rendered Rmarkdown in `html`
    |--analysis-report.pdf      # rendered Rmarkdown in `pdf`
|--explore/
    |--dead-end-analyses
    |--go-here
|--misc/
    |--reference_paper.pdf      # project-specific free-from wild-card
    |--analysis_plan.docx       # files you don't want polluting the root
````

## General notes

- code is designed to run from the project root (when cloned)
- methyl-seq data goes in the `data/` directory
- `R` code goes in the `R/` directory
- Additional reports (e.g. `*.Rmd`) files go in the `reports/` directory
- relative paths can be used when executing from within the project
