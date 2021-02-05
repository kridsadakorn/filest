<!-- README.md is generated from README.Rmd. Please edit that file -->

# FILEST <img src="man/figures/filest_logo.png" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/kridsadakorn/filest/workflows/R-CMD-check/badge.svg)](https://github.com/kridsadakorn/filest/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/FILEST)](https://CRAN.R-project.org/package=FILEST)
[![codecov](https://codecov.io/gh/kridsadakorn/filest/branch/master/graph/badge.svg?token=41WJVZOP49)](https://codecov.io/gh/kridsadakorn/filest)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.com/kridsadakorn/filest.svg?branch=master)](https://travis-ci.com/kridsadakorn/filest)
[![DOI](https://zenodo.org/badge/333290990.svg)](https://zenodo.org/badge/latestdoi/333290990)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fkridsadakorn%2Ffilest.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2Fkridsadakorn%2Ffilest?ref=badge_shield)
<!-- badges: end -->

`FILEST` (Fine-Level Structure Simulator) is a population genetic
simulator. The simulator is able to generate synthetic datasets for
single-nucleotide polymorphisms (SNP) for multiple populations. The
genetic distances among populations can be set according to the Fixation
Index (Fst). This tool is able to simulate outlying individuals and
missing SNPs can be specified. For Genome-wide association study (GWAS),
disease status can be set in desired level according risk ratio.

## Installation

You can install the released version of FILEST from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("FILEST")
```

Alternatively, you can install the dev version of FILEST from
[Github](https://github.com/kridsadakorn/filest) with

``` r
install.packages("remotes")
remotes::install_github("kridsadakorn/filest", dependencies = TRUE)
```

## Document

You can see the reference manual from:
<https://www.biostatgen.org/filest/>

## Example

This best way to understand the package is the run the demo function and
edit the parameters in the setting file:

``` r
library(FILEST)

output_dir <- demo.filest()
#> Creating a setting file ... /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//Rtmp2cYXi6/example1.txt
#> Generating the simulated data  to  ... /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//Rtmp2cYXi6
#> Start [S0] at 2021-02-04 03:13:27
#> Setting file is : /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//Rtmp2cYXi6/example1.txt
#> The simulated files will be saved in this directory: /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//Rtmp2cYXi6/example1
#> Creating data file setting #1 - rep #1
#> Done - 0.0865130424499511719 secs
#> Writing data files setting #1 - rep #1
#> Done - 0.293856143951416016 secs
#> Creating status file setting #1 - rep #1
#> Done - 0.295050144195556641 secs
#> Estimating Fst setting #1 - rep #1
#> Done - 0.317656993865966797 secs
#> Creating maker information setting #1 - rep #1
#> Done - 0.559202194213867188 secs
#> Generating PC scores #1 - rep #1
#> Generating EigenVector  #1 - rep #1
#> Done - 1.19890999794006348 secs
```

The demo function creates the setting file at a temp directory as
`example1.txt`, and you can use as a template to edit:

``` r
print(file.path(output_dir,"example1.txt"))
#> [1] "/var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//Rtmp2cYXi6/example1.txt"
```

The demo function create the simulated files at a temp directory as
\`well.`example1`.

``` r
dir(file.path(output_dir,"example1"))
#>  [1] "simSNP_rep1_data_numMark_rowInd_colVar.txt"
#>  [2] "simSNP_rep1_data_numMark_rowVar_colInd.txt"
#>  [3] "simSNP_rep1_eigenvector10.txt"             
#>  [4] "simSNP_rep1_estimated_Fst.txt"             
#>  [5] "simSNP_rep1_individuals_with_header.txt"   
#>  [6] "simSNP_rep1_individuals.txt"               
#>  [7] "simSNP_rep1_PC.pdf"                        
#>  [8] "simSNP_rep1_PC10.txt"                      
#>  [9] "simSNP_rep1.bed"                           
#> [10] "simSNP_rep1.bim"                           
#> [11] "simSNP_rep1.fam"                           
#> [12] "simSNP_rep1.RData"
```

## About

  - Prof. Kristel Van Steen, visit
    <a href="http://bio3.giga.ulg.ac.be/" border=0 style="border:0; text-decoration:none; outline:none"><img width="40px" src="man/figures/bio3_logo.png" align="center" /></a><br />
  - Kridsadakorn Chaichoompu, visit
    <a href="https://www.biostatgen.org/" border=0 style="border:0; text-decoration:none; outline:none"><img width="110px" src="man/figures/biostatgen_logo.png" align="center" /></a><br />


## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fkridsadakorn%2Ffilest.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2Fkridsadakorn%2Ffilest?ref=badge_large)