# `mar`: Mutations-Area Relationship Analysis

<!-- badges: start -->
![GitHub R package version](https://img.shields.io/github/r-package/v/meixilin/mar)
[![R-CMD-check](https://github.com/r-lib/gert/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-lib/gert/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green)
<!-- badges: end -->

## Overview
`mar` is an R package that enables the reconstruction of Mutations-Area Relationships using spatially distributed genome variation data. This tool helps researchers analyze how genetic mutations accumulate across geographic space within a species.

## Installation

To install the development version with all dependencies:

```R
library(devtools)
devtools::install_github("meixilin/mar")
```

Or in bash:

```bash
R CMD INSTALL --preclean --no-multiarch --with-keep.source mar
```

### Troubleshooting

If you encounter issues with installation above, try manually installing these dependencies before installing `mar`:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SeqArray")
install.packages(c("raster", "sars", "sads", "matrixStats"))
```

Or using `conda` / `mamba`:

```bash
mamba create --yes --name mar r-base r-devtools zlib libzlib gdal r-terra r-curl bioconda::bioconductor-seqarray
mamba activate mar
R
library(devtools)
devtools::install_github("meixilin/mar")
```

Use of `SeqArray` version >= 1.28.0 is **required** (there was a bug in the previous version that impacts plink file importing).

If the above tricks do not work, please open an issue on GitHub.

## Minimal working example

### MARPIPELINE

For running `MARPIPELINE`, example genotype and coordinate files are provided in `inst/extdata` folder.
This is a dummy data for demonstration purposes. To run the example, simply run:

```R
library(mar)
library(raster)

genofile <- system.file("extdata", "genome.tsv", package = "mar")
lonlatfile <- system.file("extdata", "lonlat.csv", package = "mar")

MARPIPELINE(name = "example", workdir = "./", genofile = genofile, lonlatfile = lonlatfile, saveobj = TRUE)
```

### Step-by-step analysis

For running individual steps, an example `genomaps` object is provided as `mar::gm1001g`.

Main analysis steps can be run step by step:

```R
library(mar)

# Fit SAD models
marsad <- MARsad(gm = gm1001g, sad_models = c("lnorm", "ls"), folded = TRUE)

# Sample to build MAR
mardf <- MARsampling(gm = gm1001g, scheme = "random")

# Calculate MAR
MARcalc(mardf = mardf)

# Simulate extinction
extdf <- MARextinction(gm = gm1001g, scheme = "random")

# Calculate extinction-based MAR
MARcalc(mardf = extdf)
```

## Preparing input data

The `mar` package performs extensive checks on the input data to ensure that the results are reliable and interpretable.
Please make sure to check the `inst/extdata` folder for example files.

Here are the list of checks that are performed, and some tips on how to prepare the data:

### Genotype Data

- The pipeline requires bi-allelic SNP genotype data.
    - Example `bcftools` command: `bcftools view -m2 -M2 -v snps ${VCF}`
- The genotype data must not contain any missing values.
    - You can use tools like [beagle](https://faculty.washington.edu/browning/beagle/beagle.html) to impute missing data.
    - You can also filter the genotype data to retain only sites without missing data, e.g., `bcftools view -i 'N_MISSING == 0' ${VCF}`.
- It works best with diploid genotype data, but can handle any ploidy as long as it is consistent across all samples. Use caution when interpreting results for non-diploid data.
    - If heterozygous genotypes are not confidently called, you can force the data to be haploid (`option_geno$ploidy = 1`) by converting heterozygous genotypes to the alternative allele. Use caution when interpreting results.
- If the reference genome is divergent from the species/population of interest, set the major allele as the reference allele to avoid issues with ancestral state identification.

### Geographic data

- Every sample with genotype data must have pairing longitude and latitude data.
- File format should be tab-delimited or comma-delimited with header (ID, LON/LONGITUDE, LAT/LATITUDE)
- Only three columns are allowed (ID, LON/LONGITUDE, LAT/LATITUDE)
- Sample IDs must be unique and in the same order as the Sample IDs provided in the genotype matrix.
