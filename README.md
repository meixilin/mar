# MAR: Mutations-Area Relationship Analysis

## Overview
MAR is an R package that enables the reconstruction of Mutations-Area Relationships using spatially distributed genome variation data. This tool helps researchers analyze how genetic mutations accumulate across geographic space within a species.

## Installation

### Development Version

To install the development version with all dependencies:

```R
library(devtools)
install_github("meixilin/mar")
```

Or in bash:

```bash
R CMD INSTALL --preclean --no-multiarch --with-keep.source mar
```

### Troubleshooting

We are working on a stable release. If you encounter issues with installation above, try manually installing these dependencies before installing `mar`:

```R
# not devtools
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SeqArray")
install.packages(c("raster", "sars", "sads", "matrixStats"))
```

Use of `SeqArray` version >= 1.28.0 is **required** (there was a bug in the previous version that impacts plink file importing).

## Minimal working example

A minimal working example is provided in the `tests/example` folder. This is a dummy data and does not contain any real data. To run the example, simply run:

```R
library(mar)
library(raster)
MARPIPELINE(name = "example", workdir = "./", genofile = "genome.tsv", lonlatfile = "lonlat.csv", saveobj = TRUE)
```

## Preparing input data

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
