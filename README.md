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

### Troubleshooting

We are working on a stable release. If you encounter issues with installation above, try manually installing these dependencies before installing `mar`:

```R
install.packages(c("magrittr", "raster", "SeqArray", "dplyr", "sars"))
```

Use of `SeqArray` version >= 1.28.0 is **required** (there was a bug in the previous version that impacts plink file importing).

## Minimal working example

A minimal working example is provided in the `tests/example` folder. This is a dummy data and does not contain any real data. To run the example, simply run:

```R
library(mar)
library(raster)
MARPIPELINE(name = "example", workdir = "./", genofile = "genome.tsv", lonlatfile = "lonlat.csv", saveobj = TRUE)
```
