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

To install the development version with all dependencies and the vignettes:

```R
library(devtools)
devtools::install_github("meixilin/mar", build_vignettes = TRUE)
```

Or in bash:

```bash
R CMD build mar && R CMD INSTALL --preclean --no-multiarch --with-keep.source mar_0.3.0.tar.gz
```

CRAN installation instruction upcoming!

### Troubleshooting

If you encounter issues with installation above, try manually installing these dependencies before installing `mar`:

```R
install.packages(c("sars", "matrixStats", "terra"))
```

Or using `conda` / `mamba`:

```bash
mamba create --yes --name mar r-base r-devtools zlib libzlib gdal r-terra r-curl
mamba activate mar
R
library(devtools)
devtools::install_github("meixilin/mar")
```

If the above tricks do not work, please open an issue on GitHub.

## Minimal working example

```R
library(mar)
vignette("mar-workflow", package = "mar")
```

Three vignettes ship with the package:

| Vignette | Topic |
| --- | --- |
| `vignette("mar-workflow", package = "mar")` | MAR Analysis Workflow |
| `vignette("mar-theory", package = "mar")` | The sampling theory behind MAR |
| `vignette("mar-paper", package = "mar")` | Worked example: global Arabidopsis 1001G dataset |
