# v0.0.1

Same code as Moi's original release.

Edits are:

* DESCRIPTION: Removed one line so library can be built.
* .gitignore
* Removed .DS_Store

# v0.0.2

Formal refactoring starts. Although code changes, it creates the same output as Moi's original release.

Edits are:

* Massive code restyling (2 space to 4 space; = to <- etc)
* Create individual files for mutdiv and MARextinction
* Add tests for package
* Remove implementation for `rasterN`
* Remove bug fix for `MARextinction_radial` in commit f1696

Test output:

* `test-MARextinction.R`: all PASS.

## 24f0 and 6289

Edits are:

* Changed `cpp` implementation for `Hn` and used R's native implementation for harmonic number.

Test output:

* `test-MARextinction.R`: all not PASS without tolerance = 1e-07. (commit id: 24f0)
* `test-MARextinction.R`: all PASS with tolerance = 1e-07. (commit id: 6289)

Cause of change:

small numeric variations.

## 2a40

Edits are:

* Put back the bug fix for `MARextinction_radial` in commit f1696

Test output:

* `test-MARextinction.R`: three PASS with tolerance = 1e-07 but `MARextinction_radial` test does not PASS.

previous **xsim-joshua-radial.rds** now named **xsim-joshua-radial.6289.rds** and becomes obselete in this commit.

new **xsim-joshua-radial.rds** generated

# v0.0.3

Put in small fix on P calculation within the genemaps framework. Used for Kristy's extinction-sim outputs for now.

Edits are:

* Count multiple samples in one raster cell. `P ... sum(cells > 0)` to `P ... sum(cells)`
* Fix harmonic numbers from `Hn(N)` to `Hn(N-1)`
* Add `n/(n-1)` term for `theta_pi` calculation
* Changed conditions for `N > 0` to `N > 1`

Test output:

* `test-MARextinction.R`: expected to pass when considering only the columns 3 to 9. The differences are plotted when running the test.
* `test-mutdiv.R`: newly added and expected to pass.

# v0.0.4

* IMPORTANT: Removed old genemap as the central data structure. Introduces `genemaps` as a S3 class.
    * `.lonlat_res` function added to automatically select the optimal resolution in generating `genemaps`
    * `.raster_lonlatr` function generates `samplemap`.
* Decided to completely discard the `raster` usage in `raster_samples` and `raster_mutmaps`.
    * `genemaps` now have a slot `samplemap` that is the same as `raster_samples`. It is used to align operations in the same grid systems.
* Fix `create_gene_maps` bug that it misses some samples at the edge. See `bugs-creategeno.R` and `mar.R`
    * Now `create_gene_maps` is obselete. Use `genemaps` constructor.
* Add a debug line in `MARsampling` to track the cells circled in each replicate runs. See `mar.R`
* Generate test file `genemaps_new-joshua.rda` from fixed `create_gene_maps` function.
* Generate test file `mares_new-joshua.rda` from `MARsampling` random scheme.
* Other:
    * added `.match` function (but not used)
    * speed profiler `profilespeed` (added but not tracked in git)

Test output:

* `test-genemaps.R`: creates `gm-joshua.rda` and `gm-arabidopsis.rda` objects. validates that samples are mapped correctedly without loss in data.

TODO:

* Sampling can be extended outside of the gridded system used by `raster_samples` and `samplemap`. But it is impossible to test for reproducibility.

# v0.0.5

* IMPORTANT: Greatly speed up `MARsampling` while ensuring reproducibility.
    * `MARsampling` in `mar.R` renamed to `MARsampling_old`
* Added two supporting function groups `mutdiv.R` and `opsraster.R` to calculate genetic diversity and perform operations on raster files.

Test output:
* `test-mutdiv`: test that new `.mutdiv.gridded` function works to reproduce results from `mutdiv` in previous version in full genemaps samples.
* `test-opsraster`: test that new `areaofraster` function works to reporduce results from `areaofraster_old`
* `test-MARsampling`: completely reproduces the results from `MARsampling_old`. The `Asq` is not matched with `a` in the last 10 replicates. Because there was a small bug in `MARsampling_old` that creates 0 as a starting point.
* `test-MARextinction`: now rolled out, awaiting new updates.
* Ran `usethis::use_pipe` and `devtools::document` to import the pipe `%>%` function.

# v0.0.6

* Removed `mutdiv_old` `areaofraster_old` and `MARsampling_old`
* Fix `MARsampling` bugs that 1. creates 0 as a starting point and 2. cannot start with only one cell selected.
* Fix `lonrange` bugs that should be `latrange`. See `xFromCol` and `yFromRow`.
* Add plotting methods and CRS for `genemaps`.
* Decided to not use `leaflet` as the main plotting method, but available in scratch.
* Add consistent `MARsampling` sampling method support modulated by the `prob` parameter in sample. (Available for extinction as well).
* Add `MARcalc`.
* All tests updated.

# v0.0.7

* Folder cleanup
* Add documentations for major functions.
* Removed Cpp codes as these are not needed to run `MARPIPELINE`. And SFS generation can be fast in R.

# v0.0.8

* Add file parser and use `SeqArray` to read genotype files.
* Change `.mutdiv.gridded` to `mutdiv.gridded`
* Update `genemaps` to `genomaps` and support `SeqArray` `margeno` format changes.
* Discontinued `padding` option in `.lonlat_raster` to avoid rounding errors.
* Update `MARPIPELINE` to use the new `genomaps` format.

# v0.0.9

* MARPIPELINE completely rewritten.
* Add `mutdiv.cells` to calculate genetic diversity in a cell list.
* Rewrite `MARextinction` to match `MARsampling`.
* Add `MARcalc_all` to calculate MAR/EMAR relationship.
* Update dependencies from `genemaps` object to `genomaps` object `gm$samplemap` to `gm$maps$samplemap` etc.
* Not fully tested yet.

# Did not use but tried

* Considered using `SeqArray` as the main framework for `genomaps` but it would require reading from disk constantly.
    * Now only use `SeqArray` to read plink and vcf files.

# TODO

* Add `sfs` functionality for genotypes
* Add missing data handling for `mutdiv`
* Update `MARextinction` to match `MARsampling`
* Update `MARPIPELINE`

# Notes

* Heterozygous genotypes can be handled normally given the definition of `pi`.
* Non-biallelic SNPs are not implemented.
* Theoretically, this pipeline works for any ploidy but the interpretation of the results is not straightforward for non-diploid organisms.
* Potentials to extrapolate genetic diversity in area of know occurrences.

