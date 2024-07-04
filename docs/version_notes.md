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

# v0.1.0

* IMPORTANT: Removed old genemap as the central data structure. Introduces `genemap` as a S3 class.
* Fix `create_gene_maps` bug that it misses some samples at the edge. See `bugs-creategeno.R`

Test output:


