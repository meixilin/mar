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

##

Edits are:

* Changed `cpp` implementation for `Hn` and used R's native implementation for harmonic number.

Test output:

* `test-MARextinction.R`: all not PASS without tolerance = 1e-07.
* `test-MARextinction.R`: all PASS with tolerance = 1e-07.

Cause of change:

small numeric variations.



