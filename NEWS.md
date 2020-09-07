# Rdimtools 1.0.4

* Methods for feature selection `do.wdfs`, `do.uwdfs`, `do.procrustes`, and `do.mifs` are added.
* `README` now shows numbers of currently available functions for DR and IDE.

# Rdimtools 1.0.3

* Fixed memory leaks in `do.sne` and `do.tsne`.
* `do.lsls` added as a supervised feature selection method.

# Rdimtools 1.0.2

* README contains minimal examples for both dimension reduction and estimation.
* Porting to pure C++ implementations started, gaining computational efficiency. 
* `do.lmds` function is fixed for its discrepancy in nested structure.


# Rdimtools 1.0.1

* NEWS reformatted and [package website](http://kyoustat.com/Rdimtools/) is now available.
* Dependency on R package [ADMM](https://CRAN.R-project.org/package=ADMM) is removed.


# Rdimtools 1.0.0

## Major changes
* LDA solves trace ratio problem directly.
* Many of dependencies are removed.
* 133 dimension reduction methods available.
* 17 intrinsic dimension estimation methods available.

## Bug fixes
* Error fixed in `do.lscore` function (thanks to Jordan Lin).

