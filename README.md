# VCM
VCM package contains three approaches for solving variance components model: PX-EM algorithm, MM algorithm and Method of Moments

Installation
===========

To install the development version of bivas, it's easiest to use the 'devtools' package. Note that bivas depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mxcai/bivas")
```

Usage
===========
[The 'bivas' vignette](https://github.com/mxcai/bivas/blob/master/vignettes/bivas_package.pdf?raw=true) provides a quick start for the usage of the package. The following help page also provides quick reference and examples:

```
library(bivas)
package?bivas
```

Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-bivas](https://github.com/mxcai/sim-bivas).

References
==========

Mingxuan Cai, Mingwei Dai, Jingsi Ming, Heng Peng, Jin Liu and Can Yang. BIVAS: A scalable Bayesian method for bi-level variable selection. Under review. arXiv link: http://arxiv.org/abs/1803.10439.


Development
==========

This R package is developed by Mingxuan Cai and Can Yang (macyang@ust.hk)
