VCM
===========
VCM package contains three approaches for solving variance components model: PX-EM algorithm, MM algorithm and Method of Moments.

Installation
===========

To install VCM, use the following code in R console:

```
#install.packages("devtools")
library(devtools)
install_github("mxcai/VCM")
```

Usage
===========
[The 'VCM' vignette](https://github.com/mxcai/VCM/blob/master/inst/doc/vignette.pdf) provides details of derivation and a quick start for the usage of the package.


References
==========

Jiang, J., Li, C., Paul, D., Yang, C., & Zhao, H. (2016). On high-dimensional misspecified mixed model analysis in genome-wide association study. The Annals of Statistics, 44(5), 2127-2160.

Liu, C., Rubin, D. B., & Wu, Y. N. (1998). Parameter expansion to accelerate EM: the PX-EM algorithm. Biometrika, 85(4), 755-770.

Foulley, J. L., & Van Dyk, D. A. (2000). The PX-EM algorithm for fast stable fitting of Henderson's mixed model. Genetics Selection Evolution, 32(2), 143.

Zhou, H., Hu, L., Zhou, J., & Lange, K. (2018). MM algorithms for variance components models. Journal of Computational and Graphical Statistics, (just-accepted), 1-30.

Wu, Y., & Sankararaman, S. (2018). A scalable estimator of SNP heritability for biobank-scale data. Bioinformatics, 34(13), i187-i194.


Development
==========

This R package is developed by Mingxuan Cai (mcaiad@ust.hk).
