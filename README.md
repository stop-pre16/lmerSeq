# lmerSeq
This is the repo for the lmerSeq R package, which has the wrapper functions to run analyses of transformed RNA-Seq counts from experiments with repeated measures or other designs where correlation between samples/observations needs to be accounted for by using linear mixed models.

To install the package, the pbapply, lmerTest, dplyr, lme4, magrittr, gtools, nlme, numDeriv, devtools, and lavaSearch2 packages need to first be installed:

install.packages(pkgs = c('pbapply', 'lmerTest', 'dplyr', 'devtools', 'magrittr', 'gtools', 'nlme', 'numDeriv', 'lavaSearch2'))

Next, the mcmseq package can be installed using devtools:

devtools::install_github("stop-pre16/lmerSeq", build_vignettes = T)

To read the detailed vignette with an example analysis, run:

vignette("lmerSeq-vignette")
