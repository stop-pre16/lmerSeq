# lmerSeq: An R package for fitting linear mixed models to RNA-Seq data
This is the repo for the lmerSeq R package, which has the wrapper functions to run analyses of transformed RNA-Seq counts from experiments with repeated measures or other designs where correlation between samples/observations needs to be accounted for by using linear mixed models.

## Installation
The lmerSeq package can easily be installed from github using devtools. Dependencies can be installed from CRAN. 
- First, the pbapply, lmerTest, dplyr, lme4, magrittr, gtools, nlme, numDeriv, devtools, and lavaSearch2 packages should be installed from CRAN:
```
install.packages(pkgs = c('pbapply', 'lmerTest', 'dplyr', 'devtools', 'magrittr', 'gtools', 'nlme', 'numDeriv', 'lavaSearch2'))
```
- Next, the lmerSeq package can be installed using devtools:
```
devtools::install_github("stop-pre16/lmerSeq", build_vignettes = T)
```

## Tutorials
Tutorials are available in the Tutorials directory of the repository. The following tutorials are available:
- [Detailed tutorial for lmerSeq analysis of a simulated dataset with detailed explanation of each function](https://htmlpreview.github.io/?https://github.com/stop-pre16/lmerSeq/blob/master/Tutorials/lmerSeq_vignette.html)
- [Simple analysis of a publicly available dataset using random effects](https://htmlpreview.github.io/?https://github.com/stop-pre16/lmerSeq/blob/master/Tutorials/CS_shock_tutorial.html)
- [Generalized least squares models and using compound symmetric and unstructured correlation structures](https://htmlpreview.github.io/?https://github.com/stop-pre16/lmerSeq/blob/master/Tutorials/lmerSeq_GLS_vignette.html)

In addition, a detailed vignette with an example analysis, can be viewed in R:
```
vignette("lmerseq-vignette")
```
