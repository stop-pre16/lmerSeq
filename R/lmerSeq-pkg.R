#' @title lmerSeq: Wrappers for Fitting Linear Mixed Models to RNA-Seq Data
#'
#' @description
#' Fitting utilities for estimating linear mixed models for use
#' with RNA-Seq data from studies with clustered/longitudinal study designs
#'
#' @importFrom lme4 fixef
#' @importFrom lme4 isSingular
#' @importFrom magrittr %>%
#' @importFrom lmerTest lmer
#' @importFrom lmerTest contest
#' @importFrom pbapply pblapply
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom gtools smartbind
#' @importFrom parallel mclapply
#' @importFrom nlme gls
#' @importFrom nlme getVarCov
#' @importFrom nlme glsEstimate
#' @importFrom numDeriv jacobian
#' @importFrom lavaSearch2 sCorrect
#' @importFrom lavaSearch2 compare2
#'
#'
#' @docType package
#' @name lmerSeq-pkg
NULL
