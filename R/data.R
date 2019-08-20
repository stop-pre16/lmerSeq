#' Example transformed RNA-Seq data from paired observations
#'
#' A list that has 2 components: A matrix with transformed RNA-Seq values
#' (simulated from a negative binomial and then transformed using the VST in the DESeq2 package) and
#' a data frame that has the sample meta data
#'
#' @format A list with 2 named objects:
#' \describe{
#'   \item{vst_expr}{VST transformed expression}
#'   \item{sample_meta_data}{data frame with variables for group, time, and subject ID.  Rows are in the same order as columns of vst_expr}
#'   ...
#' }
"expr_example"
