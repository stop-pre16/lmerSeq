#' Function to Fit linear mixed models to transformed RNA-Seq data
#'
#' Wrapper function that its linear mixed models to (transformed) RNA-Seq data using the utilities present using the gls function from the nlme package.
#'
#' @param form A one-sided linear formula describing both the fixed-effects part of the model using the syntax of the nlme package
#' @param cor_str Correlation structure defined as one of the corStruct options available in the nlme package
#' @param expr_mat A (G x N) numeric matrix or data frame of transformed RNA-seq counts (e.g. using VST from DESeq2), with genes in rows and samples in columns. G = number of genes.  N = number of samples.
#' @param gene_names An optional character vector of gene names (length G).  If unspecified, row names from the expression matrix will be used.
#' @param sample_data Data frame with N rows containing the fixed- and random-effects terms included in the formula.  The rows of the data frame must correspond (and be in the same order as) the columns of the expression matrix.
#' @param method Should the models be fit with REML or regular ML?
#' @param parallel If on Mac or linux, use forking (via mclapply) to parallelize fits
#' @param cores Number of cores to use (default is 2)
#'
#' @examples
#' data("expr_data")
#' vst_expr <- expr_example$vst_expr
#' sample_meta_data <- expr_example$sample_meta_data
#'
#' ##  Only including 10 genes in the expression matrix
#' vst_expr <- vst_expr[1:10, ]
#'
#' ##  Fit the Model
#' fit.lmerSeq.gls <- lmerSeq.fit.gls(form = ~ group * time,
#'                            cor_str = corCompSymm(form = ~ 1 | id),
#'                            expr_mat = vst_expr,
#'                            sample_data = sample_meta_data,
#'                            method = 'REML')
#'
#' @export
#'

lmerSeq.fit.gls <- function(form = NULL, # Formula for fixed effects
                            cor_str = NULL, # Correlation structure defined as one of the corStruct options available in the nlme package
                            expr_mat = NULL, # Matrix of transformed RNA-Seq counts where rows are genes and columns are samples
                            gene_names = NULL, # A vector of gene names (the length of the number of rows in the expression matrix).  If unspecified, rownames from the expression matrix will be used.
                            sample_data = NULL, # A data frame with sample meta data
                            weights = NULL, # Matrix of same dimension as expr_mat with gene-specific weights for each sample
                            method = "REML", # Fit mixed models using REML or ML
                            parallel = F,
                            cores = 2
){

  ############################################################################################################
  #Error Messages for insufficient or inconsistent information
  ############################################################################################################

  ### Insufficient Information ###

  if(is.null(expr_mat)==T ) {
    stop("An expression matrix must be provided.")}

  if(is.null(form)==T ) {
    stop("A formula must be provided.")}

  if(is.null(sample_data)==T ) {
    stop("sample_data is missing.")}

  if(is.null(cor_str)==T ) {
    stop("correlation structure is missing.")}

  ### Inconsistent information ###

  if((ncol(expr_mat)==nrow(sample_data))==F ) {
    stop("The expression matrix and sample data include differing numbers of samples.")}

  if(is.null(gene_names)==F & (nrow(expr_mat)==length(gene_names))==F ) {
    print("The expression matrix and gene_names indicate differing numbers of genes.
          Row names of the expression matrix will be used as gene names.")
    gene_names = rownames(expr_mat)
  }

  ################################################################################################
  # Calculate Default Values if none supplied
  ################################################################################################

  # Gene Names
  if(is.null(gene_names)){
    if(is.null(rownames(expr_mat))==T){rownames(expr_mat)<-seq(1,nrow(expr_mat),1)}
    gene_names = rownames(expr_mat)
  }

  ############################################################################################################
  # Begin Analysis
  ############################################################################################################
  # Make sure expr_mat is a matrix
  expr_mat <- as.matrix(expr_mat)

  # Ensure that contrast, gene and fixed effect names are supplied as characters
  gene_names <- as.character(gene_names)
  form_sub <- update(form, expr ~ .)
  n_samples = ncol(expr_mat)
  if(parallel == F){
    if(is.null(weights)){
      ret <- pbapply::pblapply(X = 1:nrow(expr_mat),
                               # ret <- lapply(X = 1:nrow(expr_mat),
                               FUN = function(i){

                                 dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ])))
                                 ret_sub <- tryCatch({
                                   tmp1 <- suppressMessages(nlme::gls(model = form_sub,
                                                                      correlation = cor_str,
                                                                      data = dat_sub,
                                                                      method = method))
                                   suppressWarnings(lavaSearch2::sCorrect(tmp1, df = TRUE, numeric.derivative=TRUE) <- FALSE)
                                   tmp1
                                 }, error = function(e) {
                                   ret_sub2 <- NULL
                                 })
                                 # suppressWarnings(lavaSearch2::sCorrect(ret_sub, df = TRUE, numeric.derivative=TRUE) <- FALSE)
                                 ret2 <- list(fit = ret_sub,
                                              gene = gene_names[i],
                                              formula = form_sub,
                                              data = dat_sub)
                               })
    }
    else{
      ret <- pbapply::pblapply(X = 1:nrow(expr_mat),
                               FUN = function(i){
                                 dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ]),
                                                                          weights_sub = weights[i, ]))

                                 ret_sub <- tryCatch({
                                   tmp1 <- suppressMessages(nlme::gls(model = form_sub,
                                                                      correlation = cor_str,
                                                                      data = dat_sub,
                                                                      REML = REML,
                                                                      weights = weights_sub))
                                 }, error = function(e) {
                                   ret_sub2 <- NULL
                                 })
                                 suppressWarnings(lavaSearch2::sCorrect(ret_sub, df = TRUE, numeric.derivative=TRUE) <- FALSE)
                                 ret2 <- list(fit = ret_sub,
                                              gene = gene_names[i],
                                              formula = form_sub,
                                              data = dat_sub)

                               })
    }
  }
  else{
    if(is.null(weights)){
      ret = parallel::mclapply(1:nrow(expr_mat), mc.silent = T, mc.cores = cores, function(i){

        dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ])))

        ret_sub <- tryCatch({
          tmp1 <- suppressMessages(nlme::gls(model = form_sub,
                                             correlation = cor_str,
                                             data = dat_sub,
                                             REML = REML))
        }, error = function(e) {
          ret_sub2 <- NULL
        })
        suppressWarnings(lavaSearch2::sCorrect(ret_sub, df = TRUE, numeric.derivative=TRUE) <- FALSE)
        ret2 <- list(fit = ret_sub,
                     gene = gene_names[i],
                     formula = form_sub,
                     data = dat_sub)
      })
    }
    else{
      df_combined = cbind(expr_mat, weights)
      ret = parallel::mclapply(1:nrow(expr_mat), mc.silent = T, mc.cores = cores, function(i){
        dat_sub <- cbind(sample_data,
                         data.frame(expr = as.numeric(df_combined[i, 1:n_samples]),
                                    weights_sub = as.numeric(df_combined[i, (n_samples+1):(2*n_samples)])))

        ret_sub <- tryCatch({
          tmp1 <- suppressMessages(nlme::gls(model = form_sub,
                                             correlation = cor_str,
                                             data = dat_sub,
                                             REML = REML,
                                             weights = weights_sub))
        }, error = function(e) {
          ret_sub2 <- NULL
        })
        suppressWarnings(lavaSearch2::sCorrect(ret_sub, df = TRUE, numeric.derivative=TRUE) <- FALSE)
        ret2 <- list(fit = ret_sub,
                     gene = gene_names[i],
                     formula = form_sub,
                     data = dat_sub)

      })
    }
  }
  # names(ret) = gene_names
  return(ret)
}
