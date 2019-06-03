#' Function to Fit MCMSeq Models
#'
#' Fits negative binomial generlized linear models and negative binomial generalized linear mixed models to RNA-Seq data using MCMC.
#'
#' @param expr_mat A (G x N) numeric matrix or data frame of transformed RNA-seq counts (e.g. using VST from DESeq2), with genes in rows and samples in columns. G = number of genes.  N = number of samples.
#' @param gene_names An optional character vector of gene names (length G).  If unspecified, row names from the expression matrix will be used.
#' @param form A one-sided linear formula describing both the fixed-effects and random-effects parts of the model using the syntax of the lme4 package
#' @param sample_data Data frame with N rows containing the fixed effects terms included in the fixed_effects formula, as well as any random effects listed in random_intercept.  The rows of the data frame must correspond (and be in the same order as) the columns of the expression matrix.
#' @param REML Should the models be fit with REML or regular ML?
#' @export
#'

lmerSeq.fit <- function(expr_mat=NULL, # matrix of transformed RNA-Seq counts where rows are genes and columns are samples
                       gene_names = NULL, # a vector of gene names (the length of the number of rows in the expression matrix).  If unspecified, rownames from the expression matrix will be used.
                       form = NULL, #formula for fixed effects
                       sample_data = NULL,
                       REML = T
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

  ### Inconsistent information ###

  # fixed_terms <- attributes(terms(form))$term.labels[attributes(terms(form))$order==1]
  #
  # if(sum(fixed_terms %in% colnames(sample_data)) != length(fixed_terms)){
  #   stop(paste0("The following fixed effects terms are missing from the sample_data:",
  #               fixed_terms[!(fixed_terms %in% colnames(sample_data))]))}

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
    gene_names = rownames(expr_mat)}


  ############################################################################################################
  # Begin Analysis
  ############################################################################################################
  # Make sure expr_mat is a matrix
  expr_mat <- as.matrix(expr_mat)

  # Ensure that contrast, gene and fixed effect names are supplied as characters
  gene_names <- as.character(gene_names)
  form_sub <- update(form, expr ~ .)
  ret <- pbapply::pblapply(X = 1:nrow(expr_mat), FUN = function(i){
    dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ])))
    ret_sub <- tryCatch({
      tmp1 <- lmerTest::lmer(formula = form_sub, data = dat_sub, REML = REML)
    }, error = function(e) {
      ret_sub2 <- NA
    })
  })
  ret$gene_names = gene_names
  return(ret)
}
