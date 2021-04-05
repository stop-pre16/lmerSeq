#' Function to summarize individual linear contrasts coefficients
#'
#' Conducts t- or F-tests on linear contrasts of regression coefficients from fits done using lmerSeq_fit function.  If the contrast matrix has only 1 row, a t-test is done.  If the contrast matrix has more than 1 row, an F-test is done
#'
#' @param lmerSeq_results Results object from running lmerSeq.fit
#' @param contrast_mat Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See p.adjust.methods
#' @param ddf Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
#' @param sort_results Should the results table be sorted by adjusted p-value?
#'
#' @examples
#' data("expr_data")
#' vst_expr <- expr_example$vst_expr
#' sample_meta_data <- expr_example$sample_meta_data
#'
#' ##Only including 10 genes in the expression matrix
#' vst_expr <- vst_expr[1:10, ]
#'
#' ##  Fit the Model
#' fit_lmerSeq <- lmerSeq.fit(form = ~ group * time + (1|ids),
#'                            expr_mat = vst_expr,
#'                            sample_data = sample_meta_data,
#'                            REML = T)
#'
#'
#' ##  1 dimensional contrast (t-test)
#' cont_mat1 <- rbind(c(0, 1, 0, 1)) # group diff. at followup
#' contrast_summary1 <- lmerSeq.contrast(lmerSeq_results = fit_lmerSeq,
#'                                       contrast_mat = cont_mat1,
#'                                       p_adj_method = 'BH',
#'                                       ddf = 'Satterthwaite',
#'                                       sort_results = T)
#' print(contrast_summary1)
#'
#' ##  multi-dimensional contrast (F-test)
#' cont_mat2 <- rbind(c(0, 1, 0, 0),
#'                    c(0, 0, 1, 0),
#'                    c(0, 0, 0, 1)) # simultaneous test of all coefficients
#' contrast_summary2 <- lmerSeq.contrast(lmerSeq_results = fit_lmerSeq,
#'                                       contrast_mat = cont_mat2,
#'                                       p_adj_method = 'BH',
#'                                       ddf = 'Satterthwaite',
#'                                       sort_results = T)
#' print(contrast_summary2)
#'
#' @export
#'

lmerSeq.contrast <- function(lmerSeq_results = NULL, # Results object from running lmerSeq.fit
                            contrast_mat = NULL, # Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
                            p_adj_method = "BH", # Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                            ddf = "Satterthwaite", # Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
                            sort_results = T, # Should the results table be sorted by adjusted p-value?
                            include_singular = F # Should genes with singular fits be included in the results?
){
  n_fits = length(lmerSeq_results)
  idx_tmp = 1
  while(is.null(lmerSeq_results[[idx_tmp]]$fit) & idx_tmp <= n_fits){
    idx_tmp = idx_tmp + 1
  }
  if(idx_tmp > n_fits){
    stop("Model fits for all genes are null")
  }
  coef_names <- names(lme4::fixef(lmerSeq_results[[idx_tmp]]$fit))
  gene_names <- do.call(c, lapply(lmerSeq_results, function(x){return(x$gene)}))
  colnames(contrast_mat) <- coef_names
  ############################################################################################################
  #Error Messages for insufficient or inconsistent information
  ############################################################################################################
  if(!(p_adj_method %in% p.adjust.methods)){
    stop("Invalid p_adj_method")
  }

  if(!(ddf %in% c("Satterthwaite", "Kenward-Roger"))){
    stop("Invalid ddf method")
  }

  if(ncol(contrast_mat) != length(coef_names)){
    stop("Number of columns in the contrast matrix does not match number of model coefficients")
  }

  joint_flag <- nrow(contrast_mat) > 1

  # idx_singular <- do.call(c, lapply(lmerSeq_results, function(x){isSingular(x$fit)}))
  idx_singular <- do.call(c, lapply(lmerSeq_results, function(x){
    if(is.null(x$fit)){
      return(T)
    }
    else{
      return(lme4::isSingular(x$fit))
    }
  }))
  genes_singular_fits <- gene_names[idx_singular]
  ret <- do.call(rbind, pblapply(lmerSeq_results, function(x){
    # x = lmerSeq_results[[1]]
    if(is.null(x)){
      return(NA)
    }
    cont_sub <- lmerTest::contest(x$fit, L = contrast_mat, joint = joint_flag, ddf = ddf)
    return(cont_sub)
  }))
  if(include_singular == F){
    ret[idx_singular, ] = NA
  }
  ret <- data.frame(gene = gene_names, ret)
  if(joint_flag){
    names(ret)[7] <- 'p_val_raw'
  }
  else{
    names(ret)[8] <- 'p_val_raw'
  }
  ret <- ret %>% mutate(p_val_adj = p.adjust(p_val_raw, method = p_adj_method))
  if(sort_results){
    ret <- ret %>%
      arrange(p_val_adj)
  }
  # ret <- gtools::smartbind(ret, data.frame(gene = genes_singular_fits))
  rownames(ret) <- NULL
  ret2 <- list(contrast_mat = contrast_mat,
               summary_table = ret,
               singular_fits = genes_singular_fits,
               ddf = ddf,
               p_adj_method = p_adj_method)
  return(ret2)
}
