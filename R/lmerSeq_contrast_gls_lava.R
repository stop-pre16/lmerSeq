#' Function to summarize individual linear contrasts coefficients
#'
#' Conducts t- or F-tests on linear contrasts of regression coefficients from fits done using lmerSeq.fit.gls function.  If the contrast matrix has only 1 row, a t-test is done.  If the contrast matrix has more than 1 row, an F-test is done
#'
#' @param lmerSeq_results Results object from running lmerSeq.fit.gls
#' @param contrast_mat Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See p.adjust.methods
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
#' fit.lmerSeq.gls <- lmerSeq.fit.gls(form = ~ group * time,
#'                            cor_str = corCompSymm(form = ~ 1 | id),
#'                            expr_mat = vst_expr,
#'                            sample_data = sample_meta_data,
#'                            method = 'REML')
#'
#'
#' ##  1 dimensional contrast (t-test)
#' cont_mat1 <- rbind(c(0, 1, 0, 1)) # group diff. at followup
#' contrast_summary1 <- lmerSeq.contrast.gls.lava(lmerSeq_results = fit.lmerSeq.gls,
#'                                       contrast_mat = cont_mat1,
#'                                       p_adj_method = 'BH',
#'                                       sort_results = T)
#' print(contrast_summary1)
#'
#' ##  multi-dimensional contrast (F-test)
#' cont_mat2 <- rbind(c(0, 1, 0, 0),
#'                    c(0, 0, 1, 0),
#'                    c(0, 0, 0, 1)) # simultaneous test of all coefficients
#' contrast_summary2 <- lmerSeq.contrast.gls.lava(lmerSeq_results = fit.lmerSeq.gls,
#'                                       contrast_mat = cont_mat2,
#'                                       p_adj_method = 'BH',
#'                                       sort_results = T)
#' print(contrast_summary2)
#'
#' @export
#'

lmerSeq.contrast.gls.lava <- function(lmerSeq_results = NULL, # Results object from running lmerSeq.fit.gls
                                      contrast_mat = NULL, # Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
                                      p_adj_method = "BH", # Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                                      sort_results = T # Should the results table be sorted by adjusted p-value?
){
  n_fits = length(lmerSeq_results)
  idx_tmp = 1
  while(is.null(lmerSeq_results[[idx_tmp]]$fit) & idx_tmp <= n_fits){
    idx_tmp = idx_tmp + 1
  }
  if(idx_tmp > n_fits){
    stop("Model fits for all genes are null")
  }
  coef_names <- names(lmerSeq_results[[idx_tmp]]$fit$coefficients)
  gene_names <- do.call(c, lapply(lmerSeq_results, function(x){return(x$gene)}))
  colnames(contrast_mat) <- coef_names
  form_sub = lmerSeq_results[[idx_tmp]]$formula
  ############################################################################################################
  #Error Messages for insufficient or inconsistent information
  ############################################################################################################
  if(!(p_adj_method %in% p.adjust.methods)){
    stop("Invalid p_adj_method")
  }

  # if(!(ddf %in% c("Satterthwaite", "Kenward-Roger"))){
  #   stop("Invalid ddf method")
  # }

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
      # return(is.character(x$fit$apVar) | abs(getVarCov(x$fit)[1, 2]) < 1e-5)
      return(is.character(x$fit$apVar))
    }
  }))
  genes_singular_fits <- gene_names[idx_singular]
  ret <- do.call(rbind, pblapply(lmerSeq_results, function(x){
    # x = lmerSeq_results[[1]]
    if(is.null(x$fit)){
      return(NA)
    }
    # else if(is.character(x$fit$apVar) | abs(getVarCov(x$fit)[1, 2]) < 1e-5){
    else if(is.character(x$fit$apVar)){
      return(NA)
    }
    else{
      cont_sub = satt_test_gls_lava(L = contrast_mat,
                                    gls.mod = x$fit,
                                    data = x$data,
                                    joint = joint_flag)
      return(cont_sub)
    }
  }))
  ret[idx_singular, ] = NA
  ret <- data.frame(gene = gene_names, ret)
  if(joint_flag){
    names(ret)[5] <- 'p_val_raw'
  }
  else{
    names(ret)[6] <- 'p_val_raw'
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
               # ddf = ddf,
               p_adj_method = p_adj_method)
  return(ret2)
}
