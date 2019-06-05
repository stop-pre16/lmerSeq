#' Function to Fit MCMSeq Models
#'
#' Fits negative binomial generlized linear models and negative binomial generalized linear mixed models to RNA-Seq data using MCMC.
#'
#' @param lmerSeq_results Results object from running lmerSeq.fit
#' @param contrast_mat Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See p.adjust.methods
#' @param ddf Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
#' @param sort_results Should the results table be sorted by adjusted p-value?
#' @export
#'

lmerSeq.contrast <- function(lmerSeq_results = NULL, # Results object from running lmerSeq.fit
                            contrast_mat = NULL, # Numeric matrix representing the contrast to be tested.  Matrix must have the same number of columns as the number of coefficients in the model.  If the matrix has multiple rows, a simultaneous F-test will be done
                            p_adj_method = "BH", #Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                            ddf = "Satterthwaite", #Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
                            sort_results = T #Should the results table be sorted by adjusted p-value?
){
  gene_names <- do.call(c, lapply(lmerSeq_results, function(x){return(x$gene)}))
  coef_names <- names(fixef(lmerSeq_results[[1]]$fit))
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

  idx_singular <- do.call(c, lapply(lmerSeq_results, function(x){isSingular(x$fit)}))
  genes_singular_fits <- gene_names[idx_singular]
  ret <- do.call(rbind, lapply(lmerSeq_results[!idx_singular], function(x){
    # x = lmerSeq_results[[1]]
    cont_sub <- lmerTest::contest(x$fit, L = contrast_mat, joint = joint_flag, ddf = ddf)
    return(cont_sub)
  }))
  ret <- data.frame(gene = gene_names[!idx_singular], ret)
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
  ret <- gtools::smartbind(ret, data.frame(gene = genes_singular_fits))
  rownames(ret) <- NULL
  ret2 <- list(contrast_mat = contrast_mat,
               summary_table = ret,
               singular_fits = genes_singular_fits,
               ddf = ddf,
               p_adj_method = p_adj_method)
  return(ret2)
}
