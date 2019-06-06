#' Function to summarize individual regression coefficients
#'
#' Conducts t-tests on individual regression coefficients from fits done using lmerSeq_fit function
#'
#' @param lmerSeq_results Results object from running lmerSeq.fit
#' @param coefficient Character string or numeric indicator of which coefficient to summarize
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See p.adjust.methods
#' @param ddf Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
#' @param sort_results Should the results table be sorted by adjusted p-value?
#' @export
#'

lmerSeq.summary <- function(lmerSeq_results = NULL, # Results object from running lmerSeq.fit
                        coefficient = NULL, # Character string or numeric indicator of which coefficient to summarize
                        p_adj_method = "BH", # Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                        ddf = "Satterthwaite", # Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
                        sort_results = T # Should the results table be sorted by adjusted p-value?
){
  coef_names <- names(lme4::fixef(lmerSeq_results[[1]]$fit))
  gene_names <- do.call(c, lapply(lmerSeq_results, function(x){return(x$gene)}))
  ############################################################################################################
  #Error Messages for insufficient or inconsistent information
  ############################################################################################################
  if(is.numeric(coefficient)){
    if((coefficient > length(coef_names)) | coefficient < 1){
      stop("Coefficient number is invalid")
    }
    coef_out <- coef_names[coefficient]
  }
  else{
    if(!(coefficient %in% coef_names)){
      stop("Coefficient name is invalid")
    }
    coef_out <- coefficient
  }

  if(!(p_adj_method %in% p.adjust.methods)){
    stop("Invalid p_adj_method")
  }

  if(!(ddf %in% c("Satterthwaite", "Kenward-Roger"))){
    stop("Invalid ddf method")
  }

  idx_singular <- do.call(c, lapply(lmerSeq_results, function(x){lme4::isSingular(x$fit)}))
  genes_singular_fits <- gene_names[idx_singular]
  ret <- do.call(rbind, lapply(lmerSeq_results[!idx_singular], function(x){
    # x = lmerSeq_results$fitted_models[[1]]
    res_sub <- summary(x$fit, ddf = ddf)$coefficients[coefficient, ]
    return(res_sub)
  }))
  ret <- data.frame(gene = gene_names[!idx_singular], ret)
  names(ret)[6] <- 'p_val_raw'
  ret <- ret %>% mutate(p_val_adj = p.adjust(p_val_raw, method = p_adj_method))
  if(sort_results){
    ret <- ret %>%
      arrange(p_val_adj)
  }
  ret <- gtools::smartbind(ret, data.frame(gene = genes_singular_fits))
  rownames(ret) <- NULL
  ret2 <- list(coefficient = coef_out,
               summary_table = ret,
               singular_fits = genes_singular_fits,
               ddf = ddf,
               p_adj_method = p_adj_method)
  return(ret2)
}
