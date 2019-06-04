#' Function to Fit MCMSeq Models
#'
#' Fits negative binomial generlized linear models and negative binomial generalized linear mixed models to RNA-Seq data using MCMC.
#'
#' @param lmerSeq_results Results object from running lmerSeq.fit
#' @param coefficient Character string or numeric indicator of which coefficient or contrast to summarize
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See p.adjust.methods
#' @param ddf Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
#' @export
#'

lmerSeq.summary <- function(lmerSeq_results = NULL, # Results object from running lmerSeq.fit
                        coefficient = NULL, # Character string or numeric indicator of which coefficient or contrast to summarize
                        p_adj_method = "BH", #Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                        ddf = "Satterthwaite" #Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
){
  coef_names <- names(fixef(lmerSeq_results$fitted_models[[1]]))
  gene_names <- lmerSeq_results$gene_names
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

  idx_singular <- do.call(c, lapply(lmerSeq_results$fitted_models, isSingular))
  genes_singular_fits <- gene_names[idx_singular]
  ret <- do.call(rbind, lapply(lmerSeq_results$fitted_models[!idx_singular], function(x){
    # x = lmerSeq_results$fitted_models[[1]]
    res_sub <- summary(x, ddf = ddf)$coefficients[coefficient, ]
    return(res_sub)
  }))
  ret <- data.frame(gene = gene_names[!idx_singular], ret)
  names(ret)[6] <- 'p_val_raw'
  ret <- ret %>% mutate(p_val_adj = p.adjust(p_val_raw, method = p_adj_method)) %>%
    arrange(p_val_adj)
  ret2 <- list(coefficient = coef_out, summary_table = ret, singular_fits = genes_singular_fits)
  return(ret2)
}
