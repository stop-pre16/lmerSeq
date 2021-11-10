
#########################################################
# Helper functions
#########################################################

# Function to calculate quad forms (variance of contrast)
qform <- function(x, A) {
  sum(x * (A %*% x)) # quadratic form: x'Ax
}

# From EMMEANS PACKAGE
.get.lt = function(X) {
  rtn = X[lower.tri(X, diag = TRUE)]
  attr(rtn, "nrow") = nrow(X)
  rtn
}

# From EMMEANS PACKAGE
.lt2mat = function(lt) {
  if (is.null(n <- attr(lt, "nrow")))
    n = (sqrt(8 * length(lt) + 1) - 1) / 2
  X = matrix(NA, nrow = n, ncol = n)
  lti = which(lower.tri(X, diag = TRUE))
  X[lti] = lt
  X = t(X)
  X[lti] = lt
  X
}

# Gradient From EMMEANS PACKAGE
gls_grad = function(object, call, data, V) {
  obj = object$modelStruct
  conLin = object
  class(conLin) = class(obj)
  X = model.matrix(eval(call$model), data = data)
  y = data[[all.vars(call)[1]]]
  conLin$Xy = cbind(X, y)
  conLin$fixedSigma = FALSE
  func = function(x) {
    obj = nlme::`coef<-`(obj, value = x)
    tmp = nlme::glsEstimate(obj, conLin)
    .get.lt(crossprod(tmp$sigma * tmp$varBeta))  # lower triangular form
  }
  res = numDeriv::jacobian(func, coef(obj))
  G = lapply(seq_len(ncol(res)), function(j) .lt2mat(res[, j]))
  G[[1 + length(G)]] = 2 * V  # gradient wrt log sigma
  G
}

gls_grad_lmerSeq = function(object, call, data, V, form_fit) {
  obj = object$modelStruct
  conLin = object
  class(conLin) = class(obj)
  X = model.matrix(form_fit, data = data)
  y = data$expr
  conLin$Xy = cbind(X, y)
  conLin$fixedSigma = FALSE
  func = function(x) {
    obj = nlme::`coef<-`(obj, value = x)
    tmp = nlme::glsEstimate(obj, conLin)
    .get.lt(crossprod(tmp$sigma * tmp$varBeta))  # lower triangular form
  }
  res = numDeriv::jacobian(func, coef(obj))
  G = lapply(seq_len(ncol(res)), function(j) .lt2mat(res[, j]))
  G[[1 + length(G)]] = 2 * V  # gradient wrt log sigma
  G
}

# Output T Table // Adapted from LMER TEST
mk_ttable_gls <- function(estimate, se, ddf) {
  tstat <- (estimate)/se
  pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
  data.frame("Estimate"=estimate, "Std. Error"=se, "df"=ddf,
             "t value"=tstat, "Pr(>|t|)"=pvalue, check.names=FALSE)
}

# Output F Table // Adapted from LMER TEST
mk_Ftable_gls <- function(Fvalue, ndf, ddf) {
  pvalue <- pf(q=Fvalue, df1=ndf, df2=ddf, lower.tail=FALSE)
  data.frame("NumDF"=ndf, "DenDF"=ddf,
             "F value"=Fvalue, "Pr(>F)"=pvalue, check.names = FALSE)
}

# Fstat ddf from T stat ddf // LMER TEST
get_Fstat_ddf_gls <- function(nu, tol=1e-8) {
  # Computes denominator df for an F-statistic that is derived from a sum of
  # squared t-statistics each with nu_m degrees of freedom.
  #
  # nu : vector of denominator df for the t-statistics
  # tol: tolerance on the consequtive differences between elements of nu to
  #      determine if mean(nu) should be returned.
  #
  # Result: a numeric scalar
  #
  # Returns nu if length(nu) == 1. Returns mean(nu) if all(abs(diff(nu)) < tol;
  # otherwise ddf appears to be downward biased.
  fun <- function(nu) {
    if(any(nu <= 2)) 2 else {
      E <- sum(nu / (nu - 2))
      2 * E / (E - (length(nu))) # q = length(nu) : number of t-statistics
    }
  }
  stopifnot(length(nu) >= 1,
            # all(nu > 0), # returns 2 if any(nu < 2)
            all(sapply(nu, is.numeric)))
  if(length(nu) == 1L) return(nu)
  if(all(abs(diff(nu)) < tol)) return(mean(nu))
  if(!is.list(nu)) fun(nu) else vapply(nu, fun, numeric(1L))
}


#########################################################
# Function for GLS Satterthwaite F and T tests
#########################################################

# L = contrast matrix to test
# gls.mod = gls model fit
# data = data frame used to fit the model
# joint = TRUE for F test, FALSE for individuval T tests

# Single DF test
satt_gls_1df <- function(L, gls.mod, G, joint){
  # Estimate
  estimate <- L %*% coef(gls.mod)

  # Covariance of contrast
  var_con <- qform(c(L),vcov(gls.mod)) # variance of contrast

  # Standard errors
  se <- sqrt(var_con)

  # t statistic
  tstat <- estimate/se

  # Calc for DF
  grad_var_con <- vapply(G, function(x) qform(c(L), x), numeric(1L))

  # Denominator of Satt. Calc
  satt_denom <- qform(grad_var_con, gls.mod$apVar)

  # DDF
  ddf <- drop(2 * var_con^2 / satt_denom) # denominator DF


  if(joint==T){
    return(mk_Ftable_gls(Fvalue=tstat^2, ndf=1L, ddf=ddf))
  }else{return(mk_ttable_gls(estimate = estimate, se = se, ddf = ddf))}

}

# Combined function for sinlge or multiple DF tests
satt_test_gls <- function(L, gls.mod, data, joint=FALSE, eps = sqrt(.Machine$double.eps)){
  # Calculate basic quantities
  # Covariance of fixed effects
  V <- vcov(gls.mod)

  # Gradient
  G <- (gls_grad(gls.mod, gls.mod$call, data, V))

  # If L is a single row / single contrast
  if (nrow(L)==1){
    if(joint==F){satt_gls_1df(L, gls.mod, G, joint = FALSE)
    }else{satt_gls_1df(L, gls.mod, G, joint = TRUE)}
  }else{# Continue with multiple row contrasts
    if(joint==FALSE){
      cont_names <- rownames(L)
      if(is.null(cont_names)){cont_names <- 1:nrow(L)}
      temp <- lapply(cont_names, function(nm){satt_gls_1df(L = L[nm,], gls.mod, G, joint = FALSE)})
      res <- do.call(rbind, temp)
      rownames(res) <- cont_names
      return(res)
    }else{# Perform multiple DF F test
      # Compute Var(L beta) and eigen-decompose:
      VLbeta <- L %*% V %*% t(L) # Var(contrast) = Var(Lbeta)
      eig_VLbeta <- eigen(VLbeta)
      P <- eig_VLbeta$vectors
      d <- eig_VLbeta$values
      tol <- max(eps * d[1], 0)
      pos <- d > tol
      q <- sum(pos) # rank(VLbeta)
      if(q < nrow(L))
        warning("Contrast is rank deficient and test may be affected")
      if(q <= 0) { # shouldn't happen if L is a proper contrast
        return(mk_Ftable_gls(NA, NA, NA))
      }else{
        PtL <- crossprod(P, L)[1:q, ]
        if(q == 1) { # 1D case:
          return(satt_gls_1df(PtL, gls.mod, G, joint = FALSE)) # Need to check on this...
        }else{# multi-D case proceeds:
          # Compute t-squared values and F-value:
          t2 <- drop(PtL %*% coef(gls.mod))^2 / d[1:q]
          Fvalue <- sum(t2) / q
          # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
          grad_PLcov <- lapply(1:q, function(m) {
            vapply(G, function(J) qform(PtL[m, ], J), numeric(1L))
          })
          # Compute degrees of freedom for the q t-statistics:
          nu_m <- vapply(1:q, function(m) {
            2*(d[m])^2 / qform(grad_PLcov[[m]], gls.mod$apVar) }, numeric(1L)) # 2D_m^2 / g'Ag
          # Compute ddf for the F-value:
          ddf <- get_Fstat_ddf_gls(nu_m, tol=1e-8)
          mk_Ftable_gls(Fvalue, ndf=q, ddf=ddf)
        }
      }
    }
  }
}

satt_test_gls_lmerSeq <- function(L, gls.mod, data, joint=FALSE, eps = sqrt(.Machine$double.eps), form_fit){
  # Calculate basic quantities
  # Covariance of fixed effects
  V <- vcov(gls.mod)

  # Gradient
  # G <- (gls_grad(gls.mod, gls.mod$call, data, V))
  G <- (gls_grad_lmerSeq(object = gls.mod,
                         call = gls.mod$call,
                         data = data,
                         V = V,
                         form_fit = form_fit))

  # If L is a single row / single contrast
  if (nrow(L)==1){
    if(joint==F){satt_gls_1df(L, gls.mod, G, joint = FALSE)
    }else{satt_gls_1df(L, gls.mod, G, joint = TRUE)}
  }else{# Continue with multiple row contrasts
    if(joint==FALSE){
      cont_names <- rownames(L)
      if(is.null(cont_names)){cont_names <- 1:nrow(L)}
      temp <- lapply(cont_names, function(nm){satt_gls_1df(L = L[nm,], gls.mod, G, joint = FALSE)})
      res <- do.call(rbind, temp)
      rownames(res) <- cont_names
      return(res)
    }else{# Perform multiple DF F test
      # Compute Var(L beta) and eigen-decompose:
      VLbeta <- L %*% V %*% t(L) # Var(contrast) = Var(Lbeta)
      eig_VLbeta <- eigen(VLbeta)
      P <- eig_VLbeta$vectors
      d <- eig_VLbeta$values
      tol <- max(eps * d[1], 0)
      pos <- d > tol
      q <- sum(pos) # rank(VLbeta)
      if(q < nrow(L))
        warning("Contrast is rank deficient and test may be affected")
      if(q <= 0) { # shouldn't happen if L is a proper contrast
        return(mk_Ftable_gls(NA, NA, NA))
      }else{
        PtL <- crossprod(P, L)[1:q, ]
        if(q == 1) { # 1D case:
          return(satt_gls_1df(PtL, gls.mod, G, joint = FALSE)) # Need to check on this...
        }else{# multi-D case proceeds:
          # Compute t-squared values and F-value:
          t2 <- drop(PtL %*% coef(gls.mod))^2 / d[1:q]
          Fvalue <- sum(t2) / q
          # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
          grad_PLcov <- lapply(1:q, function(m) {
            vapply(G, function(J) qform(PtL[m, ], J), numeric(1L))
          })
          # Compute degrees of freedom for the q t-statistics:
          nu_m <- vapply(1:q, function(m) {
            2*(d[m])^2 / qform(grad_PLcov[[m]], gls.mod$apVar) }, numeric(1L)) # 2D_m^2 / g'Ag
          # Compute ddf for the F-value:
          ddf <- get_Fstat_ddf_gls(nu_m, tol=1e-8)
          mk_Ftable_gls(Fvalue, ndf=q, ddf=ddf)
        }
      }
    }
  }
}

# Combined function for single or multiple DF tests
satt_test_gls_lava <- function(L, gls.mod, data, joint=FALSE, eps = sqrt(.Machine$double.eps)){


  object = gls.mod
  # Do not use the small sample bias correction from lavaSearch
  # suppressWarnings(sCorrect(object, df = TRUE, numeric.derivative=TRUE) <- FALSE)

  ######################################
  # Get Coefficient Names
  ######################################
  ## *** mean coefficients
  mean.coef <- stats::coef(object)

  var.coef <- c(sigma2 = stats::sigma(object)^2)
  if(!is.null(object$modelStruct$varStruct)){
    var.coef <- c(var.coef,
                  stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE)^2)
  }

  if(!is.null(object$modelStruct$corStruct)){
    cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)

    n.var <- length(var.coef)
    n.cor <- length(cor.coef)

    ## check unstructured covariance matrix
    if(!is.null(object$modelStruct$varStruct) && ((n.var*(n.var-1))/2 == n.cor)){

      vecgroup <- attr(unclass(object$modelStruct$corStruct), "group")
      veccov.cor <- unname(unlist(attr(object$modelStruct$corStruct, "covariate")))
      veccov.var <- attr(object$modelStruct$varStruct, "groups")

      table.covvar <- table(veccov.cor,veccov.var)
      newlevels.cor <- colnames(table.covvar)[apply(table.covvar, 1, which.max)]
      veccov.cor2 <- factor(veccov.cor, levels = 0:max(veccov.cor), labels = newlevels.cor)

      if(identical(as.character(veccov.cor2),as.character(veccov.var))){

        cor.coefName <- apply(utils::combn(newlevels.cor, m = 2), MARGIN = 2, FUN = paste, collapse = "")
        names(cor.coef) <- paste0("corCoef",cor.coefName)

      }else{
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
      }


    }else{
      names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }


  }else{
    cor.coef <- NULL
  }

  p <- c(mean.coef, cor.coef, var.coef)
  attr(p, "mean.coef") <- names(mean.coef)
  attr(p, "var.coef") <- names(var.coef)
  attr(p, "cor.coef") <- names(cor.coef)

  ###############################################
  # Update and create names for contrast mat
  # Create the null hypothesis vector
  ###############################################
  # lavaSearch2 requires the var-cov params to be included in the contrast mat
  L_new <- L

  for(i in 1:(length(var.coef)+length(cor.coef))){L_new <- cbind(L_new, 0)}

  colnames(L_new) <- names(p)

  rhs <- rep(0, nrow(L_new))

  ###############################################
  # Perform the test
  ##3############################################
  resTest <- lavaSearch2::compare2(object=object,
                      contrast = L_new,
                      null=rhs, df=T,
                      robust=F)


  # If L is a single row / single contrast
  if (nrow(L_new)==1){

    mk_ttable_gls(estimate = resTest$estimate[,1], se = resTest$estimate[,2], ddf = resTest$parameter)

  }else{
    # Perform multiple DF F test
    temp <- gsub('df1 = ', '', names(resTest$parameter))
    temp <- as.numeric(gsub(', df2', '', temp))

    mk_Ftable_gls(Fvalue=resTest$statistic, ndf=temp, ddf=resTest$parameter)
  }
}





