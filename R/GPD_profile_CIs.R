


#' Compute extreme quantile from GPD parameters
#'
#' @param quantile_lvl Probability level of the desired extreme quantile.
#' @param threshold_lvl Probability level of the GPD threshold.
#' @param threshold GPD threshold value.
#' @param scale Value(s) for the GPD scale parameter.
#' @param shape Value(s) for the GPD shape parameter.
#'
#' @return The quantile value at probability level `quantile_lvl`.
#' @export
#'
#' @examples GPD_quantiles(quantile_lvl=0.999, threshold_lvl=0.95, threshold=0, scale=1, shape=0.1)
GPD_quantiles <- function(quantile_lvl, threshold_lvl, threshold, scale, shape){
  if(any(shape==0)){
    gpd_qs <- rep(as.double(NA), max(length(threshold),length(scale),length(shape)))
    if(length(shape)==1){shape <- rep(shape,length(gpd_qs))}
    gpd_qs[shape!=0] <- ((((1-quantile_lvl)/(1-threshold_lvl))^{-shape} - 1) * (scale / shape) + threshold)[shape!=0]
    gpd_qs[shape==0] <- (log((1-threshold_lvl)/(1-quantile_lvl)) * scale + threshold)[shape==0]
  }else{
    gpd_qs <- (((1-quantile_lvl)/(1-threshold_lvl))^{-shape} - 1) * (scale / shape) + threshold
  }
  return(gpd_qs)
}


# #' GPD quantile
# #'
# #' @param threshold GPD threshold value.
# #' @param scale Scale parameter.
# #' @param shape Shape parameter.
# #' @param quantile_lvl Probability level for the quantile (1 - 1/obs_return_period).
# #' @param threshold_lvl Probability level of the threshold `threshold`.
# #'
# #' @return The quantile at probability level `quantile_lvl`, for the specified GPD distribution.
# #' @export
# #'
# #' @examples GPD_quantile(threshold=0, scale=1, shape=0.1, quantile_lvl=0.999, threshold_lvl=0.95)
# GPD_quantile <- function(threshold=0, scale=1, shape, quantile_lvl, threshold_lvl=0){
#   # evd::qgev(1-p, loc=loc, scale=scale, shape=shape)
#   # GPD: return(threshold - scale/shape * (1-(-log(1-quantile_lvl))^(-shape)))
#   return(GPD_quantiles(p=quantile_lvl, p0=threshold_lvl, threshold=threshold, scale=scale, shape=shape))
# }

#' GPD endpoint
#'
#' @param threshold GPD threshold value.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Other optional parameters for the internal endpoint function (Unused).
#' @param endpoint_type Whether to compute the upper endpoint (returns `Inf` if `shape>=0`),
#' or the unrestricted endpoint (returns the lower GPD endpoint if `shape>=0`).
#'
#' @return The endpoint of the specified GPD distribution.
#' @export
#'
#' @examples GPD_endpoint(threshold=0, scale=1, shape=-0.1, endpoint_type='upper')
GPD_endpoint <- function(threshold=0, scale=1, shape, ..., endpoint_type=c("unrestricted", "upper")){
  endpoint_type <- match.arg(endpoint_type)
  if(endpoint_type=="upper"){endpt_fct <- GPD_upper_endpoint}
  if(endpoint_type=="unrestricted"){endpt_fct <- GPD_unrestricted_endpoint}
  return(endpt_fct(threshold=threshold, scale=scale, shape=shape, ...))
}

#' Unrestricted GPD endpoint
#'
#' @param threshold GPD threshold value.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Unused.
#'
#' @return The unrestricted endpoint of the specified GPD distribution
#' (returns the upper GPD endpoint if `shape<0`, and the lower endpoint otherwise).
#'
#' @keywords internal
GPD_unrestricted_endpoint <- function(threshold=0, scale=1, shape, ...){
  # @examples GPD_unrestricted_endpoint(threshold=0, scale=1, shape=-0.1)
  return(threshold - scale/shape)
}

#' GPD upper endpoint
#'
#' @param threshold GPD threshold value.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Unused.
#'
#' @return The upper endpoint of the specified GPD distribution (returns `Inf` if `shape>=0`).
#'
#' @keywords internal
GPD_upper_endpoint <- function(threshold=0, scale=1, shape, ...){
  # @examples GPD_upper_endpoint(threshold=0, scale=1, shape=-0.1)
  if(length(shape==1) & (length(threshold)>1 | length(scale)>1)){
    shape <- rep(shape, max(length(threshold),length(scale)))
  }
  endpts <- rep(Inf, length(shape))
  endpts[shape<0] <- (threshold - scale/shape)[shape<0]
  return(endpts)
}

#' GPD log likelihood
#'
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#'
#' @return The GPD log-likelihood evaluated at the given parameters, given the data observations `Y`.
#' @export
GPD_log_likelihood <- function(Y, threshold, scale, shape, obs_weights=NULL, ill_defined_value=-10^6){
  
  #TODO: ill_defined_value -Inf or -1e16 -1e15 ?
  
  exceeds <- (Y > threshold)
  Z <- (Y[exceeds]-threshold)
  
  resc <- 1+shape*Z/scale
  # cat('DEBUG: min resc=', min(resc), '\n')
  
  # w <- if(is.null(obs_weights)){rep(1, length(Z))}else{obs_weights[exceeds]}
  
  if(is.null(obs_weights)){
    # Unweighted log-likelihood
    
    # if(any(resc <= 0) | any(scale <= 0)){return(ill_defined_value)}
    if(any(resc <= 1e-15) | any(scale <= 1e-15)){return(ill_defined_value)}
    
    # if(length(scale)==1 & length(shape)==1){
    #   return(-length(Z)*log(scale) - (1+1/shape)*sum(log(resc)))
    # }else if(length(scale)>1){
    #   return(-sum(log(scale)) - sum(log(resc)*(1/shape + 1)))
    # }else{
    #   return(-length(Z)*log(scale) - sum(log(resc)*(1/shape + 1)))
    # }
    
    if(length(scale)==1){
      term1 <- (-length(Z)*log(scale))
    }else if(length(scale)>1){
      term1 <- (-sum(log(scale)))
    }else{stop('Error with length(scale) in GPD_log_likelihood.')}
    
    if(length(shape)==1){
      term2 <- (-(1+1/shape)*sum(log(resc)))
    }else if(length(shape)>1){
      term2 <- (-sum(log(resc)*(1/shape + 1)))
    }else{stop('Error with length(shape) in GPD_log_likelihood.')}
    
  }else{
    # Weighted log-likelihood
    if(length(obs_weights)!=length(Y)){stop('obs_weights should be of the same length as Y.')}
    
    w <- obs_weights[exceeds]
    
    if(any(w<0)){stop('obs_weights should be non-negative.')}
    if(sum(w)==0){stop('exceedence obs_weights cannot be all zero.')}
    
    # cat('DEBUG: sum(obs_weights)=', sum(obs_weights), '\n')
    # cat('DEBUG: sum(w)=', sum(w), '\n')
    # cat('DEBUG: scale=', scale, '\n')
    # cat('DEBUG: shape=', shape, '\n')
    
    #TODO: normalize weights?
    # w <- w / sum(w) * length(w) # sum(w) = length(w)
    
    nonzero_w_ids <- (w!=0)
    w <- w[nonzero_w_ids]
    resc <- resc[nonzero_w_ids]
    
    # if(any(resc[w!=0] <= 0)){return(ill_defined_value)}
    if(any(resc <= 1e-15)){return(ill_defined_value)}
    
    
    if(length(scale)==1){
      # if(scale <= 0){return(ill_defined_value)}
      if(scale <= 1e-15){return(ill_defined_value)}
      # cat('DEBUG: => -sum(w)*log(scale)=', -sum(w)*log(scale), '\n')
      term1 <- (-sum(w)*log(scale))
    }else if(length(scale)>1){
      scale <- scale[nonzero_w_ids]
      # if(any(scale <= 0)){return(ill_defined_value)}
      if(any(scale <= 1e-15)){return(ill_defined_value)}
      term1 <- (-sum(w*log(scale)))
    }else{stop('Error with length(scale) in GPD_log_likelihood.')}
    
    if(length(shape)==1){
      # cat('DEBUG: => (1+1/shape)=', (1+1/shape), '\n')
      # cat('DEBUG: => sum(w*log(resc))=', -sum(w)*log(scale) - (1+1/shape)*sum(w*log(resc)), '\n')
      term2 <- (-(1+1/shape)*sum(w*log(resc)))
    }else if(length(shape)>1){
      shape <- shape[nonzero_w_ids]
      term2 <- (-sum(w*log(resc)*(1/shape + 1)))
    }else{stop('Error with length(shape) in GPD_log_likelihood.')}
    
  }
  # cat('DEBUG: ==> ll=', term1 + term2, '\n')
  return(term1 + term2)
}

# GEV_log_likelihood_rlvl <- function(Y, return_lvl, scale, shape, p){
#   loc <- return_lvl + scale/shape * (1-(-log(1-p))^(-shape))
#   return(GEV_log_likelihood(Y, loc, scale, shape)) #TODO: simplify?
# }
#
# GEV_log_likelihood_endpt <- function(Y, endpoint, scale, shape){
#   loc <- endpoint + scale/shape
#   return(GEV_log_likelihood(Y, loc, scale, shape))
# }

#' Internal GPD log-likelihood function format for optimisation
#'
#' @param a .
#' @param Y .
#' @param threshold .
#' @param threshold_lvl .
#' @param parametrization .
#' @param quantile_lvl .
#' @param scamat .
#' @param shamat .
#' @param negative .
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#'
#' @keywords internal
GPD_log_likelihood_optim <- function(a, Y, threshold=0, threshold_lvl=0, parametrization=c("classical", "orthogonal", "quantile", "endpoint"), quantile_lvl=1.-(1./100.), # , "othogonal"
                                     scamat=as.matrix(1), shamat=as.matrix(1), negative=FALSE, obs_weights=NULL, ill_defined_value=-10^6){
  parametrization <- match.arg(parametrization)
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  sc_pars <- a[1:nbsca]
  sh_pars <- a[nbsca+(1:nbsha)]
  
  if(parametrization=="quantile"){
    # GEV: lo_pars[1] <- a[1] + sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
    sc_pars[1] <- ((a[1] - threshold)*sh_pars[1])/(((1-threshold_lvl)/(1-quantile_lvl))^{sh_pars[1]} - 1) # r = u + s/x * (((1-p0)/(1-p))^x - 1) -> s = ((r-u)*x)/(((1-p0)/(1-p))^x - 1)
  }
  if(parametrization=="endpoint"){
    # GEV: lo_pars[1] <- a[1] + sc_pars[1]/sh_pars[1]
    sc_pars[1] = (threshold - a[1])*sh_pars[1] # e = u - s/x -> s = (u - e)*x
  }
  
  scale <- c(scamat %*% sc_pars)
  shape <- c(shamat %*% sh_pars)
  
  ll <- GPD_log_likelihood(Y=Y, threshold=threshold, scale=scale, shape=shape, 
                           obs_weights=obs_weights, ill_defined_value=ill_defined_value)
  
  if(negative){
    return(-ll)
  }else{
    return(ll)
  }
}


#' Maximum-likelihood GPD estimate
#'
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parametrization Likelihood parametrization. Alternatives to `classical` substitute for the scale parameter.
#' @param quantile_lvl Quantile probability level for the `'quantile'` parametrization.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits).
#' @param x_rlvl Covariate vector at which to reparametrize for the `'quantile'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param out_param Additional output parametrization (same as `parametrization`, by default).
#' If `out_param != parametrization`, the parameters are reparametrized from `parametrization` to `out_param` after estimation, in a separate output.
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param verbose Verbose level, as integer.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @return The fitted maximum-likelihood GPD as a `GPD_ML` object, containing:
#' \item{mle}{The estimated maximum likelihood GPD parameters, as a named vector (expressed in `parametrization`).}
#' \item{loglik}{The log-likelihood of the estimated parameters, given the data}
#' \item{conv}{Whether the optimisation procedure converged.
#' See the `convergence` output of [stats::optim()] for more details.}
#' \item{hessian}{The loglikelihood hessian evaluated at the estimated parameters, given the data.}
#' \item{parametrization}{Likelihood parametrization.}
#' \item{out_mle}{The estimated maximum likelihood GPD parameters, reparametrized in `out_param`.}
#' \item{out_parametrization}{Additional output parametrization.}
#' @export
GPD_maxlik <- function(Y, threshold=0, threshold_lvl=0, parametrization=c("classical", "orthogonal", "quantile", "endpoint"), 
                       quantile_lvl=1.-(1./100.), orthogonal=FALSE,
                       X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL, out_param=parametrization, 
                       obs_weights=NULL, ill_defined_value=-10^6,
                       hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"),
                       verbose=1, ...){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "orthogonal", "quantile", "endpoint"))
  method <- match.arg(method)
  
  if(parametrization=='orthogonal'){stop('orthogonal argument deprecated.')}
  if(out_param=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  #Dealing with covariates
  feat_mats <- GPD_covariate_matrices(X=X, x_rlvl=x_rlvl, scale_cols=scale_cols, shape_cols=shape_cols,
                                      parametrization=parametrization, out_param=out_param)
  # locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  # nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  # TODO : inner parametrization instead of out_param? Here best?
  
  # Initialization
  init <- GPD_param_init(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, nbsca=nbsca, nbsha=nbsha, 
                         parametrization=parametrization, quantile_lvl=quantile_lvl, obs_weights=obs_weights)
  
  if(verbose>=10){{cat("\nDEBUG: GPD_maxlik: init value optim =", init, "\n")}}
  
  sol <- stats::optim(par=init, fn=GPD_log_likelihood_optim,
                      Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=parametrization, quantile_lvl=quantile_lvl,
                      scamat=scamat, shamat=shamat, negative=TRUE,
                      obs_weights=obs_weights, ill_defined_value=ill_defined_value,
                      method = method, control = list(maxit = maxit, ...), hessian = hessian)
  mle <- sol$par
  out_mle <- GPD_change_parametrization(mle, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=parametrization, new_parametrization=out_param,
                                        quantile_lvl=quantile_lvl, nbsca=nbsca, nbsha=nbsha)
  
  # mle names
  names(mle) <- GPD_param_names(nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  
  out <- list(mle=mle, loglik=-sol$value, conv=sol$convergence, hessian=-sol$hessian,
              parametrization=parametrization, out_mle=out_mle, out_parametrization=out_param)
  class(out) <- c("GPD_ML")
  
  return(out)
}


#' Optimisation evaluation step for the GPD profile likelihood
#'
#' @param a .
#' @param val .
#' @param Y .
#' @param threshold .
#' @param threshold_lvl .
#' @param parametrization .
#' @param id_param .
#' @param quantile_lvl .
#' @param scamat .
#' @param shamat .
#' @param negative .
#' @param orthogonal .
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#'
#' @keywords internal
optim_step_GPD_profile <- function(a, val, Y, threshold=0, threshold_lvl=0, parametrization=c("classical", "orthogonal", "quantile", "endpoint"),
                                   id_param, quantile_lvl=1.-(1./100.),
                                   scamat=as.matrix(1), shamat=as.matrix(1), negative=TRUE, orthogonal=FALSE, 
                                   obs_weights=NULL, ill_defined_value=-10^6){
  parametrization <- match.arg(parametrization)
  
  pars_all <- vector_insert(a, val, id_param)
  
  pll <- GPD_log_likelihood_optim(a=pars_all, Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=parametrization, quantile_lvl=quantile_lvl,
                                  scamat=scamat, shamat=shamat, negative=negative, obs_weights=obs_weights, ill_defined_value=ill_defined_value)
  return(pll)
}


#' Internal function for the GPD profile likelihood
#'
#' @param val .
#' @param Y .
#' @param threshold .
#' @param threshold_lvl .
#' @param parameter .
#' @param quantile_lvl .
#' @param scamat .
#' @param shamat .
#' @param subparam_id .
#' @param orthogonal .
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param init .
#' @param hessian .
#' @param maxit .
#' @param method .
#' @param verbose Verbose level, as integer.
#' @param ... .
#'
#' @keywords internal
GPD_profile_loglik_internal <- function(val, Y, threshold=0, threshold_lvl=0, parameter=c("shape", "scale", "quantile", "endpoint"), quantile_lvl=1.-(1./100.),
                                        scamat=as.matrix(1), shamat=as.matrix(1), subparam_id=0,
                                        orthogonal=FALSE, obs_weights=NULL, ill_defined_value=-10^6, init=NULL,
                                        hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"),
                                        verbose=1, ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  parametrization <- switch(parameter, "shape"="classical", "scale"="classical",
                            "quantile"="quantile", "endpoint"="endpoint")
  id_param <- GPD_profpar_id(parameter=parameter, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  
  # Initialization
  if(is.null(init)){
    init <- GPD_param_init(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, nbsca=nbsca, nbsha=nbsha, 
                           parametrization=parametrization, quantile_lvl=quantile_lvl, obs_weights=obs_weights)
    init <- init[-id_param]
  }
  #TODO: else{ check init format/dimensions. }
  
  if(verbose>=10){{cat("\nDEBUG: GPD_profile_loglik_internal: init value optim =", init, "\n")}}
  
  sol <- stats::optim(par=init, fn=optim_step_GPD_profile,
                      val=val, Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=parametrization, 
                      id_param=id_param, quantile_lvl=quantile_lvl,
                      scamat=scamat, shamat=shamat, negative=TRUE, orthogonal=orthogonal, 
                      obs_weights=obs_weights, ill_defined_value=ill_defined_value, 
                      method = method, control = list(maxit = maxit, ...), hessian = hessian)
  mle_other <- sol$par
  
  # mle names
  params_names <- GPD_param_names(nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  names(mle_other) <- params_names[-id_param]
  
  param_val <- val
  names(param_val) <- params_names[id_param]
  
  out <- list(param_val=param_val, param_name=params_names[id_param], mle_other=mle_other,
              loglik=-sol$value, conv=sol$convergence, hessian=-sol$hessian,
              parameter=parameter, parametrization=parametrization, subparam_id=subparam_id, id_param=id_param)
  class(out) <- c("GPD_profML")
  return(out)
}


#' GPD profile log-likelihood
#'
#' @param val Parameter value at which to evaluate the GPD profile log-likelihood.
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param quantile_lvl Quantile probability level for the `'quantile'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'quantile'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param init Optional initial values for the remaining parameter's optimisation process, in the correct internal format.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @return The GPD profile log-likelihood of `parameter` evaluated at `val`, given the data,
#' as a `GPD_profML` object containing:
#' \item{param_val}{(Named) parameter value at which the GPD profile log-likelihood was evaluated.}
#' \item{param_name}{Name of the evaluated profile likelihood parameter.}
#' \item{mle_other}{Maximum-likelihood estimate of the other GPD parameters.}
#' \item{loglik}{Profile GPD log-likelihood value of `parameter` evaluated at `val`, given the data.}
#' \item{conv}{Whether the optimisation procedure converged.
#' See the `convergence` output of [stats::optim()] for more details.}
#' \item{hessian}{The loglikelihood hessian evaluated at the estimated parameters, given the data.}
#' \item{parameter}{Name of the evaluated profile likelihood parameter given as argument (redundent).}
#' \item{parametrization}{Likelihood parametrization.}
#' \item{subparam_id}{Index of the parameter coefficient for which the profile log-likelihood was computed.}
#' \item{id_param}{Index of the likelihood profile parameter, in the internal parameter vector format.}
#' @export
GPD_profile_loglik <- function(val, Y, threshold=0, threshold_lvl=0, parameter=c("shape", "scale", "quantile", "endpoint"),
                               subparam_id=0, quantile_lvl=1.-(1./100.), orthogonal=FALSE,
                               X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL, obs_weights=NULL, ill_defined_value=-10^6,
                               init=NULL, hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  
  #Dealing with covariates
  parametrization <- switch(parameter, "shape"="classical", "scale"="classical",
                            "quantile"="quantile", "endpoint"="endpoint")
  feat_mats <- GPD_covariate_matrices(X=X, x_rlvl=x_rlvl, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  # locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  # nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  out <- GPD_profile_loglik_internal(val=val, Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parameter=parameter, quantile_lvl=quantile_lvl,
                                     scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                     orthogonal=orthogonal, obs_weights=obs_weights, ill_defined_value=ill_defined_value, init=init,
                                     hessian=hessian, maxit=maxit, method=method, ...)
  return(out)
}


#' GPD profile log-likelihood curve
#'
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param alpha Confidence alpha for the profile log-likelihood confidence interval (i.e. for the confidence line on the profile plot).
#' @param quantile_lvl Quantile probability level for the `'quantile'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'quantile'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param warmstart_table Evaluation table from a previous run.
#' @param stepsize Numerical size of each evaluation step, in the profile parameter's scale.
#' @param steps_beyond_conf Number of additional steps to take (in each direction)
#' after the profile log-likelihood values reach below the confidence line.
#' @param initial_MLE_para Parametrization used for the initial maximum likelihood estimate (defaults to classical, for better stability).
#' @param max_steps Maximum number of steps taken (in each direction).
#' If the confidence line was not reached, the corresponding confidence interval endpoint will be infinite.
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used for the initial maximum likelihood optimisation. See [stats::optim()] for more details.
#' @param method_prof The optimisation method to be used for the profile likelihood optimisation. See [stats::optim()] for more details.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @return The GPD profile log-likelihood curve for the desired `parameter`,
#' with confidence line and resulting `(1-alpha)` confidence interval, as a `GPD_profileLogLik` object containing:
#' \item{mle}{The estimated maximum likelihood GPD parameters, as a named vector
#' (expressed in the profile parametrization).}
#' \item{ci}{Length-two vector containing the lower and upper endpoints of the
#' desired profile likelihood confidence interval.}
#' \item{profile_loglik}{Named matrix containing the profile loglikelihood value (Column 2)
#' for each considered profile parameter value (Column 1).}
#' \item{conf_line}{Confidence line for the desired profile likelihood confidence interval.
#' See e.g. Coles (2001) for more details.}
#' \item{eval_table}{Tibble ([tibble::tibble()]) containing the history of
#' profile log-likelihood evaluation values, and related metadata.}
#' \item{param_name}{Name of the profiled parameter (infered, for debugging purposes).}
#' \item{parameter}{Name of the profiled parameter (given).}
#' \item{parametrization}{Parametrization used for the profile likelihood.}
#' \item{subparam_id}{Index of the parameter coefficient for which the profile log-likelihood was computed.}
#' \item{id_param}{Index of the profile parameter in the GPD parameter vector.}
#' @export
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer. \doi{doi:10.1007/978-1-4471-3675-0}.
GPD_profile_loglik_curve <- function(Y, threshold=0, threshold_lvl=0, parameter=c("shape", "scale", "quantile", "endpoint"),
                                     subparam_id=0, alpha=0.05, quantile_lvl=1.-(1./100.), orthogonal=FALSE,
                                     X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL, warmstart_table=NULL,
                                     stepsize=0.1, steps_beyond_conf=5, initial_MLE_para=c("classical", "same"), max_steps=1e4, 
                                     obs_weights=NULL, ill_defined_value=-10^6,
                                     hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                                     method_prof=c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  method_prof <- match.arg(method_prof)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parametrization <- switch(parameter, "shape"="classical", "scale"="classical",
                            "quantile"="quantile", "endpoint"="endpoint")
  inner_pram <- if(initial_MLE_para=="same"){parametrization}else{initial_MLE_para}
  
  if(method_prof=="default"){
    if(method=="Nelder-Mead" & (is.null(scale_cols) & is.null(shape_cols))){
      # Nelder-Mead does not work well with one-dimensional pars
      method_prof <- "BFGS"
    } else {
      method_prof <- method
    }
  }
  
  sol_mle <- GPD_maxlik(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=inner_pram, quantile_lvl=quantile_lvl, orthogonal=orthogonal,
                        X=X, x_rlvl=x_rlvl, scale_cols=scale_cols, shape_cols=shape_cols,
                        out_param=parametrization, obs_weights=obs_weights, ill_defined_value=ill_defined_value, hessian=hessian, maxit=maxit, method=method, ...)
  # mle and its log-likelihood
  mle <- sol_mle$out_mle
  ll_mle <- sol_mle$loglik
  
  # covariate matrices
  feat_mats <- GPD_covariate_matrices(X=X, x_rlvl=x_rlvl, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  id_param <- GPD_profpar_id(parameter=parameter, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  params_names <- GPD_param_names(nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  param_name <- params_names[id_param]
  
  # set up eval table
  eval_table <- tibble::as_tibble(as.list(c(mle, log_lik=ll_mle)))
  if(!is.null(warmstart_table)){
    eval_table <- dplyr::full_join(eval_table, warmstart_table, by=colnames(eval_table)) %>%
      dplyr::arrange(tidyselect::all_of(param_name))
  }
  
  # chi-squared 1-alpha quantile (df 1)
  chi2_q <- stats::qchisq(p=(1-alpha), df=1)
  
  conf_line <- ll_mle - (chi2_q/2.0) + 1e-15 # add a small epsilon to avoid numerical issues
  # 2*(ll_mle-plli) <= chi2_q
  # plli >= conf_line
  
  vali <- mle[[id_param]]
  plli <- ll_mle
  nb_below <- 0
  going_down <- TRUE
  
  stop_thresh_down <- vali - max_steps*stepsize
  stop_thresh_up <- vali + max_steps*stepsize
  
  while(going_down | nb_below<=steps_beyond_conf){
    valim1 <- vali
    pllim1 <- plli
    if(going_down){
      vali <- vali - stepsize
    }else{
      vali <- vali + stepsize
    }
    
    queried <- eval_table[eval_table[[param_name]]==vali,]
    if(nrow(queried)==1){
      plli <- queried$log_lik
    }else{
      out <- GPD_profile_loglik_internal(val=vali, Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parameter=parameter, quantile_lvl=quantile_lvl,
                                         scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                         orthogonal=orthogonal, obs_weights=obs_weights, ill_defined_value=ill_defined_value, init=mle[-id_param],
                                         hessian=hessian, maxit=maxit, method=method_prof, ...)
      plli <- out$loglik
      
      if(going_down){
        eval_table <- dplyr::bind_rows(c(out$param_val, out$mle_other, log_lik=out$loglik), eval_table)
      }else{
        eval_table <- dplyr::bind_rows(eval_table, c(out$param_val, out$mle_other, log_lik=out$loglik))
      }
    }
    
    if(plli < conf_line){
      nb_below <- nb_below + 1
      if(nb_below==1){
        if(going_down){
          ci_down <- vali #vali + (conf_line-plli)/(pllim1-plli) * (valim1-vali)/2 #/2 for approximate concavity correction
        }else{
          ci_up <- vali #vali - (conf_line-plli)/(pllim1-plli) * (vali-valim1)/2 #/2 for approximate concavity correction
        }
      }
      if(nb_below>steps_beyond_conf){
        if(going_down){
          going_down <- FALSE
          vali <- mle[[id_param]]
          plli <- ll_mle
          nb_below <- 0
        }
      }
    }else{
      nb_below <- 0
      if(going_down & (vali < stop_thresh_down)){
        ci_down <- -Inf
        going_down <- FALSE
        vali <- mle[[id_param]]
        plli <- ll_mle
      }
      if(!going_down & (vali > stop_thresh_up)){
        ci_up <- Inf
        nb_below <- steps_beyond_conf+1 # TODO: cleaner add another bool variable?
      }
    }
  }
  eval_table %>% dplyr::arrange(tidyselect::all_of(param_name))
  profile_loglik <- as.matrix(eval_table[,c(params_names[id_param], "log_lik")])
  out <- list(mle=mle, ci=c(ci_down, ci_up), profile_loglik=profile_loglik, conf_line=conf_line,
              eval_table=eval_table, param_name=params_names[id_param],
              parameter=parameter, parametrization=parametrization, subparam_id=subparam_id, id_param=id_param)
  class(out) <- c("GPD_profileLogLik")
  return(out)
}

#' GPD profile CI using binary search
#'
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param alpha Confidence alpha for the profile log-likelihood confidence interval (i.e. for the confidence line on the profile plot).
#' @param quantile_lvl Quantile probability level for the `'quantile'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'quantile'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param warmstart_table Evaluation table from a previous run.
#' @param init_step_pos Initial numerical size of each evaluation step to the right, in the profile parameter's scale.
#' @param init_step_neg Initial numerical size of each evaluation step to the left, in the profile parameter's scale.
#' @param tol Numerical tolerance for convergence, in the profile parameter's scale.
#' @param steps_beyond_conf Number of additional steps to take (in each direction)
#' after the profile log-likelihood values reach below the confidence line.
#' @param initial_MLE_para Parametrization used for the initial maximum likelihood estimate (defaults to classical, for better stability).
#' @param max_steps Maximum number of steps taken (in each direction).
#' If the confidence line was not reached, the corresponding confidence interval endpoint will be infinite.
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used for the initial maximum likelihood optimisation. See [stats::optim()] for more details.
#' @param method_prof The optimisation method to be used for the profile likelihood optimisation. See [stats::optim()] for more details.
#' @param verbose Verbose level, as integer.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @return The GPD profile log-likelihood confidence interval for the desired `parameter`,
#' with confidence line and resulting `(1-alpha)` confidence interval, as a `GPD_profileLogLik` object containing:
#' \item{mle}{The estimated maximum likelihood GPD parameters, as a named vector
#' (expressed in the profile parametrization).}
#' \item{ci}{Length-two vector containing the lower and upper endpoints of the
#' desired profile likelihood confidence interval.}
#' \item{profile_loglik}{Named matrix containing the profile loglikelihood value (Column 2)
#' for each considered profile parameter value (Column 1).}
#' \item{conf_line}{Confidence line for the desired profile likelihood confidence interval.
#' See e.g. Coles (2001) for more details.}
#' \item{eval_table}{Tibble ([tibble::tibble()]) containing the history of
#' profile log-likelihood evaluation values, and related metadata.}
#' \item{param_name}{Name of the profiled parameter (infered, for debugging purposes).}
#' \item{parameter}{Name of the profiled parameter (given).}
#' \item{parametrization}{Parametrization used for the profile likelihood.}
#' \item{subparam_id}{Index of the parameter coefficient for which the profile log-likelihood was computed.}
#' \item{id_param}{Index of the profile parameter in the GPD parameter vector.}
#' @export
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer. \doi{doi:10.1007/978-1-4471-3675-0}.
GPD_profile_CI <- function(Y, threshold=0, threshold_lvl=0, parameter=c("shape", "scale", "quantile", "endpoint"),
                           subparam_id=0, alpha=0.05, quantile_lvl=1.-(1./100.), orthogonal=FALSE,
                           X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL, warmstart_table=NULL,
                           init_step_pos=100, init_step_neg=10, tol=0.01, steps_beyond_conf=5, initial_MLE_para=c("classical", "same"), max_steps=1e3, 
                           obs_weights=NULL, ill_defined_value=-10^6,
                           hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                           method_prof=c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), verbose=1, ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  method_prof <- match.arg(method_prof)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parametrization <- switch(parameter, "shape"="classical", "scale"="classical",
                            "quantile"="quantile", "endpoint"="endpoint")
  inner_pram <- if(initial_MLE_para=="same"){parametrization}else{initial_MLE_para}
  
  if(method_prof=="default"){
    if(method=="Nelder-Mead" & (is.null(scale_cols) & is.null(shape_cols))){
      # Nelder-Mead does not work well with one-dimensional pars
      method_prof <- "BFGS"
    } else {
      method_prof <- method
    }
  }
  
  if(verbose>=5){cat("DEBUG: START MLE\n")}
  
  sol_mle <- GPD_maxlik(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parametrization=inner_pram, quantile_lvl=quantile_lvl, orthogonal=orthogonal,
                        X=X, x_rlvl=x_rlvl, scale_cols=scale_cols, shape_cols=shape_cols,
                        out_param=parametrization, obs_weights=obs_weights, ill_defined_value=ill_defined_value, hessian=hessian, maxit=maxit, method=method,
                        verbose=verbose, ...)
  # mle and its log-likelihood
  mle <- sol_mle$out_mle
  ll_mle <- sol_mle$loglik
  
  # covariate matrices
  feat_mats <- GPD_covariate_matrices(X=X, x_rlvl=x_rlvl, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  id_param <- GPD_profpar_id(parameter=parameter, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  params_names <- GPD_param_names(nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  param_name <- params_names[id_param]
  
  # set up eval table
  eval_table <- tibble::as_tibble(as.list(c(mle, log_lik=ll_mle)))
  if(!is.null(warmstart_table)){
    eval_table <- dplyr::full_join(eval_table, warmstart_table, by=colnames(eval_table)) %>%
      dplyr::arrange(tidyselect::all_of(param_name))
  }
  
  # chi-squared 1-alpha quantile (df 1)
  chi2_q <- stats::qchisq(p=(1-alpha), df=1)
  
  conf_line <- ll_mle - (chi2_q/2.0) + 1e-15 # add a small epsilon to avoid numerical issues
  # 2*(ll_mle-plli) <= chi2_q
  # plli >= conf_line
  
  vali <- mle[[id_param]]
  plli <- ll_mle
  nb_below <- 0
  going_down <- TRUE
  final_phase <- FALSE
  stepsize <- init_step_neg
  
  stop_thresh_down <- vali - max_steps*init_step_neg
  stop_thresh_up <- vali + max_steps*init_step_pos
  
  if(verbose>=5){cat("DEBUG: MLE=", vali, "; ML=", plli, "; conf_line=", conf_line, "\n ================ \n")}
  
  while(going_down | nb_below<=steps_beyond_conf){
    valim1 <- vali
    pllim1 <- plli
    final_phase <- (abs(stepsize)<tol & stepsize>0)
    
    if(going_down){
      vali <- vali - stepsize
    }else{
      vali <- vali + stepsize
    }
    
    # queried <- eval_table[eval_table[[param_name]]==vali,]
    queried <- eval_table[abs(eval_table[[param_name]]-vali)<1e-14,]
    if(nrow(queried)==1){
      plli <- queried$log_lik
    }else{
      out <- GPD_profile_loglik_internal(val=vali, Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parameter=parameter, quantile_lvl=quantile_lvl,
                                         scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                         orthogonal=orthogonal, obs_weights=obs_weights, ill_defined_value=ill_defined_value, init=mle[-id_param],
                                         hessian=hessian, maxit=maxit, method=method_prof, verbose=verbose, ...)
      plli <- out$loglik
      # TODO: catch errors and set to -Inf or check if plli is Inf/-Inf/NA/NaN/NULL
      
      if(going_down){
        eval_table <- dplyr::bind_rows(c(out$param_val, out$mle_other, log_lik=out$loglik), eval_table)
      }else{
        eval_table <- dplyr::bind_rows(eval_table, c(out$param_val, out$mle_other, log_lik=out$loglik))
      }
    }
    
    if((going_down & vali > mle[[id_param]]+1e-14) | (!going_down & vali < mle[[id_param]]-1e-14)){
      stop('In "GPD_profile_CI": CI search got back to MLE. This might be due to numerical instability, or to a large "tol" value.')
    }
    
    if(plli < conf_line){
      if(verbose>=10){cat(vali, "* (", plli, ") | ", sep = '')}
      if(final_phase){
        nb_below <- nb_below + 1
        if(nb_below==1){
          if(going_down){
            ci_down <- vali #vali + (conf_line-plli)/(pllim1-plli) * (valim1-vali)/2 #/2 for approximate concavity correction
          }else{
            ci_up <- vali #vali - (conf_line-plli)/(pllim1-plli) * (vali-valim1)/2 #/2 for approximate concavity correction
          }
        }
        if(nb_below>steps_beyond_conf){
          if(going_down){
            if(verbose>=10){cat("\n ==== DEBUG: GOING UP at vali=", vali, " ==== \n")}
            going_down <- FALSE
            vali <- mle[[id_param]]
            plli <- ll_mle
            nb_below <- 0
            stepsize <- init_step_pos
          }
        }
      }else{
        if(stepsize < 0){
          if(abs(stepsize)>=tol){
            stepsize <- stepsize/2
          }else{
            stepsize <- stepsize
          }
        }else if (stepsize > 0){
          if(abs(stepsize)>=tol){
            stepsize <- -stepsize/2
          }else{
            stop('In "GPD_profile_CI": algorithm should be in final phase but is not.')
            # final_phase <- TRUE
          }
        }else{
          stop('In "GPD_profile_CI": stepsize got to zero.')
        }
      }
    }else{
      if(verbose>=10){cat(vali, " (", plli, ") | ", sep = '')}
      nb_below <- 0
      if(stepsize < 0){
        if(abs(stepsize)>=tol){
          stepsize <- -stepsize/2
        }else{
          stepsize <- -stepsize
          final_phase <- TRUE
        }
      }
      if(going_down & (vali < stop_thresh_down)){
        if(verbose>=10){cat("\nDEBUG: MAX STEPS DOWN at vali=", vali, "\n")}
        ci_down <- -Inf
        going_down <- FALSE
        vali <- mle[[id_param]]
        plli <- ll_mle
        stepsize <- init_step_pos
      }
      if(!going_down & (vali > stop_thresh_up)){
        if(verbose>=10){cat("\nDEBUG: MAX STEPS UP at vali=", vali, "\n")}
        ci_up <- Inf
        nb_below <- steps_beyond_conf+1 # TODO: cleaner add another bool variable?
      }
    }
  }
  eval_table %>% dplyr::arrange(tidyselect::all_of(param_name))
  profile_loglik <- as.matrix(eval_table[,c(params_names[id_param], "log_lik")])
  out <- list(mle=mle, ci=c(ci_down, ci_up), profile_loglik=profile_loglik, conf_line=conf_line,
              eval_table=eval_table, param_name=params_names[id_param],
              parameter=parameter, parametrization=parametrization, subparam_id=subparam_id, id_param=id_param)
  class(out) <- c("GPD_profileLogLik")
  return(out)
}


# find_profile_ci_GEV <- function(Y, parameter=c("shape", "location", "scale", "quantile", "endpoint"),
#                                 subparam_id=0, alpha=0.05, p=1./100., orthogonal=FALSE,
#                                 X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
#                                 warmstart_table=NULL,
#                                 hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
#   parameter <- match.arg(parameter)
#   method <- match.arg(method)
#   parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
#                             "quantile"="quantile", "endpoint"="endpoint")
#
#   #TODO
#   return(list(mle=mle, ci=c(ci_down, ci_up), eval_table=eval_table))
# }


#' Multi-value conditional GPD profile likelihood confidence intervals
#'
#' @description
#' For non-stationary models, the quantile reparametrization depends on covariate values.
#' This function repeats the profile likelihood procedure for several covariate values.
#' It enables obtaining a quantile (or endpoint) curve, with profile-likelihood confidence bands,
#' as a function of the covariate values.
#'
#'
#' @param Y Data observations.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parameter Parameter for which to compute the profile likelihood confidence intervals.
#' @param alpha Confidence alpha for the profile likelihood confidence intervals.
#' @param quantile_lvl Quantile probability level for the `'quantile'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param X_rlvl Covariate matrix at which to reparametrize for the `'quantile'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' Columns should be variables, and each row should give one covariate realization at which to reparametrize and obtain a CI.
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param init_step_pos Initial numerical size of each evaluation step to the right, in the profile parameter's scale.
#' @param init_step_neg Initial numerical size of each evaluation step to the left, in the profile parameter's scale.
#' @param tol Numerical tolerance for convergence, in the profile parameter's scale.
#' @param steps_beyond_conf Number of additional steps to take (in each direction)
#' after the profile log-likelihood values reach below the confidence line.
#' @param initial_MLE_para Parametrization used for the initial maximum likelihood estimate (defaults to classical, for better stability).
#' @param max_steps Maximum number of steps taken (in each direction).
#' If the confidence line was not reached, the corresponding confidence interval endpoint will be infinite.
#' @param obs_weights Optional observation weights for weighted likelihood.
#' @param ill_defined_value Value to return if the arguments are out of support (e.g. negative scale, or non-positive arguments to logarithms).
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used for the initial maximum likelihood optimisation. See [stats::optim()] for more details.
#' @param method_prof The optimisation method to be used for the profile likelihood optimisation. See [stats::optim()] for more details.
#' @param parallel_strat Parallel strategy. One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures
#' that can be active at the same time before blocking.
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers.
#' The function should return a numeric scalar.
#' Defaults to [future::availableCores()]`-1` if `NULL` (default), with `"multicore"` constraint in the relevant case.
#' Ignored if `strategy=="sequential"`.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @return The GPD profile log-likelihood `(1-alpha)` confidence intervals for the desired `parameter`,
#' for each desired covariate values, as a [tibble::tibble()], with columns:
#' \item{obs}{Index of the observation (i.e. row) of X_rlvl for which the parameter estimate and CI was computed.}
#' \item{`<parameter name>`}{Conditional estimate of the `parameter`.}
#' \item{ci_down}{Lower endpoint of the conditional `(1-alpha)` profile-likelihood confidence interval for `parameter`.}
#' \item{ci_up}{Upper endpoint of the conditional `(1-alpha)` profile-likelihood confidence interval for `parameter`.}
#' \item{parameter}{Name of the parameter for which the estimates and CIs were computed.}
#' \item{alpha}{Confidence alpha for the profile likelihood confidence intervals.}
#' \item{quantile_lvl}{Quantile probability level for `'quantile'` parameter (only if `parameter==quantile`).}
#' @export
GPD_profile_CIs_multiple <- function(Y, threshold=0, threshold_lvl=0, parameter=c("quantile", "endpoint"),
                                     alpha=0.05, quantile_lvl=1.-(1./100.), orthogonal=FALSE,
                                     X=NULL, X_rlvl=NULL, scale_cols=NULL, shape_cols=NULL,
                                     init_step_pos=100, init_step_neg=10, tol=0.01, steps_beyond_conf=5, 
                                     initial_MLE_para=c("classical", "same"), max_steps=1e4, 
                                     obs_weights=NULL, ill_defined_value=-10^6,
                                     hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                                     method_prof=c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                                     parallel_strat=c("none", "multisession", "sequential", "multicore"), n_workers=NULL, ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  method_prof <- match.arg(method_prof)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parallel_strat <- match.arg(parallel_strat)
  
  `%fun%` <- set_doFuture_strategy(strategy=parallel_strat, n_workers=n_workers)
  
  i <- NULL
  results <- foreach::foreach(i=(1:nrow(X_rlvl)), .errorhandling="stop", .combine=dplyr::bind_rows) %fun% {
    # out <- GPD_profile_loglik_curve(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, parameter=parameter, subparam_id=0, alpha=alpha, quantile_lvl=quantile_lvl, orthogonal=orthogonal,
    #                                 X=X, x_rlvl=X_rlvl[i,], scale_cols=scale_cols, shape_cols=shape_cols,
    #                                 warmstart_table=NULL, stepsize=stepsize, steps_beyond_conf=steps_beyond_conf,
    #                                 initial_MLE_para=initial_MLE_para, max_steps=max_steps,
    #                                 hessian=hessian, maxit=maxit, method=method, method_prof=method_prof, ...)
    out <- GPD_profile_CI(Y=Y, threshold=threshold, threshold_lvl=threshold_lvl, 
                          parameter=parameter, subparam_id=0, alpha=alpha, quantile_lvl=quantile_lvl, orthogonal=orthogonal,
                          X=X, x_rlvl=X_rlvl[i,], scale_cols=scale_cols, shape_cols=shape_cols,
                          warmstart_table=NULL, init_step_pos=init_step_pos, init_step_neg=init_step_neg, tol=tol, steps_beyond_conf=steps_beyond_conf,
                          initial_MLE_para=initial_MLE_para, max_steps=max_steps, obs_weights=obs_weights, ill_defined_value=ill_defined_value,
                          hessian=hessian, maxit=maxit, method=method, method_prof=method_prof, ...)
    tibble::as_tibble(as.list(c(obs=i, out$mle[parameter], ci_down=out$ci[1], ci_up=out$ci[2])))
  }
  end_doFuture_strategy()
  
  results <- results %>% dplyr::bind_cols(parameter=parameter, alpha=alpha)
  if(parameter=="quantile"){
    results$quantile_lvl <- quantile_lvl
  }
  
  return(results)
}



# GPD_paraboot_CIs_multiple <- function(...){ #TODO:?
# }




#' GPD covariate matrices
#'
#' Creates separate GPD covariate matrices from the main UI specification.
#'
#' @param X .
#' @param x_rlvl .
#' @param scale_cols .
#' @param shape_cols .
#' @param parametrization .
#' @param out_param .
#'
#' @returns The necessary internal matrices in a list.
#'
#' @keywords internal
GPD_covariate_matrices <- function(X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL,
                                   parametrization=c("classical", "orthogonal", "quantile", "endpoint"), out_param=parametrization){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "quantile", "endpoint"))
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  if(out_param=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  center_rlvl_X <- parametrization=="quantile"|parametrization=="endpoint"|out_param=="quantile"|out_param=="endpoint"
  # center_rlvl_X <- !is.null(x_rlvl)
  if(!is.null(X)){
    if(center_rlvl_X){
      X_feats <- t(t(X) - c(x_rlvl)) #scale(X, center=c(x_rlvl))#t(t(X) - c(x_rlvl))
    }else{
      X_feats <- X #scale(X)
    }
  } else {
    if(!all(is.null(scale_cols),is.null(shape_cols))){
      stop("Must provide covariate matrix X for any non-stationary specification.")
    }
    X_feats <- NULL
  }
  if(!is.null(scale_cols)){
    scamat <- cbind(rep(1, nrow(X)), X_feats[, scale_cols, drop=F])
  }else{
    scamat <- as.matrix(1)
  }
  if(!is.null(shape_cols)){
    shamat <- cbind(rep(1, nrow(X)), X_feats[, shape_cols, drop=F])
  }else{
    shamat <- as.matrix(1)
  }
  return(list(X_feats=X_feats, scale_mat=scamat, shape_mat=shamat))
}

#' Space-efficient GPD covariate matrix with indices
#'
#' CURRENTLY UNUSED
#'
#' @param X .
#' @param x_rlvl .
#' @param scale_cols .
#' @param shape_cols .
#' @param parametrization .
#' @param out_param .
#'
#' @returns The space-efficient GPD covariate matrix with indices in a named list.
#'
#' @keywords internal
GPD_covariate_matrix <- function(X=NULL, x_rlvl=NULL, scale_cols=NULL, shape_cols=NULL,
                                 parametrization=c("classical", "orthogonal", "quantile", "endpoint"), out_param=parametrization){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "quantile", "endpoint"))
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  if(out_param=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  center_rlvl_X <- parametrization=="quantile"|parametrization=="endpoint"|out_param=="quantile"|out_param=="endpoint"
  # center_rlvl_X <- !is.null(x_rlvl)
  if(!is.null(X)){
    if(center_rlvl_X){
      X_feats <- t(t(X) - c(x_rlvl)) #scale(X, center=c(x_rlvl))#t(t(X) - c(x_rlvl))
    }else{
      X_feats <- X #scale(X)
    }
    X_feats <- cbind(rep(1, nrow(X)), X_feats)
  } else {
    if(!all(is.null(scale_cols),is.null(shape_cols))){
      stop("Must provide covariate matrix X for any non-stationary specification.")
    }
    # X_feats <- NULL
    X_feats <- as.matrix(1)
  }
  indsca <- c(1,scale_cols+1)
  indsha <- c(1,shape_cols+1)
  return(list(X_feats=X_feats, indsca=indsca, indsha=indsha))
}


#' Initial GPD parameter vector defaults for profile optimization
#'
#' @param Y .
#' @param threshold .
#' @param threshold_lvl .
#' @param nbsca .
#' @param nbsha .
#' @param parametrization .
#' @param quantile_lvl .
#'
#' @returns The initial parameter values as a vector, in the correct internal format.
#'
#' @keywords internal
GPD_param_init <- function(Y, threshold=0, threshold_lvl=0, nbsca, nbsha, 
                           parametrization=c("classical", "orthogonal", "quantile", "endpoint"), 
                           quantile_lvl, obs_weights=NULL){
  parametrization <- match.arg(parametrization)
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  # Z <- Y[Y>threshold]
  exceeds <- (Y > threshold)
  Z <- c(Y[exceeds]-threshold) #TODO check if -threshold or not
  
  if(is.null(obs_weights)){
    # Unweighted mean and variance
    meanZ <- mean(Z, na.rm = TRUE)
    varZ <- stats::var(Z, na.rm = TRUE)
    
  }else{
    if(length(obs_weights)!=length(Y)){stop('obs_weights should be of the same length as Y.')}
    
    w <- obs_weights[exceeds]
    
    if(any(w<0)){stop('obs_weights should be non-negative.')}
    if(sum(w)==0){stop('exceedence obs_weights cannot be all zero.')}
    
    # normalize weights
    w <- c(w/sum(w))
    # w <- c(w/sum(w)*length(w)) # to sum to 1
    
    # Weighted mean and variance
    meanZ <- stats::weighted.mean(Z, w, na.rm=TRUE)
    varZ <- stats::weighted.mean((Z - meanZ)^2, w, na.rm=TRUE)
    # wvZ <- stats::weighted.mean((Z - wmZ)^2, w, na.rm=TRUE) * (sum(w)/(sum(w)^2 - sum(w^2))) # weighted var
    # wvZ <- sum(w * (Z - wmZ)^2, na.rm=TRUE)
    # wvZ <- sum(w * (Z - mZ)^2, na.rm=TRUE) / (1 - sum(w^2)) # weighted var
  }
  # GEV: inisc <- c(sqrt(6 * stats::var(Y))/pi, rep(0,nbsca-1))
  # GEV: inilo <- c(mean(Y) - 0.57722 * inisc[1], rep(0,nbloc-1))
  
  inisc <- c(sqrt(6 * varZ)/pi, rep(0,nbsca-1))
  in1 <- meanZ - 0.57722 * inisc[1]
  
  if(parametrization=="endpoint"){
    inish <- c(-0.1, rep(0,nbsha-1))
    iniep <- c(GPD_endpoint(threshold, inisc[1], inish[1]), rep(0,nbsca-1))
    init = c(iniep, inish)
  } else {
    inish <- c(0.1, rep(0,nbsha-1))
    if(parametrization=="quantile"){
      # GEV: inirl <- c(GPD_quantile(inilo[1], inisc[1], inish[1], p), rep(0,nbloc-1))
      inirl <- c(GPD_quantiles(quantile_lvl=quantile_lvl, threshold_lvl=threshold_lvl, threshold=threshold, scale=inisc[1], shape=inish[1]), rep(0,nbsca-1))
      init = c(inirl, inish)
    } else {
      init = c(inisc, inish)
    }
  }
  
  return(init)
}


#' GPD profile parameter vector id
#'
#' @param parameter .
#' @param nbsca .
#' @param nbsha .
#' @param subparam_id .
#'
#' @returns The id of the desired GPD parameter in the specified internal vector format.
#'
#' @keywords internal
GPD_profpar_id <- function(parameter=c("shape", "scale", "quantile", "endpoint"),
                           nbsca=1, nbsha=1, subparam_id=0){
  parameter <- match.arg(parameter)
  if((parameter=="scale" & subparam_id>=nbsca)|(parameter=="shape" & subparam_id>=nbsha)|
     (parameter=="quantile" & subparam_id!=0)|(parameter=="endpoint" & subparam_id!=0)|(subparam_id<0)){
    stop("'subparam_id' is out of bounds.")
  }
  id <- switch(parameter,
               "shape" = nbsca+1+subparam_id,
               "scale" = 1+subparam_id,
               "quantile" = 1,
               "endpoint" = 1)
  return(id)
}


#' GPD parameter vector names
#'
#' @param nbsca .
#' @param nbsha .
#' @param parametrization .
#'
#' @returns The names of the GPD parameters in the specified internal vector format, as a vector.
#'
#' @keywords internal
GPD_param_names <- function(nbsca=1, nbsha=1, parametrization=c("classical", "orthogonal", "quantile", "endpoint")){
  parametrization <- match.arg(parametrization)
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  scanames <- if(nbsca==1) "scale" else paste0("scale", 0:(nbsca-1))
  shanames <- if(nbsha==1) "shape" else paste0("shape", 0:(nbsha-1))
  if(parametrization=="quantile" | parametrization=="endpoint"){
    scanames[1] <- parametrization
  }
  return(c(scanames, shanames))
}


#' GPD parameter vector reparametrization
#'
#' @param parameters Vector of GPD parameters in the internal format.
#' @param threshold GPD threshold value.
#' @param threshold_lvl Probability level of the threshold `threshold`.
#' @param parametrization Current parametrization of `parameters`.
#' @param new_parametrization Desired new parametrization.
#' @param quantile_lvl Quantile probability level for the `'quantile'` parameter.
#' @param nbsca Number of scale parameter coefficients (i.e. one plus the number of scale parameter covariates).
#' @param nbsha Number of shape parameter coefficients (i.e. one plus the number of shape parameter covariates).
#'
#' @returns The vector of GPD parameters reparametrized from `parametrization` to `new_parametrization`.
#' @export
GPD_change_parametrization <- function(parameters, threshold=0, threshold_lvl=0, parametrization=c("classical", "orthogonal", "quantile", "endpoint"),
                                       new_parametrization=c("classical", "orthogonal", "quantile", "endpoint"), quantile_lvl=1.-(1./100.),
                                       nbsca=1, nbsha=1){
  parametrization <- match.arg(parametrization)
  new_parametrization <- match.arg(new_parametrization)
  
  if(parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  if(new_parametrization=='orthogonal'){stop('orthogonal parametrization not implemented yet.')}
  
  sc_pars <- parameters[1:nbsca]
  sh_pars <- parameters[nbsca+(1:nbsha)]
  
  
  if(parametrization=="classical"){
    if(new_parametrization=="quantile"){
      # GEV: lo_pars[1] <- lo_pars[1] - sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
      sc_pars[1] <- threshold + sc_pars[1]/sh_pars[1] * (((1-threshold_lvl)/(1-quantile_lvl))^{sh_pars[1]} - 1)
    }
    if(new_parametrization=="endpoint"){
      # GEV: lo_pars[1] <- lo_pars[1] - sc_pars[1]/sh_pars[1]
      sc_pars[1] <- threshold - sc_pars[1]/sh_pars[1]
    }
  }
  if(parametrization=="quantile"){
    if(new_parametrization=="classical"){
      # GEV: lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
      sc_pars[1] <- ((sc_pars[1] - threshold)*sh_pars[1])/(((1-threshold_lvl)/(1-quantile_lvl))^{sh_pars[1]} - 1)
    }
    if(new_parametrization=="endpoint"){
      # GEV: lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (-(-log(1-p))^(-sh_pars[1]))
      sc_pars[1] <- threshold - (((sc_pars[1] - threshold)*sh_pars[1])/(((1-threshold_lvl)/(1-quantile_lvl))^{sh_pars[1]} - 1))/sh_pars[1]
    }
  }
  if(parametrization=="endpoint"){
    if(new_parametrization=="classical"){
      # GEV: lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1]
      sc_pars[1] <- (threshold - sc_pars[1])*sh_pars[1]
    }
    if(new_parametrization=="quantile"){
      # GEV: lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (-log(1-p))^(-sh_pars[1])
      sc_pars[1] <- threshold + ((threshold - sc_pars[1])*sh_pars[1])/sh_pars[1] * (((1-threshold_lvl)/(1-quantile_lvl))^{sh_pars[1]} - 1)
    }
  }
  new_pars <- c(sc_pars, sh_pars)
  names(new_pars) <- GPD_param_names(nbsca=nbsca, nbsha=nbsha, parametrization=new_parametrization)
  return(new_pars)
}

