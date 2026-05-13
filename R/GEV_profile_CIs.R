

#' Compute return level from GEV parameters
#'
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param return_period Return period for the desired return level.
#'
#' @returns The return level of the specified GEV distribution with a return period of `return_period`.
#' In other terms, the quantile of the GEV distribution at probability level `1 - 1/return_period`.
#' @export
#'
#' @examples GEV_return_level(loc=0, scale=1, shape=0.1, return_period=100)
GEV_return_level <- function(loc=0, scale=1, shape, return_period){
  p <- 1/return_period
  # evd::qgev(1-p, loc=loc, scale=scale, shape=shape)
  return(loc - scale/shape * (1-(-log(1-p))^(-shape)))
}

#' GEV endpoint
#'
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Other optional parameters for the internal endpoint function (Unused).
#' @param endpoint_type Whether to compute the upper endpoint (returns `Inf` if `shape>=0`),
#' or the unrestricted endpoint (returns the lower GEV endpoint if `shape>=0`).
#'
#' @returns The endpoint of the specified GEV distribution.
#' @export
#'
#' @examples GEV_endpoint(loc=0, scale=1, shape=-0.1, endpoint_type='upper')
GEV_endpoint <- function(loc=0, scale=1, shape, ..., endpoint_type=c("unrestricted", "upper", "lower")){
  endpoint_type <- match.arg(endpoint_type)
  if(endpoint_type=="upper"){endpt_fct <- GEV_upper_endpoint}
  if(endpoint_type=="lower"){endpt_fct <- GEV_lower_endpoint}
  if(endpoint_type=="unrestricted"){endpt_fct <- GEV_unrestricted_endpoint}
  return(endpt_fct(loc=loc, scale=scale, shape=shape, ...))
}

#' Unrestricted GEV endpoint
#'
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Unused.
#'
#' @returns The unrestricted endpoint of the specified GEV distribution
#' (returns the upper GEV endpoint if `shape<0`, and the lower endpoint otherwise).
#'
#' @keywords internal
GEV_unrestricted_endpoint <- function(loc=0, scale=1, shape, ...){
  # @examples GEV_unrestricted_endpoint(loc=0, scale=1, shape=-0.1)
  return(loc - scale/shape)
}

#' GEV upper endpoint
#'
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Unused.
#'
#' @returns The upper endpoint of the specified GEV distribution (returns `Inf` if `shape>=0`).
#'
#' @keywords internal
GEV_upper_endpoint <- function(loc=0, scale=1, shape, ...){
  # @examples GEV_upper_endpoint(loc=0, scale=1, shape=-0.1)
  if(length(shape==1) & (length(loc)>1 | length(scale)>1)){
    shape <- rep(shape, max(length(loc),length(scale)))
  }
  endpts <- rep(Inf, length(shape))
  endpts[shape<0] <- (loc - scale/shape)[shape<0]
  return(endpts)
}

#' GEV lower endpoint
#'
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#' @param ... Unused.
#'
#' @returns The lower endpoint of the specified GEV distribution (returns `Inf` if `shape>=0`).
#'
#' @keywords internal
GEV_lower_endpoint <- function(loc=0, scale=1, shape, ...){
  # @examples GEV_lower_endpoint(loc=0, scale=1, shape=0.1)
  if(length(shape==1) & (length(loc)>1 | length(scale)>1)){
    shape <- rep(shape, max(length(loc),length(scale)))
  }
  endpts <- rep(-Inf, length(shape))
  endpts[shape>0] <- (loc - scale/shape)[shape>0]
  return(endpts)
}

#' GEV log likelihood
#'
#' @param Z Block maxima observations.
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape parameter.
#'
#' @returns The GEV log-likelihood evaluated at the given parameters, given the data observations `Y`.
#' @export
GEV_log_likelihood <- function(Z, loc, scale, shape){
  
  resc <- 1+shape*(Z-loc)/scale
  
  if(any(resc <= 0) | any(scale <= 0)){return(-10^6)} #TODO: -Inf or -1e16 -1e15 ?
  
  if(length(scale)==1 & length(shape)==1){
    return(-length(Z)*log(scale) - (1+1/shape)*sum(log(resc)) - sum(resc^(-1/shape)))
  }else if(length(scale)>1){
    return(-sum(log(scale)) - sum(log(resc)*(1/shape + 1)) - sum(resc^(-1/shape)))
  }else{
    return(-length(Z)*log(scale) - sum(log(resc)*(1/shape + 1)) - sum(resc^(-1/shape)))
  }
}

# GEV_log_likelihood_rlvl <- function(Z, return_lvl, scale, shape, return_period){
#   p <- 1/return_period
#   loc <- return_lvl + scale/shape * (1-(-log(1-p))^(-shape))
#   return(GEV_log_likelihood(Z, loc, scale, shape)) #TODO: simplify?
# }
#
# GEV_log_likelihood_endpt <- function(Z, endpoint, scale, shape){
#   loc <- endpoint + scale/shape
#   return(GEV_log_likelihood(Z, loc, scale, shape))
# }

#' Internal GEV log-likelihood function format for optimisation
#'
#' @param a .
#' @param Z .
#' @param parametrization .
#' @param return_period .
#' @param locmat .
#' @param scamat .
#' @param shamat .
#' @param negative .
#'
#' @keywords internal
GEV_log_likelihood_optim <- function(a, Z, parametrization=c("classical", "return_level", "endpoint"), return_period=100., # , "othogonal"
                                     locmat=as.matrix(1), scamat=as.matrix(1), shamat=as.matrix(1), negative=FALSE) {
  parametrization <- match.arg(parametrization)
  p <- 1/return_period
  
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  lo_pars <- a[1:nbloc]
  sc_pars <- a[nbloc+(1:nbsca)]
  sh_pars <- a[nbloc+nbsca+(1:nbsha)]
  
  if(parametrization=="return_level"){
    lo_pars[1] <- a[1] + sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
  }
  if(parametrization=="endpoint"){
    lo_pars[1] <- a[1] + sc_pars[1]/sh_pars[1]
  }
  
  location <- c(locmat %*% lo_pars)
  scale <- c(scamat %*% sc_pars)
  shape <- c(shamat %*% sh_pars)
  
  ll <- GEV_log_likelihood(Z, location, scale, shape)
  
  if(negative){
    return(-ll)
  }else{
    return(ll)
  }
}


#' Maximum-likelihood GEV estimate
#'
#' @param Z Block maxima observations.
#' @param parametrization Likelihood parametrization. Alternatives to `classical` substitute for the location parameter.
#' @param return_period Return period for the `'return_level'` parametrization.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits).
#' @param x_rlvl Covariate vector at which to reparametrize for the `'return_level'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param loc_cols Column indices of `X` to use as covariate for the (conditional) location parameter (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param out_param Additional output parametrization (same as `parametrization`, by default).
#' If `out_param != parametrization`, the parameters are reparametrized from `parametrization` to `out_param` after estimation, in a separate output.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @returns The fitted maximum-likelihood GEV as a `GEV_ML` object, containing:
#' \item{mle}{The estimated maximum likelihood GEV parameters, as a named vector (expressed in `parametrization`).}
#' \item{loglik}{The log-likelihood of the estimated parameters, given the data}
#' \item{conv}{Whether the optimisation procedure converged.
#' See the `convergence` output of [stats::optim()] for more details.}
#' \item{hessian}{The loglikelihood hessian evaluated at the estimated parameters, given the data.}
#' \item{parametrization}{Likelihood parametrization.}
#' \item{out_mle}{The estimated maximum likelihood GEV parameters, reparametrized in `out_param`.}
#' \item{out_parametrization}{Additional output parametrization.}
#' @export
GEV_maxlik <- function(Z, parametrization=c("classical", "return_level", "endpoint"), return_period=100., orthogonal=FALSE,
                       X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL, out_param=parametrization,
                       hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "return_level", "endpoint"))
  method <- match.arg(method)
  
  #Dealing with covariates
  feat_mats <- GEV_covariate_matrices(X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization, out_param=out_param)
  locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  # TODO : inner parametrization instead of out_param? Here best?
  
  # Initialization
  init <- GEV_param_init(Z=Z, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization, return_period=return_period)
  
  sol <- stats::optim(par=init, fn=GEV_log_likelihood_optim,
                      Z=Z, parametrization=parametrization, return_period=return_period,
                      locmat=locmat, scamat=scamat, shamat=shamat, negative=TRUE,
                      method = method, control = list(maxit = maxit, ...), hessian = hessian)
  mle <- sol$par
  out_mle <- GEV_change_parametrization(mle, parametrization=parametrization, new_parametrization=out_param, return_period=return_period,
                                        nbloc=nbloc, nbsca=nbsca, nbsha=nbsha)
  
  # mle names
  names(mle) <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  
  out <- list(mle=mle, loglik=-sol$value, conv=sol$convergence, hessian=-sol$hessian,
              parametrization=parametrization, out_mle=out_mle, out_parametrization=out_param)
  class(out) <- c("GEV_ML")
  
  return(out)
}


#' Optimisation evaluation step for the GEV profile likelihood
#'
#' @param a .
#' @param val .
#' @param Z .
#' @param parametrization .
#' @param id_param .
#' @param return_period .
#' @param locmat .
#' @param scamat .
#' @param shamat .
#' @param negative .
#' @param orthogonal .
#'
#' @keywords internal
optim_step_GEV_profile <- function(a, val, Z, parametrization=c("classical", "return_level", "endpoint"), id_param, return_period=100.,
                                   locmat=as.matrix(1), scamat=as.matrix(1), shamat=as.matrix(1), negative=TRUE, orthogonal=FALSE){
  parametrization <- match.arg(parametrization)
  
  pars_all <- vector_insert(a, val, id_param)
  
  pll <- GEV_log_likelihood_optim(a=pars_all, Z=Z, parametrization=parametrization, return_period=return_period,
                                  locmat=locmat, scamat=scamat, shamat=shamat, negative=negative)
  return(pll)
}


#' Internal function for the GEV profile likelihood
#'
#' @param val .
#' @param Z .
#' @param parameter .
#' @param return_period .
#' @param locmat .
#' @param scamat .
#' @param shamat .
#' @param subparam_id .
#' @param orthogonal .
#' @param init .
#' @param hessian .
#' @param maxit .
#' @param method .
#' @param ... .
#'
#' @keywords internal
GEV_profile_loglik_internal <- function(val, Z, parameter=c("shape", "location", "scale", "return_level", "endpoint"), return_period=100.,
                                        locmat=as.matrix(1), scamat=as.matrix(1), shamat=as.matrix(1), subparam_id=0,
                                        orthogonal=FALSE, init=NULL,
                                        hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
                            "return_level"="return_level", "endpoint"="endpoint")
  id_param <- GEV_profpar_id(parameter=parameter, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  
  # Initialization
  if(is.null(init)){
    init <- GEV_param_init(Z=Z, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization, return_period=return_period)
    init <- init[-id_param]
  }
  
  sol <- stats::optim(par=init, fn=optim_step_GEV_profile,
                      val=val, Z=Z, parametrization=parametrization, id_param=id_param, return_period=return_period,
                      locmat=locmat, scamat=scamat, shamat=shamat, negative=TRUE,
                      method = method, control = list(maxit = maxit, ...), hessian = hessian)
  mle_other <- sol$par
  
  # mle names
  params_names <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  names(mle_other) <- params_names[-id_param]
  
  param_val <- val
  names(param_val) <- params_names[id_param]
  
  out <- list(param_val=param_val, param_name=params_names[id_param], mle_other=mle_other,
              loglik=-sol$value, conv=sol$convergence, hessian=-sol$hessian,
              parameter=parameter, parametrization=parametrization, subparam_id=subparam_id, id_param=id_param)
  class(out) <- c("GEV_profML")
  return(out)
}


#' GEV profile log-likelihood
#'
#' @param val Parameter value at which to evaluate the GEV profile log-likelihood.
#' @param Z Block maxima observations.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param return_period Return period for the `'return_level'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'return_level'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param loc_cols Column indices of `X` to use as covariate for the (conditional) location parameter (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param init Optional initial values for the remaining parameter's optimisation process, in the correct internal format.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @returns The GEV profile log-likelihood of `parameter` evaluated at `val`, given the data,
#' as a `GEV_profML` object containing:
#' \item{param_val}{(Named) parameter value at which the GEV profile log-likelihood was evaluated.}
#' \item{param_name}{Name of the evaluated profile likelihood parameter.}
#' \item{mle_other}{Maximum-likelihood estimate of the other GEV parameters.}
#' \item{loglik}{Profile GEV log-likelihood value of `parameter` evaluated at `val`, given the data.}
#' \item{conv}{Whether the optimisation procedure converged.
#' See the `convergence` output of [stats::optim()] for more details.}
#' \item{hessian}{The loglikelihood hessian evaluated at the estimated parameters, given the data.}
#' \item{parameter}{Name of the evaluated profile likelihood parameter given as argument (redundent).}
#' \item{parametrization}{Likelihood parametrization.}
#' \item{subparam_id}{Index of the parameter coefficient for which the profile log-likelihood was computed.}
#' \item{id_param}{Index of the likelihood profile parameter, in the internal parameter vector format.}
#' @export
GEV_profile_loglik <- function(val, Z, parameter=c("shape", "location", "scale", "return_level", "endpoint"),
                               subparam_id=0, return_period=100., orthogonal=FALSE,
                               X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                               init=NULL, hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  
  #Dealing with covariates
  parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
                            "return_level"="return_level", "endpoint"="endpoint")
  feat_mats <- GEV_covariate_matrices(X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  out <- GEV_profile_loglik_internal(val=val, Z=Z, parameter=parameter, return_period=return_period,
                                     locmat=locmat, scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                     orthogonal=orthogonal, init=init,
                                     hessian=hessian, maxit=maxit, method=method, ...)
  return(out)
}


#' GEV profile log-likelihood curve
#'
#' @param Z Block maxima observations.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param alpha Confidence alpha for the profile likelihood confidence interval (i.e. for the confidence line on the profile plot).
#' @param return_period Return period for the `'return_level'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'return_level'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param loc_cols Column indices of `X` to use as covariate for the (conditional) location parameter (for conditional/non-stationary fits).
#' @param scale_cols Column indices of `X` to use as covariate for the (conditional) scale parameter (for conditional/non-stationary fits).
#' @param shape_cols Column indices of `X` to use as covariate for the (conditional) shape parameter (for conditional/non-stationary fits).
#' @param warmstart_table Evaluation table from a previous run.
#' @param stepsize Numerical size of each evaluation step, in the profile parameter's scale.
#' @param steps_beyond_conf Number of additional steps to take (in each direction)
#' after the profile log-likelihood values reach below the confidence line.
#' @param initial_MLE_para Parametrization used for the initial maximum likelihood estimate (defaults to classical, for better stability).
#' @param max_steps Maximum number of steps taken (in each direction).
#' If the confidence line was not reached, the corresponding confidence interval endpoint will be infinite.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @returns The GEV profile log-likelihood curve for the desired `parameter`,
#' with confidence line and resulting `(1-alpha)` confidence interval, as a `GEV_profileLogLik` object containing:
#' \item{mle}{The estimated maximum likelihood GEV parameters, as a named vector
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
#' \item{id_param}{Index of the profile parameter in the GEV parameter vector.}
#' @export
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer. \doi{doi:10.1007/978-1-4471-3675-0}.
GEV_profile_loglik_curve <- function(Z, parameter=c("shape", "location", "scale", "return_level", "endpoint"),
                                     subparam_id=0, alpha=0.05, return_period=100., orthogonal=FALSE,
                                     X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                                     warmstart_table=NULL, stepsize=0.1, steps_beyond_conf=5,
                                     initial_MLE_para=c("classical", "same"), max_steps=1e4, hessian=TRUE,
                                     maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
                            "return_level"="return_level", "endpoint"="endpoint")
  inner_pram <- if(initial_MLE_para=="same"){parametrization}else{initial_MLE_para}
  
  #TODO: method_prof like for GPD? not necessary here as >1 optim parameters in profile likelihood.
  
  sol_mle <- GEV_maxlik(Z=Z, parametrization=inner_pram, return_period=return_period, orthogonal=orthogonal,
                        X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
                        out_param=parametrization, hessian=hessian, maxit=maxit, method=method, ...)
  # mle and its log-likelihood
  mle <- sol_mle$out_mle
  ll_mle <- sol_mle$loglik
  
  # covariate matrices
  feat_mats <- GEV_covariate_matrices(X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  id_param <- GEV_profpar_id(parameter=parameter, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  params_names <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
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
      out <- GEV_profile_loglik_internal(val=vali, Z=Z, parameter=parameter, return_period=return_period,
                                         locmat=locmat, scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                         orthogonal=orthogonal, init=mle[-id_param],
                                         hessian=hessian, maxit=maxit, method=method, ...)
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
  class(out) <- c("GEV_profileLogLik")
  return(out)
}


#' GEV profile CI using binary search
#'
#' @param Z Block maxima observations.
#' @param parameter Parameter for which to compute the profile log-likelihood.
#' @param subparam_id Index of the parameter coefficient for which to compute the profile log-likelihood (for conditional/non-stationary fits).
#' @param alpha Confidence alpha for the profile likelihood confidence interval (i.e. for the confidence line on the profile plot).
#' @param return_period Return period for the `'return_level'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param x_rlvl Covariate vector at which to reparametrize for the `'return_level'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' @param loc_cols Column indices of `X` to use as covariate for the (conditional) location parameter (for conditional/non-stationary fits).
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
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param verbose Verbose level, as integer.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @returns The GEV profile log-likelihood confidence interval for the desired `parameter`,
#' with confidence line and resulting `(1-alpha)` confidence interval, as a `GEV_profileLogLik` object containing:
#' \item{mle}{The estimated maximum likelihood GEV parameters, as a named vector
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
#' \item{id_param}{Index of the profile parameter in the GEV parameter vector.}
#' @export
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer. \doi{doi:10.1007/978-1-4471-3675-0}.
GEV_profile_CI <- function(Z, parameter=c("shape", "location", "scale", "return_level", "endpoint"),
                           subparam_id=0, alpha=0.05, return_period=100., orthogonal=FALSE,
                           X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL, warmstart_table=NULL,
                           init_step_pos=100, init_step_neg=10, tol=0.01, steps_beyond_conf=5,
                           initial_MLE_para=c("classical", "same"), max_steps=1e3, hessian=TRUE,
                           maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), verbose=1, ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
                            "return_level"="return_level", "endpoint"="endpoint")
  inner_pram <- if(initial_MLE_para=="same"){parametrization}else{initial_MLE_para}
  
  #TODO: method_prof like for GPD? not necessary here as >1 optim parameters in profile likelihood.
  
  if(verbose>=5){cat("DEBUG: START MLE\n")}
  
  sol_mle <- GEV_maxlik(Z=Z, parametrization=inner_pram, return_period=return_period, orthogonal=orthogonal,
                        X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
                        out_param=parametrization, hessian=hessian, maxit=maxit, method=method, ...)
  # mle and its log-likelihood
  mle <- sol_mle$out_mle
  ll_mle <- sol_mle$loglik
  
  # covariate matrices
  feat_mats <- GEV_covariate_matrices(X=X, x_rlvl=x_rlvl, loc_cols=loc_cols, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  id_param <- GEV_profpar_id(parameter=parameter, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  params_names <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
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
      out <- GEV_profile_loglik_internal(val=vali, Z=Z, parameter=parameter, return_period=return_period,
                                         locmat=locmat, scamat=scamat, shamat=shamat, subparam_id=subparam_id,
                                         orthogonal=orthogonal, init=mle[-id_param],
                                         hessian=hessian, maxit=maxit, method=method, ...)
      plli <- out$loglik
      # TODO: catch errors and set to -Inf or check if plli is Inf/-Inf/NA/NaN/NULL
      
      if(going_down){
        eval_table <- dplyr::bind_rows(c(out$param_val, out$mle_other, log_lik=out$loglik), eval_table)
      }else{
        eval_table <- dplyr::bind_rows(eval_table, c(out$param_val, out$mle_other, log_lik=out$loglik))
      }
    }
    
    if((going_down & vali > mle[[id_param]]+1e-14) | (!going_down & vali < mle[[id_param]]-1e-14)){
      stop('In "GEV_profile_CI": CI search got back to MLE. This might be due to numerical instability, or to a large "tol" value.')
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
            stop('In "GEV_profile_CI": algorithm should be in final phase but is not.')
            # final_phase <- TRUE
          }
        }else{
          stop('In "GEV_profile_CI": stepsize got to zero.')
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
  class(out) <- c("GEV_profileLogLik")
  return(out)
}


# find_profile_ci_GEV <- function(Z, parameter=c("shape", "location", "scale", "return_level", "endpoint"),
#                                 subparam_id=0, alpha=0.05, return_period=100., orthogonal=FALSE,
#                                 X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
#                                 warmstart_table=NULL,
#                                 hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"), ...){
#   parameter <- match.arg(parameter)
#   method <- match.arg(method)
#   parametrization <- switch(parameter, "shape"="classical", "location"="classical", "scale"="classical",
#                             "return_level"="return_level", "endpoint"="endpoint")
#
#   #TODO
#   return(list(mle=mle, ci=c(ci_down, ci_up), eval_table=eval_table))
# }

#' Multi-value conditional GEV profile likelihood confidence intervals
#'
#' @description
#' For non-stationary models, the return level reparametrization depends on covariate values.
#' This function repeats the profile likelihood procedure for several covariate values.
#' It enables obtaining a return-level (or endpoint) curve, with profile-likelihood confidence bands,
#' as a function of the covariate values.
#'
#'
#' @param Z Block maxima observations.
#' @param parameter Parameter for which to compute the profile likelihood confidence intervals.
#' @param alpha Confidence alpha for the profile likelihood confidence intervals.
#' @param return_period Return period for the `'return_level'` parameter.
#' @param orthogonal DEPRECATED.
#' @param X Covariate matrix (for conditional/non-stationary fits). Columns should be variables, and rows should be observations matching `Y`.
#' @param X_rlvl Covariate matrix at which to reparametrize for the `'return_level'` or `'endpoint'` parametrizations (for conditional/non-stationary fits).
#' Columns should be variables, and each row should give one covariate realization at which to reparametrize and obtain a CI.
#' @param loc_cols Column indices of `X` to use as covariate for the (conditional) location parameter (for conditional/non-stationary fits).
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
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned? See [stats::optim()] for more details.
#' @param maxit The maximum number of iterations. See [stats::optim()] for more details.
#' @param method The optimisation method to be used. See [stats::optim()] for more details.
#' @param parallel_strat Parallel strategy. One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures
#' that can be active at the same time before blocking.
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers.
#' The function should return a numeric scalar.
#' Defaults to [future::availableCores()]`-1` if `NULL` (default), with `"multicore"` constraint in the relevant case.
#' Ignored if `strategy=="sequential"`.
#' @param ... Other arguments passed to the `control` argument of [stats::optim()].
#'
#' @returns The GEV profile log-likelihood `(1-alpha)` confidence intervals for the desired `parameter`,
#' for each desired covariate values, as a [tibble::tibble()], with columns:
#' \item{obs}{Index of the observation (i.e. row) of X_rlvl for which the parameter estimate and CI was computed.}
#' \item{`<parameter name>`}{Conditional estimate of the `parameter`.}
#' \item{ci_down}{Lower endpoint of the conditional `(1-alpha)` profile-likelihood confidence interval for `parameter`.}
#' \item{ci_up}{Upper endpoint of the conditional `(1-alpha)` profile-likelihood confidence interval for `parameter`.}
#' \item{parameter}{Name of the parameter for which the estimates and CIs were computed.}
#' \item{alpha}{Confidence alpha for the profile likelihood confidence intervals.}
#' \item{return_period}{Return period for the `'return_level'` parameter (only if `parameter==return_level`).}
#' @export
GEV_profile_CIs_multiple <- function(Z, parameter=c("return_level", "endpoint"),
                                     alpha=0.05, return_period=100., orthogonal=FALSE,
                                     X=NULL, X_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                                     init_step_pos=100, init_step_neg=10, tol=0.01, steps_beyond_conf=5,
                                     initial_MLE_para=c("classical", "same"), max_steps=1e4,
                                     hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"),
                                     parallel_strat=c("none", "multisession", "sequential", "multicore"), n_workers=NULL, ...){
  parameter <- match.arg(parameter)
  method <- match.arg(method)
  initial_MLE_para <- match.arg(initial_MLE_para)
  parallel_strat <- match.arg(parallel_strat)
  
  `%fun%` <- set_doFuture_strategy(strategy=parallel_strat, n_workers=n_workers)
  
  i <- NULL
  results <- foreach::foreach(i=(1:nrow(X_rlvl)), .errorhandling="stop", .combine=dplyr::bind_rows) %fun% {
    # out <- GEV_profile_loglik_curve(Z=Z, parameter=parameter, subparam_id=0, alpha=alpha, return_period=return_period, orthogonal=orthogonal,
    #                                 X=X, x_rlvl=X_rlvl[i,], loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
    #                                 warmstart_table=NULL, stepsize=stepsize, steps_beyond_conf=steps_beyond_conf,
    #                                 initial_MLE_para=initial_MLE_para, max_steps=max_steps,
    #                                 hessian=hessian, maxit=maxit, method=method, ...)
    out <- GEV_profile_CI(Z=Z, parameter=parameter, subparam_id=0, alpha=alpha, return_period=return_period, orthogonal=orthogonal,
                          X=X, x_rlvl=X_rlvl[i,], loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols, warmstart_table=NULL,
                          init_step_pos=init_step_pos, init_step_neg=init_step_neg, tol=tol, steps_beyond_conf=steps_beyond_conf,
                          initial_MLE_para=initial_MLE_para, max_steps=max_steps,
                          hessian=hessian, maxit=maxit, method=method, ...)
    tibble::as_tibble(as.list(c(obs=i, out$mle[parameter], ci_down=out$ci[1], ci_up=out$ci[2])))
  }
  end_doFuture_strategy()
  
  results <- results %>% dplyr::bind_cols(parameter=parameter, alpha=alpha)
  if(parameter=="return_level"){
    results$return_period <- return_period
  }
  
  return(results)
}

#' GEV_paraboot_CIs_multiple (IN DEVELOPPEMENT)
#'
#' @param Z Block maxima observations.
#' @param alpha .
#' @param return_period .
#' @param R .
#' @param endpoint_type .
#' @param orthogonal .
#' @param X .
#' @param X_rlvl .
#' @param loc_cols .
#' @param scale_cols .
#' @param shape_cols .
#' @param bootstrap_X .
#' @param hessian .
#' @param maxit .
#' @param method .
#' @param parallel_strat .
#' @param n_workers .
#' @param ... .
#'
#' @returns Parametric bootstrap CIs
#'
#' @keywords internal
GEV_paraboot_CIs_multiple <- function(Z, alpha=0.05, return_period=100., R=1000, endpoint_type=c("upper", "lower", "unrestricted"), orthogonal=FALSE,
                                      X=NULL, X_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                                      bootstrap_X=FALSE, hessian=TRUE, maxit=1e6, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"),
                                      parallel_strat=c("none", "multisession", "sequential", "multicore"), n_workers=NULL, ...){
  method <- match.arg(method)
  parallel_strat <- match.arg(parallel_strat)
  endpoint_type <- match.arg(endpoint_type)
  
  n <- length(c(Z))
  
  parametrization <- "classical"
  inner_pram <- "classical"
  
  gev_model0 <- GEV_maxlik(Z=Z, parametrization=inner_pram, return_period=return_period, orthogonal=orthogonal,
                           X=X, x_rlvl=NULL, loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
                           out_param=parametrization, hessian=hessian, maxit=maxit, method=method, ...)
  #mle=mle, loglik=-sol$value, conv=sol$convergence, hessian=-sol$hessian, parametrization=parametrization, out_mle=out_mle, out_parametrization=out_param
  
  # mle and its log-likelihood
  mle0 <- gev_model0$out_mle
  ll_mle0 <- gev_model0$loglik
  
  # covariate matrices
  feat_mats <- GEV_covariate_matrices(X=X, x_rlvl=NULL, loc_cols=loc_cols, scale_cols=scale_cols,
                                      shape_cols=shape_cols, parametrization=parametrization)
  locmat <- feat_mats$loc_mat
  scamat <- feat_mats$scale_mat
  shamat <- feat_mats$shape_mat
  nbloc <- ncol(locmat)
  nbsca <- ncol(scamat)
  nbsha <- ncol(shamat)
  
  #id_param <- GEV_profpar_id(parameter=parameter, nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, subparam_id=subparam_id)
  params_names <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=parametrization)
  #param_name <- params_names[id_param]
  
  # for each x_rlvl value, simulate GEV at that loc and git a GEV? YES and option "bootstrap_X"
  # Simulate new Z given X datasets to refit GEV and then predict all rlvls with x_rlvl each time?
  
  lo_pars <- mle0[1:nbloc]
  sc_pars <- mle0[nbloc+(1:nbsca)]
  sh_pars <- mle0[nbloc+nbsca+(1:nbsha)]
  
  location <- c(locmat %*% lo_pars)
  scale <- c(scamat %*% sc_pars)
  shape <- c(shamat %*% sh_pars)
  
  #TODO: for multiple X_rlvl...?
  rls <- GEV_return_level(loc=location, scale=scale, shape=shape, return_period=return_period)
  eps <- GEV_endpoint(loc=location, scale=scale, shape=shape, endpoint_type=endpoint_type)
  
  `%fun%` <- set_doFuture_strategy(strategy=parallel_strat, n_workers=n_workers)
  
  i <- NULL
  par_results <- foreach::foreach(i=seq(R), .errorhandling="stop", .combine=dplyr::bind_rows) %fun% {
    if(bootstrap_X){
      stop("In 'GEV_paraboot_CIs_multiple': 'bootstrap_X' not implemented yet.")
    }else{
      sample <- evd::rgev(n, location, scale, shape)
      gev_modeli <- GEV_maxlik(Z=sample, parametrization=inner_pram, return_period=return_period, orthogonal=orthogonal,
                               X=X, x_rlvl=NULL, loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
                               out_param=parametrization, hessian=hessian, maxit=maxit, method=method, ...)
      mle <- gev_modeli$out_mle
      # ll_mle <- gev_modeli$loglik
      
      lo_parsi <- mle[1:nbloc]
      sc_parsi <- mle[nbloc+(1:nbsca)]
      sh_parsi <- mle[nbloc+nbsca+(1:nbsha)]
      
      locationi <- c(locmat %*% lo_parsi)
      scalei <- c(scamat %*% sc_parsi)
      shapei <- c(shamat %*% sh_parsi)
      
      #compute RL EP TODO for multiple x_rl?
      rlsi <- GEV_return_level(loc=locationi, scale=scalei, shape=shapei, return_period=return_period)
      epsi <- GEV_endpoint(loc=locationi, scale=scalei, shape=shapei, endpoint_type=endpoint_type)
    }
    
    # out <- GEV_profile_loglik_curve(Z=Z, parameter=parameter, subparam_id=0, alpha=alpha, return_period=return_period, orthogonal=orthogonal,
    #                                 X=X, x_rlvl=X_rlvl[i,], loc_cols=loc_cols, scale_cols=scale_cols, shape_cols=shape_cols,
    #                                 warmstart_table=NULL, stepsize=stepsize, steps_beyond_conf=steps_beyond_conf,
    #                                 initial_MLE_para=initial_MLE_para, hessian=hessian, maxit=maxit, method=method, ...)
    
    names(mle) <- params_names
    # tibble::as_tibble(as.list(c(obs=i, out$mle[parameter], ci_down=out$ci[1], ci_up=out$ci[2])))
    # TODO: c(mle, rl1, rl2, ..., ep1, ep2, ...)
    mle
  }
  
  # rl_results <- foreach::foreach(i=(1:nrow(X_rlvl)), .errorhandling="stop", .combine=dplyr::bind_rows) %fun% {
  #
  # }
  
  end_doFuture_strategy()
  
  # results <- results %>% dplyr::bind_cols(parameter=parameter, alpha=alpha)
  # if(parameter=="return_level"){
  #   results$return_period <- return_period
  # }
  results <- results %>% dplyr::bind_cols(alpha=alpha)
  results$return_period <- return_period
  
  #return(results)
  return(list(mle=mle0, loglik=ll_mle0, return_level=rls, endpoint=eps, CIs=results, GEV_object=gev_model0))
}





#' GEV covariate matrices
#'
#' Creates separate GEV covariate matrices from the main UI specification.
#'
#' @param X .
#' @param x_rlvl .
#' @param loc_cols .
#' @param scale_cols .
#' @param shape_cols .
#' @param parametrization .
#' @param out_param .
#'
#' @returns The necessary internal matrices in a list.
#'
#' @keywords internal
GEV_covariate_matrices <- function(X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                                   parametrization=c("classical", "return_level", "endpoint"), out_param=parametrization){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "return_level", "endpoint"))
  center_rlvl_X <- parametrization=="return_level"|parametrization=="endpoint"|out_param=="return_level"|out_param=="endpoint"
  # center_rlvl_X <- !is.null(x_rlvl)
  if(!is.null(X)){
    if(center_rlvl_X){
      X_feats <- t(t(X) - c(x_rlvl)) #scale(X, center=c(x_rlvl))#t(t(X) - c(x_rlvl))
    }else{
      X_feats <- X #scale(X)
    }
  } else {
    if(!all(is.null(loc_cols),is.null(scale_cols),is.null(shape_cols))){
      stop("Must provide covariate matrix X for any non-stationary specification.")
    }
    X_feats <- NULL
  }
  if(!is.null(loc_cols)){ #TODO: fct
    locmat <- cbind(rep(1, nrow(X)), X_feats[, loc_cols, drop=F])
  }else{
    locmat <- as.matrix(1)
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
  return(list(X_feats=X_feats, loc_mat=locmat, scale_mat=scamat, shape_mat=shamat))
}

#' Space-efficient GEV covariate matrix with indices
#'
#' CURRENTLY UNUSED
#'
#' @param X .
#' @param x_rlvl .
#' @param loc_cols .
#' @param scale_cols .
#' @param shape_cols .
#' @param parametrization .
#' @param out_param .
#'
#' @returns The space-efficient GEV covariate matrix with indices in a named list.
#'
#' @keywords internal
GEV_covariate_matrix <- function(X=NULL, x_rlvl=NULL, loc_cols=NULL, scale_cols=NULL, shape_cols=NULL,
                                 parametrization=c("classical", "return_level", "endpoint"), out_param=parametrization){
  parametrization <- match.arg(parametrization)
  out_param <- match.arg(out_param, c("classical", "return_level", "endpoint"))
  center_rlvl_X <- parametrization=="return_level"|parametrization=="endpoint"|out_param=="return_level"|out_param=="endpoint"
  # center_rlvl_X <- !is.null(x_rlvl)
  if(!is.null(X)){
    if(center_rlvl_X){
      X_feats <- t(t(X) - c(x_rlvl)) #scale(X, center=c(x_rlvl))#t(t(X) - c(x_rlvl))
    }else{
      X_feats <- X #scale(X)
    }
    X_feats <- cbind(rep(1, nrow(X)), X_feats)
  } else {
    if(!all(is.null(loc_cols),is.null(scale_cols),is.null(shape_cols))){
      stop("Must provide covariate matrix X for any non-stationary specification.")
    }
    # X_feats <- NULL
    X_feats <- as.matrix(1)
  }
  indloc <- c(1,loc_cols+1)
  indsca <- c(1,scale_cols+1)
  indsha <- c(1,shape_cols+1)
  return(list(X_feats=X_feats, indloc=indloc, indsca=indsca, indsha=indsha))
}


#' Initial GEV parameter vector defaults for profile optimization
#'
#' @param Z .
#' @param nbloc Number of location parameter coefficients (i.e. one plus the number of shape parameter covariates).
#' @param nbsca Number of scale parameter coefficients (i.e. one plus the number of scale parameter covariates).
#' @param nbsha Number of shape parameter coefficients (i.e. one plus the number of shape parameter covariates).
#' @param parametrization .
#' @param return_period .
#'
#' @returns The initial parameter values as a vector, in the correct internal format.
#'
#' @keywords internal
GEV_param_init <- function(Z, nbloc, nbsca, nbsha, parametrization=c("classical", "return_level", "endpoint"), return_period){
  parametrization <- match.arg(parametrization)
  p <- 1/return_period
  inisc <- c(sqrt(6 * stats::var(Z))/pi, rep(0,nbsca-1))
  inilo <- c(mean(Z) - 0.57722 * inisc[1], rep(0,nbloc-1))
  if(parametrization=="endpoint"){
    inish <- c(-0.1, rep(0,nbsha-1))
    iniep <- c(GEV_endpoint(inilo[1], inisc[1], inish[1]), rep(0,nbloc-1))
    init = c(iniep, inisc, inish)
  } else {
    inish <- c(0.1, rep(0,nbsha-1))
    if(parametrization=="return_level"){
      inirl <- c(GEV_return_level(inilo[1], inisc[1], inish[1], p), rep(0,nbloc-1))
      init = c(inirl, inisc, inish)
    } else {
      init = c(inilo, inisc, inish)
    }
  }
  return(init)
}


#' GEV profile parameter vector id
#'
#' @param parameter .
#' @param nbloc Number of location parameter coefficients (i.e. one plus the number of shape parameter covariates).
#' @param nbsca Number of scale parameter coefficients (i.e. one plus the number of scale parameter covariates).
#' @param nbsha Number of shape parameter coefficients (i.e. one plus the number of shape parameter covariates).
#' @param subparam_id .
#'
#' @returns The id of the desired GEV parameter in the specified internal vector format.
#'
#' @keywords internal
GEV_profpar_id <- function(parameter=c("shape", "location", "scale", "return_level", "endpoint"),
                           nbloc=1, nbsca=1, nbsha=1, subparam_id=0){
  parameter <- match.arg(parameter)
  if((parameter=="location" & subparam_id>=nbloc)|(parameter=="scale" & subparam_id>=nbsca)|(parameter=="shape" & subparam_id>=nbsha)|
     (parameter=="return_level" & subparam_id!=0)|(parameter=="endpoint" & subparam_id!=0)|(subparam_id<0)){
    stop("'subparam_id' is out of bounds.")
  }
  id <- switch(parameter,
               "shape" = nbloc+nbsca+1+subparam_id,
               "location" = 1+subparam_id,
               "scale" = nbloc+1+subparam_id,
               "return_level" = 1,
               "endpoint" = 1)
  return(id)
}


#' GEV parameter vector names
#'
#' @param nbloc .
#' @param nbsca .
#' @param nbsha .
#' @param parametrization .
#'
#' @returns The names of the GEV parameters in the specified internal vector format, as a vector.
#'
#' @keywords internal
GEV_param_names <- function(nbloc=1, nbsca=1, nbsha=1, parametrization=c("classical", "return_level", "endpoint")){
  parametrization <- match.arg(parametrization)
  locnames <- if(nbloc==1) "location" else paste0("location", 0:(nbloc-1))
  scanames <- if(nbsca==1) "scale" else paste0("scale", 0:(nbsca-1))
  shanames <- if(nbsha==1) "shape" else paste0("shape", 0:(nbsha-1))
  if(parametrization=="return_level" | parametrization=="endpoint"){
    locnames[1] <- parametrization
  }
  return(c(locnames, scanames, shanames))
}


#' GEV parameter vector reparametrization
#'
#' @param parameters Vector of GEV parameters in the internal format.
#' @param parametrization Current parametrization of `parameters`.
#' @param new_parametrization Desired new parametrization.
#' @param return_period Return period for the `'return_level'` parameter.
#' @param nbloc Number of location parameter coefficients (i.e. one plus the number of shape parameter covariates).
#' @param nbsca Number of scale parameter coefficients (i.e. one plus the number of scale parameter covariates).
#' @param nbsha Number of shape parameter coefficients (i.e. one plus the number of shape parameter covariates).
#'
#' @returns The vector of GEV parameters reparametrized from `parametrization` to `new_parametrization`.
#' @export
GEV_change_parametrization <- function(parameters, parametrization=c("classical", "return_level", "endpoint"),
                                       new_parametrization=c("classical", "return_level", "endpoint"), return_period=100.,
                                       nbloc=1, nbsca=1, nbsha=1){
  parametrization <- match.arg(parametrization)
  new_parametrization <- match.arg(new_parametrization)
  p <- 1/return_period
  
  lo_pars <- parameters[1:nbloc]
  sc_pars <- parameters[nbloc+(1:nbsca)]
  sh_pars <- parameters[nbloc+nbsca+(1:nbsha)]
  
  
  if(parametrization=="classical"){
    if(new_parametrization=="return_level"){
      lo_pars[1] <- lo_pars[1] - sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
    }
    if(new_parametrization=="endpoint"){
      lo_pars[1] <- lo_pars[1] - sc_pars[1]/sh_pars[1]
    }
  }
  if(parametrization=="return_level"){
    if(new_parametrization=="classical"){
      lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (1-(-log(1-p))^(-sh_pars[1]))
    }
    if(new_parametrization=="endpoint"){ # TODO: double-check correct, with GPD
      lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (-(-log(1-p))^(-sh_pars[1]))
    }
  }
  if(parametrization=="endpoint"){
    if(new_parametrization=="classical"){
      lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1]
    }
    if(new_parametrization=="return_level"){
      lo_pars[1] <- lo_pars[1] + sc_pars[1]/sh_pars[1] * (-log(1-p))^(-sh_pars[1])
    }
  }
  new_pars <- c(lo_pars, sc_pars, sh_pars)
  names(new_pars) <- GEV_param_names(nbloc=nbloc, nbsca=nbsca, nbsha=nbsha, parametrization=new_parametrization)
  return(new_pars)
}

