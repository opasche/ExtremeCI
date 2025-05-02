
#' Insert value in vector
#'
#' @param vect A 1-D vector.
#' @param val A value to insert in the vector.
#' @param ind The index at which to insert the value in the vector,
#' must be an integer between `1` and `length(vect) + 1`.
#'
#' @return A 1-D vector of length `length(vect) + 1`,
#' with `val` inserted at position `ind` in the original `vect`.
#'
#' @keywords internal
vector_insert <- function(vect, val, ind){
  # @examples vector_insert(c(2, 7, 3, 8), val=5, ind=3)
  n <- length(vect)
  if(ind<1 | ind>(n+1)){
    stop("In 'vector_insert': 'ind' must be an integer between 1 and (length(vect) + 1).")
  }
  if(ind == 1){
    return(c(val, vect))
  }
  if(ind == (n+1)){
    return(c(vect, val))
  }
  return(c(vect[1:(ind-1)], val, vect[ind:n]))
}


#' Get doFuture operator
#'
#' @param strategy One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#'
#' @return Returns the appropriate operator to use in a [foreach::foreach()] loop.
#' The \code{\link[foreach]{\%do\%}} operator is returned if `strategy=="sequential"`.
#' Otherwise, the \code{\link[foreach]{\%dopar\%}} operator is returned.
#' @importFrom foreach %do% %dopar%
#'
#' @keywords internal
get_doFuture_operator <- function(strategy=c("sequential", "multisession", "multicore", "mixed")){
  # @examples `%fun%` <- get_doFuture_operator("sequential")
  
  strategy <- match.arg(strategy)
  
  if(strategy == "sequential"){
    return(foreach::`%do%`)
  } else {
    return(foreach::`%dopar%`)
  }
}


#' Set a doFuture execution strategy
#'
#' @param strategy One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures
#' that can be active at the same time before blocking.
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers.
#' The function should return a numeric scalar.
#' Defaults to [future::availableCores()]`-1` if `NULL` (default), with `"multicore"` constraint in the relevant case.
#' Ignored if `strategy=="sequential"`.
#'
#' @return The corresponding [get_doFuture_operator()] operator to use in a [foreach::foreach()] loop.
#' @importFrom foreach %do% %dopar%
#' @importFrom future availableCores plan sequential multisession multicore tweak
#' @importFrom doFuture registerDoFuture
#'
#' @keywords internal
set_doFuture_strategy <- function(strategy=c("sequential", "multisession", "multicore", "mixed"),
                                  n_workers=NULL){
  # @examples \dontrun{
  # `%fun%` <- set_doFuture_strategy("multisession", n_workers=3)
  # # perform foreach::foreach loop
  # end_doFuture_strategy()
  # }
  strategy <- match.arg(strategy)
  
  doFuture::registerDoFuture()
  if(strategy == "sequential"){
    future::plan(future::sequential)
    
  } else if (strategy == "multisession"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    future::plan(future::multisession, workers = n_workers)
    
  } else if (strategy == "multicore"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores(constraints = "multicore") - 1, 1)
    }
    future::plan(future::multicore, workers = n_workers)
    
  } else if (strategy == "mixed"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    strategy_1 <- future::tweak(future::sequential)
    strategy_2 <- future::tweak(future::multisession, workers = n_workers)
    future::plan(list(strategy_1, strategy_2))
  }
  return(get_doFuture_operator(strategy))
}


#' End the currently set doFuture strategy
#'
#' @description Resets the default strategy using `future::plan("default")`.
#'
#' @importFrom future plan
#'
#' @keywords internal
end_doFuture_strategy <- function(){
  # @examples \dontrun{
  # `%fun%` <- set_doFuture_strategy("multisession", n_workers=3)
  # # perform foreach::foreach loop
  # end_doFuture_strategy()
  # }
  
  future::plan("default")
}


