
#' Profile likelihood curve plot
#'
#' Plot the profile likelihood curve for the desired parameter with the confidence line.
#'
#' @param profile_fct_object Object returned by either [GPD_profile_loglik_curve()] or [GEV_profile_loglik_curve()].
#' @param prop_below (Optional) Proportion of the distance below the confidence line to show in the plot, to crop y-axis.
#' @param legend.position Position of the legend to the side of the plot.
#' Can be one of `"bottom"` (default), `"right"`, `"top"`, `"left"`, or `"none"`.
#'
#' @returns A [ggplot2::ggplot()] object showing the profile likelihood curve for
#' the desired parameter with the confidence line.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point geom_vline labs scale_color_manual scale_linetype_manual scale_x_continuous scale_y_continuous theme coord_cartesian
#' @importFrom rlang .data
#' @export
plot_profile_loglik_curve <- function(profile_fct_object, prop_below=NULL, legend.position="bottom"){
  
  df <- profile_fct_object$eval_table
  param_name <- profile_fct_object$param_name
  mle_par <- profile_fct_object$mle[[param_name]]
  mle_pll <- df$log_lik[df[[param_name]]==mle_par]
  conf_line <- profile_fct_object$conf_line
  CI <- profile_fct_object$ci
  
  df$conf_line <- conf_line
  df["Profile log-likelihood"] <- df$log_lik
  df_mle <- data.frame(x=mle_par, LogLik=mle_pll)
  
  # qo_plot <- df %>% tidyr::gather(key="Curve", value="LogLik", "conf_line", "Profile log-likelihood", factor_key=TRUE) %>%
  #   ggplot2::ggplot(ggplot2::aes_string(x=param_name, y="LogLik")) +
  #   ggplot2::geom_line(ggplot2::aes(color=Curve, linetype=Curve), size=1) +
  #   ggplot2::geom_point(data=df_mle, ggplot2::aes(x=x, y=LogLik), color="darkred", size=3, shape=18, inherit.aes=FALSE) +
  #   ggplot2::geom_vline(xintercept=c(mle_par, CI), linetype=c("dotdash","dotted","dotted"), color="darkgrey", size=1) +
  #   ggplot2::labs(y="Profile log-likelihood", x=param_name, color=NULL, linetype=NULL) +
  #   ggplot2::scale_color_manual(values=c("Profile log-likelihood"="blue", "conf_line"="orange"), na.value="black", guide="legend") +
  #   ggplot2::scale_linetype_manual(values=c("Profile log-likelihood"="solid", "conf_line"="dashed"), na.value="solid", guide="legend") +
  #   ggplot2::scale_x_continuous(expand=c(0.02,0)) + ggplot2::scale_y_continuous(expand=c(0.05,0)) +
  #   ggplot2::theme(legend.position=legend.position) # panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5),
  
  qo_plot <- df %>% tidyr::gather(key="Curve", value="LogLik", "conf_line", "Profile log-likelihood", factor_key=TRUE) %>%
    ggplot2::ggplot(ggplot2::aes(x=.data[['param_name']], y=.data[['LogLik']])) +
    ggplot2::geom_line(ggplot2::aes(color=.data[['Curve']], linetype=.data[['Curve']]), size=1) +
    ggplot2::geom_point(data=df_mle, ggplot2::aes(x=.data[['x']], y=.data[['LogLik']]), color="darkred", size=3, shape=18, inherit.aes=FALSE) +
    ggplot2::geom_vline(xintercept=c(mle_par, CI), linetype=c("dotdash","dotted","dotted"), color="darkgrey", size=1) +
    ggplot2::labs(y="Profile log-likelihood", x=param_name, color=NULL, linetype=NULL) +
    ggplot2::scale_color_manual(values=c("Profile log-likelihood"="blue", "conf_line"="orange"), na.value="black", guide="legend") +
    ggplot2::scale_linetype_manual(values=c("Profile log-likelihood"="solid", "conf_line"="dashed"), na.value="solid", guide="legend") +
    ggplot2::scale_x_continuous(expand=c(0.02,0)) + ggplot2::scale_y_continuous(expand=c(0.05,0)) +
    ggplot2::theme(legend.position=legend.position) # panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5),
  if(!is.null(prop_below)){
    qo_plot <- qo_plot + ggplot2::coord_cartesian(ylim=c(conf_line-prop_below*(max(df$log_lik)-conf_line), max(df$log_lik)))
  }
  
  return(qo_plot)
}



#' Non-stationary confidence bands plot
#'
#' Plot the estimated conditional quantile, return level, or endpoint, with confidence bands for non-stationary models.
#'
#' @param quantiles Vector of estimated quantiles, return levels, or endpoints.
#' @param q_down Vector of lower confidence band values for `quantiles`.
#' @param q_up Vector of upper confidence band values for `quantiles`.
#' @param time_index Vector of time indices corresponding to the `quantiles`. Defaults to `seq_along(quantiles)`.
#' @param x_label Label for the x-axis. Defaults to `"X"`.
#' @param y_label Label for the y-axis. Defaults to `"Y"`.
#' @param Y (Optional) Vector of observations to add to the plot and compare to the `quantiles.`
#' @param event_index (Optional) Index of an event (observation) to be highlighted.
#' @param obs_label Label for the fitting observations. Defaults to `"Fitting obs."`.
#' @param event_label (Optional) Label for the highlighted event observation. Defaults to `"Event"`.
#' @param legend.position Position of the legend to the side of the plot.
#' Can be one of `"bottom"` (default), `"right"`, `"top"`, `"left"`, or `"none"`.
#'
#' @returns A [ggplot2::ggplot()] object showing the estimated conditional quantile with confidence bands.
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs scale_x_continuous scale_color_manual scale_shape_manual scale_size_manual scale_y_continuous theme
#' @importFrom rlang .data
#' @export
plot_data_quantile_ci <- function(quantiles, q_down, q_up, time_index=seq_along(quantiles),
                                  x_label="X", y_label="Y", # x_label="Year", y_label="Yearly temperature maxima",
                                  Y=NULL, event_index=NULL, obs_label="Fitting obs.", event_label="Event",
                                  legend.position="bottom"){
  if(is.null(Y)){
    Y <- rep(NA, length(quantiles))
  }
  
  df <- data.frame(id=time_index, Y=Y, quantiles=quantiles, q_up=q_up, q_down=q_down, Obs_type=rep(obs_label, length(Y)))
  if(!is.null(event_index)){
    df[which(time_index==event_index),]$Obs_type <- event_label
  }
  
  preds_df <- data.frame(id=time_index, Prediction=quantiles)
  
  qo_plot <- df %>% ggplot2::ggplot(ggplot2::aes(x=.data[['id']])) +
    ggplot2::geom_line(ggplot2::aes(y=.data[['quantiles']]), color="orange", size=1.5) + # , inherit.aes=FALSE
    ggplot2::geom_ribbon(ggplot2::aes(y=.data[['quantiles']], ymin=.data[['q_down']], ymax=.data[['q_up']]), fill="orange", color="orange", alpha=0.4) +
    ggplot2::geom_point(ggplot2::aes(y=.data[['Y']], shape=.data[['Obs_type']], size=.data[['Obs_type']], color=.data[['Obs_type']])) +
    ggplot2::labs(y=y_label, x=x_label, color=NULL, linetype=NULL, shape=NULL, size=NULL) +
    ggplot2::scale_x_continuous(expand=c(0.02,0)) +
    ggplot2::scale_color_manual(values=c("Fitting obs."="black", "Event"="darkred"), na.value="black", guide="legend") +
    ggplot2::scale_shape_manual(values=c("Fitting obs."=19, "Event"=18), na.value=19, guide="legend") +
    ggplot2::scale_size_manual(values=c("Fitting obs."=2, "Event"=3), na.value=2, guide="legend") +
    ggplot2::scale_y_continuous(expand=c(0.05,0)) +
    ggplot2::theme(legend.position=legend.position) # panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5),
  
  # qo_plot <- df %>% ggplot2::ggplot(ggplot2::aes(x=id)) +
  #   ggplot2::geom_line(ggplot2::aes(y=quantiles), color="orange", size=1.5) + # , inherit.aes=FALSE
  #   ggplot2::geom_ribbon(ggplot2::aes(y=quantiles, ymin=q_down, ymax=q_up), fill="orange", color="orange", alpha=0.4) +
  #   ggplot2::geom_point(ggplot2::aes(y=Y, shape=Obs_type, size=Obs_type, color=Obs_type)) +
  #   ggplot2::labs(y=y_label, x=x_label, color=NULL, linetype=NULL, shape=NULL, size=NULL) +
  #   ggplot2::scale_x_continuous(expand=c(0.02,0)) +
  #   ggplot2::scale_color_manual(values=c("Fitting obs."="black", "Event"="darkred"), na.value="black", guide="legend") +
  #   ggplot2::scale_shape_manual(values=c("Fitting obs."=19, "Event"=18), na.value=19, guide="legend") +
  #   ggplot2::scale_size_manual(values=c("Fitting obs."=2, "Event"=3), na.value=2, guide="legend") +
  #   ggplot2::scale_y_continuous(expand=c(0.05,0)) +
  #   ggplot2::theme(legend.position=legend.position) # panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5),
  
  return(qo_plot)
}



