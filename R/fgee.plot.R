#' Plot coefficient estimates from a fitted `fastFGEE` model
#'
#' Produces coefficient plots with pointwise and joint confidence intervals
#' when available.
#'
#' @param fit A fitted object returned by \code{\link{fgee}}.
#' @param num_row Number of rows used when arranging plots.
#' @param xlab X-axis label.
#' @param title_names Optional replacement titles for coefficient plots.
#' @param ylim Optional y-axis limits.
#' @param align_x Optional value used to re-center the x-axis.
#' @param x_rescale Optional x-axis rescaling factor.
#' @param y_val_lim Expansion factor for the upper y-axis limit.
#' @param y_scal_orig Expansion factor for the lower y-axis limit.
#' @param return Logical; if `TRUE`, return the plotting data instead of only
#'   drawing the plots.
#' @param terms.plot Logical; retained for backwards compatibility.
#' @param all.terms Logical; retained for backwards compatibility.
#' @param int.uncertainty Logical; retained for backwards compatibility.
#'
#' @return Invisibly returns a list of plotting data frames when
#'   `return = TRUE`; otherwise draws plots.
#'
#' @export
fgee.plot <- function(fit,
                      num_row = NULL,
                      xlab = "Functional Domain",
                      title_names = NULL,
                      ylim = NULL,
                      align_x = NULL,
                      x_rescale = 1,
                      y_val_lim = 1.1,
                      y_scal_orig = 0.05,
                      return = FALSE,
                      terms.plot = TRUE,
                      all.terms = TRUE,
                      int.uncertainty = FALSE
){

  # Load required package
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Install it with: install.packages('gridExtra')")
  }

  # ================================================================
  # Extract coefficient estimates and CIs from new fgee structure
  # ================================================================

  # Extract functional coefficients
  if (!is.null(fit$crit) && !is.null(fit$crit$ci) && !is.null(fit$crit$ci$fit)) {
    fit$betaHat <- do.call(rbind, fit$crit$ci$fit)
  } else if (!is.null(fit$betaHat)) {
    # Already exists
    fit$betaHat <- fit$betaHat
  } else {
    stop("Cannot find functional coefficients in fit object. Expected fit$crit$ci$fit or fit$betaHat.")
  }

  num_var <- p <- nrow(fit$betaHat)  # number of variables to plot
  L <- ncol(fit$betaHat)              # number of functional domain points

  plot_list <- res_list <- vector(length = num_var, "list")
  if(is.null(num_row)) num_row <- ceiling(num_var/2)

  fit$argvals <- 1:L

  # Get variable names
  if (!is.null(rownames(fit$betaHat))) {
    var.names <- rownames(fit$betaHat)
  } else if (!is.null(fit$xnames)) {
    var.names <- fit$xnames
  } else {
    var.names <- paste0("beta", 0:(p-1))
  }

  name = NULL

  align <- ifelse(is.null(align_x), 0, align_x * x_rescale)

  # ================================================================
  # Extract standard errors and construct CIs
  # ================================================================

  # Check if CIs are available
  has_pointwise_ci <- !is.null(fit$crit) && !is.null(fit$crit$ci) &&
    !is.null(fit$crit$ci$ci_pointwise)
  has_joint_ci <- !is.null(fit$crit) && !is.null(fit$crit$ci) &&
    !is.null(fit$crit$ci$ci_joint)

  if (has_pointwise_ci || has_joint_ci) {
    # Extract pointwise CIs
    lower.pt <- matrix(NA, nrow = p, ncol = L)
    upper.pt <- matrix(NA, nrow = p, ncol = L)
    lower.joint <- matrix(NA, nrow = p, ncol = L)
    upper.joint <- matrix(NA, nrow = p, ncol = L)

    for(r in 1:p) {
      if (has_pointwise_ci) {
        lower.pt[r,] <- fit$crit$ci$ci_pointwise[[r]][, "lower"]
        upper.pt[r,] <- fit$crit$ci$ci_pointwise[[r]][, "upper"]
      }

      if (has_joint_ci) {
        lower.joint[r,] <- fit$crit$ci$ci_joint[[r]][, "lower"]
        upper.joint[r,] <- fit$crit$ci$ci_joint[[r]][, "upper"]
      }
    }

    fit$betaHat.se <- (upper.pt - lower.pt) / (2 * 1.96)  # Back-calculate SE from CI
    fit$has_ci <- TRUE
    fit$lower.pt <- lower.pt
    fit$upper.pt <- upper.pt
    fit$lower.joint <- lower.joint
    fit$upper.joint <- upper.joint

  } else {
    fit$betaHat.se <- NULL
    fit$has_ci <- FALSE
  }

  # ================================================================
  # Set title names
  # ================================================================
  if(is.null(title_names)) {
    title_names <- var.names
  }
  if(length(title_names) != p) {
    title_names <- var.names
  }
  names(res_list) <- var.names

  # ================================================================
  # Create plots for each coefficient
  # ================================================================
  for(r in 1:num_var){

    if(!fit$has_ci){
      # Plot without CIs
      beta.hat.plt <- data.frame(s = fit$argvals,
                                 beta = fit$betaHat[r,])
      plot_list[[r]] <- ggplot() +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        geom_line(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale, y = beta, color = "Estimate"),
                  data = beta.hat.plt, alpha = 1, linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_colour_manual(name="", values=c("Estimate"="black")) +
        labs(x = xlab, y = bquote(paste(beta[.(r-1)], "(s)")),
             title = title_names[r]) +
        theme(legend.position = "none")

    } else {
      # Plot with CIs
      beta.hat.plt <- data.frame(
        s = fit$argvals,
        beta = fit$betaHat[r,],
        lower = fit$lower.pt[r,],
        upper = fit$upper.pt[r,],
        lower.joint = fit$lower.joint[r,],
        upper.joint = fit$upper.joint[r,]
      )

      plot_list[[r]] <- ggplot() +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        geom_ribbon(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale,
                        ymax = upper.joint, ymin = lower.joint),
                    data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
        geom_ribbon(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale,
                        ymax = upper, ymin = lower),
                    data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
        geom_line(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale,
                      y = beta, color = "Estimate"),
                  data = beta.hat.plt, alpha = 1, linewidth = 1) +
        scale_colour_manual(name="", values=c("Estimate"="black")) +
        labs(x = xlab, y = bquote(paste(beta[.(r-1)], "(s)")),
             title = title_names[r]) +
        theme(legend.position = "none")

    }

    # Set y-axis limits
    if(!is.null(ylim)){
      plot_list[[r]] <- plot_list[[r]] + coord_cartesian(ylim = ylim)
      ylimit <- ylim
    } else {
      if(!fit$has_ci){
        ylimit <- c(min(beta.hat.plt$beta), max(beta.hat.plt$beta))
        y_adjust <- y_scal_orig * (max(beta.hat.plt$beta) - min(beta.hat.plt$beta))
      } else {
        ylimit <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint))
        y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint))
      }
      ylimit[1] <- ylimit[1] - y_adjust
    }

    xlim <- ggplot2::layer_scales(plot_list[[r]])$x$range$range

    x_range <- diff(xlim) * 0.1
    y_range <- diff(ylimit) * 0.1
    y_range_up <- diff(ylimit) * 0.02

    # Extend upper limit
    ylimit.max <- max(ylimit)
    if(ylimit.max < 0) y_val_lim <- 1 / y_val_lim
    y_val_lim_vec <- c(1, y_val_lim)
    y_top <- (0.975) * diff(ylimit * y_val_lim_vec) + ylimit[1] * y_val_lim_vec[1]

    plot_list[[r]] <- plot_list[[r]] +
      coord_cartesian(ylim = ylimit * y_val_lim_vec, xlim = xlim)

    # Add vertical line at align_x if specified
    if(!is.null(align_x)){
      plot_list[[r]] <- plot_list[[r]] +
        geom_segment(aes(y = ylimit[1] - y_range, yend = y_top,
                         x = 0, xend = 0), inherit.aes = TRUE,
                     color = "black", lwd = 0.5, alpha = 0.75, linetype = "dashed")
    }

    # Add horizontal line at 0 if CI crosses zero
    if(fit$has_ci){
      if(max(beta.hat.plt$upper.joint) > 0 & min(beta.hat.plt$lower.joint) < 0){
        plot_list[[r]] <- plot_list[[r]] +
          geom_segment(aes(x = xlim[1] - x_range, xend = xlim[2] + x_range,
                           y = 0, yend = 0), inherit.aes = TRUE,
                       color = "black", lwd = 0.5, alpha = 0.75, linetype = "dashed")
      }
      colnames(beta.hat.plt) <- c("s", "beta.hat", "CI.lower.pointwise",
                                  "CI.upper.pointwise", "CI.lower.joint", "CI.upper.joint")
    } else {
      colnames(beta.hat.plt) <- c("s", "beta.hat")
    }

    res_list[[r]] <- beta.hat.plt
  }

  # Arrange and display plots using gridExtra
  plot_return <- do.call(gridExtra::grid.arrange, c(plot_list, nrow = num_row))

  if(return == TRUE){
    res_list$plot <- plot_return
    return(invisible(res_list))
  } else {
    return(invisible(plot_return))
  }
}
