#' Default fastFGEE plotting
#'
#' Take a fitted \code{fgee} object produced by \code{fastFGEE::fgee()},
#' plots the point estimates, as well as 95\% pointwise and joint confidence intervals.
#'
#' @param fgeeobj A object returned from the \code{fgee} function
#' @param num_row An integer that specifies the number of rows the plots will be displayed on. Defaults to p/2, where p is the number of predictors.
#' @param xlab A string that specifies the x-axis title (i.e., for the functional domain). Defaults to ``Functional Domain''
#' @param title_names A vector of strings that has length equal to number of covariates (plus intercept if relevant). Allows one to change the titles of the plots. Defaults to NULL which uses the variable names in the dataframe for titles.
#' @param ylim A 2-dimensional vector that specifies limits of the y-axis in plots.
#' @param align_x A scalar: aligns the plot to a certain point on the functional domain and sets this as 0. This is particularly useful if the functional domain is time. Defaults to 0.
#' @param x_rescale A scalar: rescales the x-axis of plots which is especially useful if time is the functional domain and one wishes to, for example, account for the sampling rate. Defaults to 1.
#' @param y_val_lim A positive scalar that extends the y-axis by a factor for visual purposes. Defaults to 1.10. Typically does not require adjustment.
#' @param y_scal_orig A positive scalar that determines how much to reduce the length of the y-axis on the bottom. Defaults to 0.05. Typically does not require adjustment.
#' @param return Logical, indicating whether to return the data frame with the coefficient estimates and 95% confidence intervals (CIs). Defaults to \code{FALSE}.
#' @param terms.plot Logical. Uses mgcv::plot.gam output to produce coefficient plots.
#' @param all.terms Logical. Defaults to \code{TRUE}. If set to TRUE then the partial effects of parametric model components are also plotted. See ?mgcv::plot.gam for more
#' @param int.uncertainty Logical. Defaults to \code{FALSE}. Uses mgcv::predict.gam() to account for uncertainty in the itercept.
#' @return Plots of point estimates and CIs. If \code{return = TRUE}, also returns
#' a list where each element is a data frame with the coefficient estimates and 95% CIs.
#'
#' @author Gabriel Loewinger \email{gloewinger@@gmail.com}
#'
#' @references Loewinger, G., Levis, L., Cui, E., Pereira, F. (2025). Fast Penalized 
#' Generalized Estimating Equations for Large Longitudinal Functional Datasets. \emph{arXiv}.
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' data("DTI", package = "refund")
#' set.seed(1)
#' DTI_use <- DTI[DTI$ID %in% sample(DTI$ID, 25),]
#' DTI_use <- data.frame(cca = I(DTI_use$cca),
#'                       case = DTI_use$case,
#'                       visit = 1*(DTI_use$visit),
#'                       sex = DTI_use$sex,
#'                       ID = DTI_use$ID)
#' 
#' fit_dti <- fastFGEE::fgee(formula = cca ~ case + visit,
#'                           cluster = "ID",
#'                           family = "gaussian", 
#'                           cov.type = "exchangeable"
#'                           data = DTI_use) # fit model
#'                           
#' fastFGEE::fgee.plot(fit_dti) # plot

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
  
  num_var <- p <- length(fit$model$sp) ## number of variables to plot
  plot_list <- res_list <- vector(length = num_var, "list")
  if(is.null(num_row))  num_row <- ceiling(num_var/2)
  L <- max(fit$model$var.summary$yindex.vec)
  fit$argvals  <- 1:L
  var.names <- names(fit$model$var.summary)[-1]
  # if(length(var.names) != p){
  #   var.names <- names(fit$model$sp)
  # }
  name = NULL
  
  align <- ifelse(is.null(align_x), 0, align_x*x_rescale)
  # ----------------------------------------------------------------- 
  # construct betahat and variance matrices
  # ----------------------------------------------------------------- 
  # dd <- data.frame(yindex.vec = 1:L)
  # for(r in 1:(p-1)){
  #   dd[, var.names[r]] <- 1
  # }   
  
  if(terms.plot){
    # if use plot()
    preds <- data.frame(matrix(NA, nrow = 100, ncol = 2))
    fit$argvals <- 1:100 # need this to make consistent
    colnames(preds) <- c("fit", "se")
    fit$Vp <- fit$vb
    pdf(file = NULL)
    model_plot <- plot.gam(fit$model, se = TRUE, seWithMean = TRUE, all.terms = all.terms)
    invisible(dev.off()) # prevent showing plots
    preds$se <- do.call(cbind, lapply(model_plot, function(mm) mm$se / 2 ) ) # mgcv appears to multiply se by 2 (maybe for CIs?)
    preds$fit <- do.call(cbind, lapply(model_plot, function(mm) mm$fit ) ) # coefs
    fit$betaHat <- t(preds$fit)
    
    # chunk below for adding uncertainty likely unecessary. yields different results for:
    # cbind(preds0$se.fit[,1], preds$se[,1])
    if(all.terms){
      preds$fit[,1] <- preds$fit[,1] + fit$model$coefficients[1] # add in constant term -- uncertainty accounted for with all.terms argument
      
      if(int.uncertainty){
        # add in uncertainty from intercept
        d0 <- data.frame( matrix(1, nrow = nrow(preds$fit), ncol = ncol(preds$fit)) )
        d0[,1] <- 1:nrow(d0)
        colnames(d0) <- c("yindex.vec", var.names)
        preds0 <- suppressWarnings( mgcv::predict.gam(fit$model, 
                                                      d0, 
                                                      type="iterms", 
                                                      se.fit=TRUE) )
        preds$se[,1] <- preds0$se.fit[,1]
      }
      
    }   
    
    fit$betaHat.se <- preds$se
    
  }else{
    # always plots parameteric component (with uncertainty) in functional intercept
    # calculates fitted values with predict.gam in a different way
    dd <- data.frame(matrix(1, nrow = L, ncol = p))
    dd[,1] <- 1:L
    colnames(dd) <- c("yindex.vec", var.names)
    
    preds <- mgcv::predict.gam(fit$model, newdata = dd, type="iterms", se.fit=TRUE)
    preds$fit[,1]=preds$fit[,1] + fit$beta[1]
    fit$betaHat <- t(preds$fit)
    
    # intercept
    d0 <- data.frame(matrix(1, nrow = L, ncol = p)) # may need to fill matrix with 0s
    d0[,1] <- 1:L
    colnames(d0) <- c("yindex.vec", var.names)
    preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
    preds$se.fit[,1] <- preds0$se.fit[,1]
    fit$betaHat.se <- preds$se.fit
  }

  # ----------------------------------------------------------------- 
  if(is.null(title_names))    title_names <- rownames(fit$betaHat)
  if(nrow(fit$betaHat) != length(title_names) )  title_names <- rownames(fit$betaHat)
  names(res_list) <- rownames(fit$betaHat)
  
  for(r in 1:num_var){
    
    if(is.null(fit$betaHat.se)){
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
      
    }else{
      
      beta.hat.plt <- data.frame(s = fit$argvals,
                                 beta = fit$betaHat[r,],
                                 lower = fit$betaHat[r,] - 1.96*fit$betaHat.se[,r],
                                 upper = fit$betaHat[r,] + 1.96*fit$betaHat.se[,r],
                                 lower.joint = fit$betaHat[r,] - fit$qn[r]*fit$betaHat.se[,r],
                                 upper.joint = fit$betaHat[r,] + fit$qn[r]*fit$betaHat.se[,r])
      
      plot_list[[r]] <- ggplot() +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        geom_ribbon(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale, ymax = upper.joint, ymin = lower.joint),
                    data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
        geom_ribbon(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale, ymax = upper, ymin = lower),
                    data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
        geom_line(aes(x = s / x_rescale - align/x_rescale - 1/x_rescale, y = beta, color = "Estimate"),
                  data = beta.hat.plt, alpha = 1, linewidth = 1) +
        scale_colour_manual(name="", values=c("Estimate"="black")) +
        labs(x = xlab, y = bquote(paste(beta[.(r-1)], "(s)")),
             title = title_names[r]) +
        theme(legend.position = "none")
      
    }
    
    # make x and y intercepts
    if(!is.null(ylim)){
      plot_list[[r]] <- plot_list[[r]] + coord_cartesian(ylim = ylim)
      ylimit <- ylim
    }else{
      if(is.null(fit$betaHat.se)){
        ylimit <- c(min(beta.hat.plt$beta), max(beta.hat.plt$beta))
        y_adjust <- y_scal_orig * (max(beta.hat.plt$beta) - min(beta.hat.plt$beta))
      }else{
        ylimit <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint)) #layer_scales(plot_list[[r]])$y$range$range
        y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint)) #layer_scales(plot_list[[r]])$y$range$range
      }
      ylimit[1] <- ylimit[1] - y_adjust # just scale bottom because top is scaled below
    }
    
    xlim <- layer_scales(plot_list[[r]])$x$range$range
    
    x_range <- diff(xlim) * 0.1
    y_range <- diff(ylimit) * 0.1
    y_range_up <- diff(ylimit) * 0.02
    
    # extend upper limit
    ylimit.max <- max(ylimit) # find largest
    if(ylimit.max < 0)  y_val_lim <- 1 / y_val_lim # if negative, need to invert
    y_val_lim_vec <- c(1, y_val_lim)
    y_top <- (0.975) * diff(ylimit*y_val_lim_vec) + ylimit[1]*y_val_lim_vec[1]
    plot_list[[r]] <- plot_list[[r]] + 
                        coord_cartesian(ylim = ylimit*y_val_lim_vec,
                                        xlim = xlim)
    
    if(!is.null(align_x)){
      
      plot_list[[r]] <- plot_list[[r]] +
        geom_segment(aes_string(y=ylimit[1] - y_range, yend=y_top,
                                x=0,xend=0), inherit.aes = TRUE, # don't extend up
                     color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed")
      
    }
    
    if(!is.null(fit$betaHat.se)){
      if(max(beta.hat.plt$upper.joint) > 0 & min(beta.hat.plt$lower.joint) < 0){
        plot_list[[r]] <- plot_list[[r]] +
          geom_segment(aes_string(x=xlim[1] - x_range, xend=xlim[2] + x_range,
                                  y=0,yend=0), inherit.aes = TRUE,
                       color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed")
      }
      colnames(beta.hat.plt) <- c("s", "beta.hat", "CI.lower.pointwise", "CI.upper.pointwise", "CI.lower.joint", "CI.upper.joint")
    }else{
      colnames(beta.hat.plt) <- c("s", "beta.hat")
    }
    
    res_list[[r]] <- beta.hat.plt # save data frame associated with rth fixed effect
    
  }
  
  plot_return <- do.call("grid.arrange", c(plot_list, nrow = num_row)) # plot
  plot_return # show plots
  
  if(return == TRUE){
    res_list$plot <- plot_return # save to returned object
    return(res_list)
  }else{
    return(plot_return)
  }
  
}
