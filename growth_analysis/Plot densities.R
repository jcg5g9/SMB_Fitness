
#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated posterior densities of parameters from  the VBGF  model.
#'
#' @param model A fit model
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param probs Lower and upper quantiles to use for plot limits if lower and upper are not specified.
#'
#' @return Returns and saves a figure with the posterior densities of parameters.
#' @export
plot_density <- function(model, file_name = NULL, probs = c(0.0125, 0.9875) ){
  
  library(latex2exp)
  
  
  # Line specifications
  posteriors_lwd <- rep(3, 2)
  posteriors_lty <- rep(1, 2)
  posteriors_col <- c("grey", 1)
  
  
  # Vars of interest
  vars <- c(
    "mu_linf", "mu_k", "mu_t0", 
    "linf_lineage[1]", "linf_lineage[2]", "post_linf_diff", 
    "k_lineage[1]", "k_lineage[2]", "post_k_diff", 
    "t0_lineage[1]", "t0_lineage[2]", "post_t0_diff", 
    "beta_linf[1]", "beta_linf[2]", 
    "beta_k[1]", "beta_k[2]", 
    "beta_t0[1]", "beta_t0[2]"
  )
  
  prior_vars <- c(
    "prior_mu_linf", "prior_mu_k", "prior_mu_t0", 
    "prior_linf_lineage[1]", "prior_linf_lineage[2]", "prior_linf_diff", 
    "prior_k_lineage[1]", "prior_k_lineage[2]", "prior_k_diff", 
    "prior_t0_lineage[1]", "prior_t0_lineage[2]", "prior_t0_diff", 
    "prior_beta_linf[1]", "prior_beta_linf[2]", 
    "prior_beta_k[1]", "prior_beta_k[2]", 
    "prior_beta_t0[1]", "prior_beta_t0[2]"
  )
  
  vars_latex <- list(TeX("$\\bar{L}_{\\infty}$"), TeX("$\\bar{K}$"), TeX("\\bar{a0}"), 
                     TeX("$L_{\\infty} \\, SMB$"), TeX("$L_{\\infty} \\, NB$"), TeX("\\Delta $L_{\\infty}$"),
                     TeX("$K \\, SMB$"), TeX("$K \\, NB$"), TeX("\\Delta K"),
                     TeX("$a0 \\, SMB$"), TeX("$a0 \\, NB$"), TeX("\\Delta a0"),
                     
                     TeX("$\\Beta_{L_{\\infty}} \\, Male$"), TeX("$\\Beta_{L_{\\infty}} \\, ER$"),
                     TeX("$\\Beta_{K} \\, Male$"), TeX("$\\Beta_{K} \\, ER$"),
                     TeX("$\\Beta_{a0} \\, Male$"), TeX("$\\Beta_{a0} \\, ER$")
  )
  
  
  ## Extract posterior densities
  draws <- as.data.frame(model)
  posterior_dens <- list()
  posterior_dens[[1]] <- apply(draws[,vars], 2, function(x) density(as.numeric(as.character(x)))) 
  posterior_dens[[2]] <- apply(draws[,prior_vars], 2, function(x) density(as.numeric(as.character(x)))) 
  
  
  ## X-limits
  xlow <- sapply(draws[,vars], function(x) quantile(x, probs= probs[1]))
  xup <- sapply(draws[,vars], function(x) quantile(x, probs= probs[2]))
  
  # - Linf range
  linf_ind <- c(1,4,5)
  xlow[linf_ind] <- min(xlow[linf_ind])
  xup[linf_ind] <- max(xup[linf_ind])
  
  # - K range
  k_ind <- c(2,7,8)
  xlow[k_ind] <- min(xlow[k_ind])
  xup[k_ind] <- max(xup[k_ind])
  
  # - A0 range
  a0_ind <- c(3,10,11)
  xlow[a0_ind] <- min(xlow[a0_ind])
  xup[a0_ind] <- max(xup[a0_ind])
  
  
  
  # Plot
  for(j in 1:(1 + as.numeric(!is.null(file_name)))){
    
    # PNG
    if(j == 2){
      filename <- paste0(file_name, "_posterior_density", ".png")
      png( file = filename , width=10, height = 110 / 25.4, family = "serif", units = "in", res = 300)
    }
    
    par(mfrow = c(2,ceiling(length(vars)/2) + 2))
    par( mar=c(3, 0.05 , 0.5 , 0.55) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    
    plot.new()
    
    
    # Loop through vars
    for(i in 1:length(vars)){
      
      
      # Plot them
      plot(NA,
           xlim = c(xlow[i], xup[i]),
           ylim = c(0, range(posterior_dens[[1]][[i]]$y)[2]),
           ylab = NA, xlab = vars_latex[[i]], yaxt = "n")
      
      
      lines(posterior_dens[[2]][[i]], lwd = posteriors_lwd[1], lty = posteriors_lty[1], col = posteriors_col[1]) # Prior
      lines(posterior_dens[[1]][[i]], lwd = posteriors_lwd[2], lty = posteriors_lty[2], col = posteriors_col[2]) # Posterior
      
      if(i == ceiling(length(vars)/2))  {
        plot.new()
        plot.new()
      }
      
      
      if(i %in% c(1, ceiling(length(vars)/2) + 1) ) {
        mtext(side = 2, "Density", line = 1, cex= 0.75)
      }
      
      
    }
    
    if(j > 1){ dev.off()}
  }
}

