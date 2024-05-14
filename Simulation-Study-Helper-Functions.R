# Generate fixed and random effects ----------------------------------------

set_fixed_effects <- function() {
  return( list( Beta = c(seq(from = 0.48, to = -0.48, length = dim(X)[2])), 
                S_coeff = c(1.0, 0.5, 0, 0), 
                gamma_const = log(1.0325) ) )
}


set_random_effects <- function() {
  
  # Ut - the temporal random effect 
  Ut <- 2.5 / (1 + exp( 0.0153 - 0.003455*num_obs ) )
  Ut <- Ut - Ut[which(day_constraint==1)]
  # plot(exp(Ut), type = "l")
  
  # Vt - the i.i.d. error term 
  Vt <- rnorm(n = length(num_obs), mean = 0, sd = 0.30) 
  
  # Mt - the temperature random effect  
  quad_effect <- function(x) { 5 * (x - 18)^2 * ( 4 * 10^4 + ifelse(x>18, 1, 0) * 5.33 * 10^3 + 
                                                    ifelse(x<18, 1, 0) * 4 * 10^4 ) / 10^9 }
  ## Ascending list of temperature values 
  Mt_val <- quad_effect(temp_bins + min( round( sim_data$kTemp.24hm.lag0 ) ) - 1)
  
  ## Duplicate simulated temperature values across time 
  Mt_bin <- seq( from = 1, to = length(Mt_val), by = 1)
  Mt <- vector( length = length(num_obs) )
  Mt_idx <- vector( length = length(num_obs) )
  for ( ii in 1:length(Mt) ) {
    for ( jj in 1:length(Mt_val)) { # if we are in the same bin, assign the same value  
      if ( isTRUE(all.equal( temp_values[ii] - min(temp_bins) + 1, Mt_bin[jj] )) ) { 
        Mt_idx[ii] <- jj # track bin
        Mt[ii] <- Mt_val[jj] # assign value
      }
    }
  }
  # plot(exp(Mt_val), type = "l")
  
  # Time-varying gamma: monotonically increasing (a logistic curve)
  logistic_func <- function(x, denom) ( 1 / (1 + exp(-x)) ) / denom
  gamma_logistic <- 1 * logistic_func( 3 - 0.0015 * (5211-num_obs), 10 )
  gamma_logistic <- gamma_logistic - gamma_logistic[which(day_constraint==1)] 
  # plot(exp(gamma_logistic)+log(1.0325), type = "l")
  
  # Time-varying gamma: quadratic function
  gamma_quad <- (num_obs - length(num_obs)/2)^2 / ( 1.5 * 250000000 )
  gamma_quad <- gamma_quad - gamma_quad[which(day_constraint==1)]
  # plot(exp(gamma_quad+log(1.0325)), type = "l")
  
  # Time-varying gamma: translated and scaled skew-normal density
  gamma_skew <- 100 * dsn( seq(from = 1, to = length(num_obs), by = 1), 
                           xi = length(num_obs) / 8, 
                           omega = length(num_obs) / 4, 
                           alpha = 2, log = FALSE)
  gamma_skew <- gamma_skew - gamma_skew[which(day_constraint==1)]
  # plot(exp(gamma_skew+log(1.0325)), type = "l")
  
  # Time-varying gamma: translated and scaled cubic function
  gamma_cubic <- (num_obs+1100-3000) * (num_obs+1100-3000) * (num_obs-2700-3000) / 4e11 - 
    ( (num_obs+1100-3000) * (num_obs+1100-3000) * (num_obs-2700-3000) / 4e11 )[which(day_constraint==1)]
  # plot(exp(gamma_cubic+log(1.0325)), type ="l")
  
  # Time varying gamma: step function 
  gamma_step <- vector( length = length(num_obs) )
  interval_length = length(num_obs) / 2
  gamma_step[ 1:(interval_length+180) ] <- 0
  gamma_step[ (interval_length+180+1):(2*interval_length) ]  <- 125 / interval_length 
  gamma_step <- 0.75 * ( gamma_step - gamma_step[which(day_constraint==1)] ) 
  # plot(exp(gamma_step+log(1.0325)), type = "l")
  
  return( list(Ut = Ut, Vt = Vt, 
               Mt = Mt, Mt_val = Mt_val, 
               gamma_logistic = gamma_logistic, 
               gamma_quad = gamma_quad,
               gamma_skew = gamma_skew, 
               gamma_cubic = gamma_cubic, 
               gamma_step = gamma_step ) )
}




# Simulate mortality counts ------------------------------------------------
mortality_const <- function() {
  lambda_const <- exp( X %*% fixed_effects$Beta + 
                         S %*% fixed_effects$S_coeff + 
                         pollutant * fixed_effects$gamma_const + 
                         random_effects$Vt + random_effects$Ut + random_effects$Mt ) 
  return( rpois( n = length(num_obs), 
                 lambda = lambda_const ) )  
}

# Time-varying gamma - logistic curve 
mortality_logistic <- function() {
  lambda_logistic <- exp( X %*% fixed_effects$Beta + 
                            S %*% fixed_effects$S_coeff + 
                            pollutant * fixed_effects$gamma_const + 
                            pollutant * random_effects$gamma_logistic  +
                            random_effects$Vt + random_effects$Ut + random_effects$Mt )  
  return( rpois( n = length(num_obs), 
                 lambda = lambda_logistic ) )
}


# Time-varying gamma - quadratic curve 
mortality_quad <- function() {
  lambda_quad <- exp( X %*% fixed_effects$Beta + 
                        S %*% fixed_effects$S_coeff + 
                        pollutant * fixed_effects$gamma_const +
                        pollutant * random_effects$gamma_quad +
                        random_effects$Vt + random_effects$Ut + random_effects$Mt )  
  return( rpois( n = length(num_obs), 
                 lambda = lambda_quad ) )
}

mortality_skew <- function() {
  lambda_skew <- exp( X %*% fixed_effects$Beta + 
                        S %*% fixed_effects$S_coeff + 
                        pollutant * fixed_effects$gamma_const + 
                        pollutant * random_effects$gamma_skew +
                        random_effects$Vt + random_effects$Ut + random_effects$Mt )  
  return( rpois( n = length(num_obs), 
                 lambda = lambda_skew ) )
}

mortality_cubic <- function() {
  lambda_cubic <- exp( X %*% fixed_effects$Beta + 
                         S %*% fixed_effects$S_coeff + 
                         pollutant * fixed_effects$gamma_const +
                         pollutant * random_effects$gamma_cubic +
                         random_effects$Vt + random_effects$Ut + random_effects$Mt )  
  return( rpois( n = length(num_obs), 
                 lambda = lambda_cubic ) )
}

mortality_step <- function() {
  lambda_step <- exp( X %*% fixed_effects$Beta + 
                        S %*% fixed_effects$S_coeff + 
                        pollutant * fixed_effects$gamma_const + 
                        pollutant * random_effects$gamma_step +
                        random_effects$Vt + random_effects$Ut + random_effects$Mt )  
  return(  rpois( n = length(num_obs), 
                  lambda = lambda_step ) )
}




# Summary table helper functions ------------------------------------------
fixed_effects_summary <- function(INLA_model) {
  return( t( sapply( INLA_model$marginals.fixed, function(xx) { 
    brinla::bri.density.summary(inla.tmarginal(xx, fun = exp)) } ) ) )
}

random_effects_summary <- function(INLA_model) {
  return( brinla::bri.hyperpar.summary(INLA_model) )
}

sim_table_fixed_effect <- function(INLA_model) {
  fe_summary <- fixed_effects_summary(INLA_model)
  return( c( fe_summary["Pollutant","mean"],
             fe_summary["Pollutant","q0.025"],
             fe_summary["Pollutant","q0.975"] ) )
}

sim_table_random_effect <- function(INLA_model) {
  re_summary <- random_effects_summary(INLA_model)
  return( c( re_summary["SD for Pollutant_RW2","mean"],
             re_summary["SD for Pollutant_RW2","q0.025"],
             re_summary["SD for Pollutant_RW2","q0.975"] ) )
}

sim_table_mk <- function(mk_sample) { 
  c( mean(mk_sample), 
     quantile( mk_sample, probs = 0.025 ), 
     quantile( mk_sample, probs = 0.975) ) 
}


table_entry <- function(mean, lower, upper, ii) { 
  mean <- format( round( mean[ii], 3), nsmall = 3)
  lower <- format( round( lower[ii], 3), nsmall = 3)
  upper <- format( round( upper[ii], 3), nsmall = 3)
  
  if ( mean < 0.001 ) mean <- "<0.001"
  if ( lower < 0.001 ) lower <- "<0.001"
  if ( upper < 0.001 ) upper <- "<0.001"
  
  return( paste0( mean, " (", lower, ", ", upper, ")") )
}

table_entry_momentum <- function(summary_vector) { 
  mean <- round( summary_vector[1], 4)
  lower <- round( summary_vector[2], 4)
  upper <- round( summary_vector[3], 4)
  
  if ( mean < 0.0001 ) mean <- "<0.0001"
  if ( lower < 0.0001 ) lower <- "<0.0001"
  if ( upper < 0.0001 ) upper <- "<0.0001"
  
  return( paste0( mean, " (", lower, ", ", upper, ")") )
}

table_entry_stability <- function(summary_vector) { 
  mean <- round( summary_vector[1], 4)
  lower <- round( summary_vector[2], 4)
  upper <- round( summary_vector[3], 4)
  
  if ( mean < 0.0001 ) mean <- "<0.0001"
  if ( lower < 0.0001 ) lower <- "<0.0001"
  if ( upper < 0.0001 ) upper <- "<0.0001"
  
  return( paste0( mean, " (", lower, ", ", upper, ")") )
}


# Fixed Effects tables ---------------------------------------------------
summary_FE <- function(fe, re, stat) {
  Table <- c( fe["Pollutant", stat],
              fe['(Intercept)', stat], 
              fe['Tuesday', stat], 
              fe['Wednesday', stat], 
              fe['Thursday', stat], 
              fe['Friday', stat], 
              fe['Saturday', stat], 
              fe['Sunday', stat], 
              fe['cos12', stat], 
              fe['cos6', stat], 
              fe['sin12', stat], 
              fe['sin6', stat], 
              10 * re['SD for Temperature_RW2', stat],
              365.25 * re['SD for daily_time', stat],
              re['SD for daily_iid', stat] )
}

column_FE <- function(mean, lower, upper) {
  t( tibble(
    "Ozone: $\\beta_{O3}$" = table_entry( mean, lower, upper, 1 ),
    "Intercept: $\\beta_{\\text{0}}$" = table_entry( mean, lower, upper, 2 ),
    "Tuesday: $\\beta_{\\text{Tues}}$" = table_entry( mean, lower, upper, 3 ),
    "Wednesday: $\\beta_{\\text{Wed}}$" = table_entry( mean, lower, upper, 4 ),
    "Thursday: $\\beta_{\\text{Thurs}}$" = table_entry( mean, lower, upper, 5 ),
    "Friday: $\\beta_{\\text{Fri}}$" = table_entry( mean, lower, upper, 6 ),
    "Saturday: $\\beta_{\\text{Sat}}$" = table_entry( mean, lower, upper, 7 ),
    "Sunday: $\\beta_{\\text{Sun}}$" = table_entry( mean, lower, upper, 8 ),
    "$\\beta_{\\text{cos12}}$" = table_entry( mean, lower, upper, 9 ),
    "$\\beta_{\\text{cos6}}$" = table_entry( mean, lower, upper, 10 ),
    "$\\beta_{\\text{sin12}}$" = table_entry( mean, lower, upper, 11 ),
    "$\\beta_{\\text{sin6}}$" = table_entry( mean, lower, upper, 12 ),
    "Temperature: 10$\\nu_{1}$" = table_entry( mean, lower, upper, 13 ),
    "Time: 365.25$\\nu_{2}$" = table_entry( mean, lower, upper, 14 ),
    "IID Error Term: $\\tau$" = table_entry( mean, lower, upper, 15 ) ) )
}

Appendix_Table_SimStudy_FE <- function( Table, caption, label ) {
  knitr::kable( Table, caption = caption, escape = FALSE, longtable = TRUE,
                row.names = TRUE, booktabs = TRUE, format = "latex", 
                label = label ) %>%
    kableExtra::kable_styling("striped", full_width = FALSE, latex_options = "hold_position") %>%
    pack_rows("Pollutant Effects", start_row = 1, end_row = 1) %>%
    pack_rows("Day-of-the-Week Effects", start_row = 2, end_row = 8) %>%
    pack_rows("Seasonal Effects", start_row = 9, end_row = 12) %>%
    pack_rows("Standard Deviation", start_row = 13, end_row = 15) %>%
    add_header_above( c(" " = 1, "Constant Gamma" = 1, "Logistic Gamma" = 1, "Quadratic Gamma" = 1) )
}




# Random walk tables -------------------------------------------------------

summary_RW <- function(fe, re, stat) {
  Table <- c( fe['Pollutant', stat], 
              365.25 * re['SD for Pollutant_RW2', stat],
              fe['(Intercept)', stat], 
              fe['Tuesday', stat], 
              fe['Wednesday', stat], 
              fe['Thursday', stat], 
              fe['Friday', stat], 
              fe['Saturday', stat], 
              fe['Sunday', stat], 
              fe['cos12', stat], 
              fe['cos6', stat], 
              fe['sin12', stat], 
              fe['sin6', stat], 
              10 * re['SD for Temperature_RW2', stat],
              365.25 * re['SD for daily_time', stat],
              re['SD for daily_iid', stat] )
}

column_RW <- function(mean, lower, upper) {
  t( tibble( 
    "Ozone: $\\beta_{O3}$" = table_entry( mean, lower, upper, 1 ),
    "Ozone: 365.25$\\sigma$" = table_entry( mean, lower, upper, 2 ),
    "Intercept: $\\beta_{0}$" = table_entry( mean, lower, upper, 3 ),
    "Tuesday: $\\beta_{\\text{Tues}}$" = table_entry( mean, lower, upper, 4 ),
    "Wednesday: $\\beta_{\\text{Wed}}$" = table_entry( mean, lower, upper, 5 ),
    "Thursday: $\\beta_{\\text{Thurs}}$" = table_entry( mean, lower, upper, 6 ),
    "Friday: $\\beta_{\\text{Fri}}$" = table_entry( mean, lower, upper, 7 ),
    "Saturday: $\\beta_{\\text{Sat}}$" = table_entry( mean, lower, upper, 8 ),
    "Sunday: $\\beta_{\\text{Sun}}$" = table_entry( mean, lower, upper, 9 ),
    "$\\beta_{cos12}$" = table_entry( mean, lower, upper, 10 ),
    "$\\beta_{cos6}$" = table_entry( mean, lower, upper, 11 ),
    "$\\beta_{sin12}$" = table_entry( mean, lower, upper, 12 ),
    "$\\beta_{sin6}$" = table_entry( mean, lower, upper, 13 ),
    "Temperature: 10$\\nu_{1}$" = table_entry( mean, lower, upper, 14 ),
    "Time: 365.25$\\nu_{2}$" = table_entry( mean, lower, upper, 15 ),
    "IID Error Term: $\\tau$" = table_entry( mean, lower, upper, 16 ) ) )
}

Appendix_Table_SimStudy_RW <- function( Table, caption, label ) {
  knitr::kable( Table, caption = caption, escape = FALSE, longtable = TRUE,
                row.names = TRUE, booktabs = TRUE, format = "latex", 
                label = label ) %>%
    kableExtra::kable_styling("striped", full_width = FALSE, latex_options = "hold_position") %>%
    pack_rows("Pollution Effects", start_row = 1, end_row = 2) %>%
    pack_rows("Day-of-the-Week Effects", start_row = 3, end_row = 9) %>%
    pack_rows("Seasonal Effects", start_row = 10, end_row = 13) %>%
    pack_rows("Standard Deviation", start_row = 14, end_row = 16) %>%
    add_header_above( c(" " = 1, "Constant Gamma" = 1, "Logistic Gamma" = 1, "Quadratic Gamma" = 1) )
}




# Helper Functions: Posterior densities --------------------------------------------------


s.weights <- function(inla_model) {
  weights <- c()
  for (ii in 1:inla_model$misc$configs$nconfig) {
    weights[ii] <- inla_model$misc$configs$config[[ii]]$log.posterior
  }
  return( exp(weights) / sum(exp(weights)) )
}

sampSizes <- function(inla_model, n = 1000) {
  return(drop(rmultinom(1, n, s.weights(inla_model))))
}

extractAllMeans <- function(inla_model) {
  return(purrr::map(inla_model$misc$configs$config, function(xx) xx$mean))
}

extractAllPrecMat <- function(inla_model) {
  return(purrr::map(inla_model$misc$configs$config, function(xx) xx$Q))
}

sample_constr_rw <- function(INLA_model) {
  # Number of: configurations, constraints, sample size per configuration
  N <- INLA_model$misc$configs$nconfig # N = number of configurations
  N_constraint <- INLA_model$misc$configs$constr$nc # N_constrained = number of constraints
  n <- sampSizes(INLA_model); n_samp <- n[which(n>0)]
  
  # Remove configurations with 0 samples 
  config_indices <- which(n > 0)
  n <- n[config_indices] # only sample configurations with >0 samples 
  N <- length(n) # N = number of configurations with >0 samples being drawn 
  
  # Set linear constraints (Ax = e)
  model_A <- INLA_model$misc$configs$constr$A # k x n matrix of linear constraints
  model_e <- INLA_model$misc$configs$constr$e # size k vector setting linear constraints
  model_constr = list(A = model_A, e = model_e)
  
  # Extract posterior means and (full) precision matrices from INLA model
  model_mu <- extractAllMeans(INLA_model)[config_indices] # posterior means 
  model_Q <- extractAllPrecMat(INLA_model)[config_indices] # posterior precision matrices 
  
  # Sample from the constrained random walk
  constrained_samples <- vector("list", length = N)
  for (ii in 1:N) {
    constrained_samples[[ii]] <- INLA:::inla.qsample( n = n[ii],  mu = model_mu[[ii]], 
                                                      Q = model_Q[[ii]], constr = model_constr, 
                                                      num.threads = 6 )
  }
  return(constrained_samples)
}


posterior_sample_function <- function(INLA_model){
  constrained_sample <- sample_constr_rw(INLA_model)
  posterior_sample <- do.call(cbind, constrained_sample)
  return(posterior_sample)
}

density_sigma <- function(INLA_model, post_sample = posterior_sample) { 
  
  gamma_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")])
  posterior_sample_gamma <- post_sample[gamma_indices-6210,] # edit
  
  posterior_sample_sd <- apply( posterior_sample_gamma, 2, sd )  
  return( posterior_sample_sd )  
}


Stability_Summary <- function(INLA_model, post_sample = posterior_sample) {  
  c( mean( density_sigma( INLA_model, post_sample ) ),
     quantile( density_sigma( INLA_model, post_sample ), c(0.025, 0.975) ) )
}

sim_table_mk <- function(mk_sample) { 
  c( mean(mk_sample), 
     quantile( mk_sample, probs = 0.025 ), 
     quantile(mk_sample, probs = 0.975) ) 
}

kernel_dens <- function( sample, set_bandwidth ) { 
  return( bkde( sample, range.x = c(0,1.005), 
                kernel = "normal", gridsize = 5 * 1000,
                bandwidth = set_bandwidth ) )
}

kernel_dens_mom <- function( sample, set_bandwidth ) { 
  return( bkde( sample, range.x = c(0,1.005), 
                kernel = "normal", gridsize = 1 * 1000,
                bandwidth = set_bandwidth ) )
}


# Posterior Plots: Summary Statistics ---------------------------------------


sample_posterior_momentum <- function(INLA_model, post_sample = posterior_sample) { 
  
  gamma_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")])
  posterior_sample_gamma <- post_sample[gamma_indices-6210,]
  
  Nsim <- dim(posterior_sample_gamma)[2]
  N <- dim(posterior_sample_gamma)[1]
  Momentum_Posterior <- Momentum( Nsim = Nsim, N = N, theSim = posterior_sample_gamma )
  
  return(Momentum_Posterior)
}

plot_posterior_sigma <- function( model_1, model_2, model_3, 
                                  legend_position, legend_size, 
                                  first_simStudy, legend_include = FALSE,
                                  x_min = 0, x_max = 0.02, offset = 0 ) {
  
  colors <- brewer.pal(n = 8, name = "Dark2")[5:8]
  
  model_1_sigma <- kernel_dens( density_sigma( INLA_model = model_1, 
                                               post_sample = posterior_sample_function( model_1 ) ),
                                set_bandwidth = 0.0002 ) 
  model_2_sigma <- kernel_dens( density_sigma( INLA_model = model_2, 
                                               post_sample = posterior_sample_function( model_2 ) ),
                                set_bandwidth = 0.0002 )
  model_3_sigma <- kernel_dens( density_sigma( INLA_model = model_3, 
                                               post_sample = posterior_sample_function( model_3 ) ),
                                set_bandwidth = 0.0002 )
  
  x_min <- x_min; x_max = x_max
  y_max <- max( model_1$y, model_2_sigma$y, model_3_sigma$y )
  
  ### Plot 3 posterior densities
  par( mar = c(3,2,2,0.1), xaxs = "r", yaxs = "r" )
  plot( x = model_1_sigma$x, y = model_1_sigma$y,
        xlim = c(x_min, x_max), ylim = c(0, y_max + offset), xlab = "", ylab = "Density", 
        type = 'l', col = colors[1], lty = 'solid', lwd = 3 )
  lines( x = model_2_sigma$x, y = model_2_sigma$y,
         type = 'l', col = colors[2], lty = 'solid', lwd = 3 )
  lines( x = model_3_sigma$x, y = model_3_sigma$y,
         type = 'l', col = colors[3], lty = 'solid', lwd = 3 )
  
  if ( legend_include ) {  
    legend( legend_position, 
            legend = c( ifelse( first_simStudy, "Constant Effect", "Skewed Effect"),
                        ifelse( first_simStudy, "Logistic Effect", "Cubic Effect"),
                        ifelse( first_simStudy, "Quadratic Effect", "Step Function Effect") ),
            col = c(colors[1], colors[2], colors[3]), 
            lty = 1, lwd = 3, cex = legend_size, box.lty = 0)
  }
}


plot_posterior_momentum <- function( momentum_sample_1, momentum_sample_2, 
                                     momentum_sample_3, legend_position, 
                                     legend_size, legend_include = FALSE, 
                                     offset, first_simStudy ) {
  
  colors <- brewer.pal(n = 8, name = "Dark2")[5:8]
  
  density_1 <- kernel_dens_mom(momentum_sample_1, 0.05)
  density_2 <- kernel_dens_mom(momentum_sample_2, 0.0015)
  density_3 <- kernel_dens_mom(momentum_sample_3, 0.05)
  
  x_min <- min( density_1$x, density_2$x, density_3$x )
  x_max <- max( density_1$x, density_2$x, density_3$x )
  y_max <- max( density_1$y, density_2$y, density_3$y )
  
  ### Plot posteriors
  par( mar = c(3,2,2,0.1), xaxs = "r", yaxs = "r" )
  plot( x = density_1$x, y = density_1$y, 
        xlab = "", ylab = "", main = "", 
        lwd = 3, col = colors[1], type = "l",
        ylim = c(0, y_max+offset), xlim = c(0,1) )
  lines( x = density_2$x, y = density_2$y, 
         lwd = 3, col = colors[2], main = "" )
  lines( x = density_3$x, y = density_3$y, 
         lwd = 3, col = colors[3], main = "" )
  
  ### Add legend 
  if ( legend_include ) {
    legend( legend_position, 
            legend = c( ifelse( first_simStudy, "Constant Effect", "Skewed Effect"),
                        ifelse( first_simStudy, "Logistic Effect", "Cubic Effect"),
                        ifelse( first_simStudy, "Quadratic Effect", "Step Function Effect") ),
            col = c(colors[1], colors[2], colors[3]), 
            lty = 1, lwd = 3, cex = legend_size, box.lty = 0)
  }
}





# Plot RW2 Gamma ------------------------------------------------
calculateMeanRandomEffect <- function(INLA_model) {
  
  mean_re <- list( Gamma_t = rep(1, length = length(num_obs)) )
  for ( kk in 1:length(num_obs) ) {
    if ( kk != which(day_constraint == 1) ) {
      mean_re$Gamma_t[kk] <- inla.emarginal( INLA_model$marginals.random$Pollutant_RW2[[kk]], fun = exp ) 
    }
  }
  
  return(mean_re)
}


RE_indices_for_plotting <- function(INLA_model) {
  time_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="daily_time")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="daily_time")])
  temp_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Temperature_RW2")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="Temperature_RW2")])
  pollutant_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")])
  
  time_indices <- time_indices; temp_indices <- temp_indices; pollutant_indices <- pollutant_indices
  
  ### Error checking - is the constrained index an index of that random effect?
  model_A <- INLA_model$misc$configs$constr$A 
  constraint_idx_time <- which(model_A[1,] == 1)
  constraint_idx_temp <- which(model_A[2,] == 1)
  constraint_idx_pollutant <- which(model_A[3,] == 1)
  # Give an error if it is not.
  # if ( !(constraint_idx_time %in% time_indices) ) { stop("Check constraint_idx of 'daily_time'.") } 
  # if ( !(constraint_idx_temp %in% temp_indices) ) { stop("Check constraint_idx of 'Temperature_RW2'.") } 
  # if ( !(constraint_idx_pollutant %in% pollutant_indices) ) { stop("Check constraint_idx of 'Pollutant_RW2'.") } 
  
  return( list("Time" = time_indices,
               "Temp" = temp_indices,
               "Pollutant" = pollutant_indices) )
}


Calc_Global_CI_lincomb <- function(constr_samples, indices, fe_sample) {
  ### Concatenate samples for 'indices' into one matrix + 
  ### add sample from fixed effect to each column of said matrix 
  sim_matrix_pollutant <- matrix(NA, nrow = nrow(constr_samples[indices,]),
                                 ncol = ncol(constr_samples[indices,]) )
  
  for (ii in 1:ncol(sim_matrix_pollutant)) {
    sim_matrix_pollutant[,ii] <- constr_samples[indices,ii] + fe_sample[ii]
  }
  
  cset <- create_curve_set(list( r = 1:dim(sim_matrix_pollutant)[1], 
                                 obs = rowMeans(sim_matrix_pollutant),
                                 sim_m = sim_matrix_pollutant ) )
  cr <- central_region(cset, coverage = 0.80) 
  return( list("Global_LCL" = cr$lo,
               "Global_UCL" = cr$hi) )
}


plot_Gamma_lincomb_RW <- function( INLA_model, Global_CI_spacing = 10, mean_spacing = 10,
                                   first_sample_index = 0, cex_size = 1, threads = 6, 
                                   sim_const = FALSE, sim_logistic = FALSE, sim_quad = FALSE, 
                                   sim_skew = FALSE, sim_cubic = FALSE, sim_step = FALSE, 
                                   full_legend = FALSE,  legend_position = "top", 
                                   legend_mean_only = FALSE, main_paper = FALSE, 
                                   mean_color = "grey25", post_sample = posterior_sample ) { 
  
  ### Read in INLA object and calculate mean random effect 
  Gammat <- sim_data$Date
  meanRandomEffect_RW <- calculateMeanRandomEffect(INLA_model)
  
  ### Sample constrained random walk for pollutant (from INLA model)
  gamma_indices <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")]-1 +
    (1:INLA_model$misc$configs$contents$length[which(INLA_model$misc$configs$contents$tag=="Pollutant_RW2")])
  posterior_sample_gamma <- post_sample[gamma_indices,]
  
  ### Posterior samples of random walk sorted in increasing order of Momentum statistic
  Nsim <- dim(posterior_sample_gamma)[2]
  N <- dim(posterior_sample_gamma)[1]
  real_mk_sample <- Momentum( Nsim = Nsim, N = N, theSim = posterior_sample_gamma )
  ps_gamma_mk <- rbind(posterior_sample_gamma, real_mk_sample)
  ps_gamma_mk_sorted <- ps_gamma_mk[, order(ps_gamma_mk[6211, ])]
  
  ### Posterior samples of random walk sorted in increasing order of squared second differences
  ps_gamma_sdsq <- colMeans(apply(posterior_sample_gamma, 2, function(x) { diff(diff(x)) } ))
  ps_gamma_sdsq <- rbind(posterior_sample_gamma, ps_gamma_sdsq)
  ps_gamma_sdsq_sorted <- ps_gamma_sdsq[, order(ps_gamma_sdsq[6211, ])]
  
  gamma_idx <- INLA_model$misc$configs$contents$start[which(INLA_model$misc$configs$contents$tag=="Pollutant")]
  gamma_fixed_effect <- post_sample[gamma_idx,] 
  
  Global_CI <- Calc_Global_CI_lincomb(constr_samples = posterior_sample,
                                      indices = RE_indices_for_plotting(INLA_model)$Pollutant,
                                      fe_sample = gamma_fixed_effect)
  indices_of_quantiles <- seq( from = first_sample_index, to = dim(posterior_sample_gamma)[2],
                               by = dim(posterior_sample_gamma)[2] / 5 )[2:5]
  
  if ( main_paper ) { range_min = 0.98; range_max <- 1.02 }
  
  if ( sim_const ) { range_min <- 1.02; range_max <- 1.05 }
  if ( sim_logistic ) { range_min <- 0.96; range_max <- 1.10 }
  if ( sim_quad ) { range_min <- 1.02; range_max <- 1.07 }
  
  if ( sim_skew ) { range_min <- 1.00; range_max <- 1.09 }
  if ( sim_cubic ) { range_min <- 0.97; range_max <- 1.08 }
  if ( sim_step ) { range_min <- 1.01; range_max <- 1.11 }
  
  
  ### Plot posterior samples of Gamma(t), based on the quartiles/quintiles of:
  ### their posterior MK statistic + second differences 
  post_samples_to_plot <- cbind( exp( ps_gamma_sdsq_sorted[1:dim(posterior_sample_gamma)[1],indices_of_quantiles] ) * 
                                   inla.zmarginal(inla.tmarginal(INLA_model$marginals.fixed[[12]], fun = exp), silent = TRUE)[[1]], 
                                 exp( ps_gamma_mk_sorted[1:dim(posterior_sample_gamma)[1],indices_of_quantiles] ) * 
                                   inla.zmarginal(inla.tmarginal(INLA_model$marginals.fixed[[12]], fun = exp), silent = TRUE)[[1]] )
  
  par( mar = c(3,2,3,0.1), xaxs = "r", yaxs = "r" )
  matplot( x = Gammat, y = post_samples_to_plot, col = "grey80", xaxt = "n",
           type = 'l', lwd = 1,  lty = 1, ylab = "Relative Risk", 
           ylog = TRUE, cex.main = cex_size, cex.axis = cex_size, 
           cex.lab = cex_size, ylim = c( range_min, range_max ) )
  axis.Date(side = 1, Gammat, format = "%Y")
  
  ### Plot true effect (if simulated)
  if (sim_const) {
    sim_const_effect <- exp( fixed_effects$gamma_const + rep(0,6210) )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_const_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_logistic) {
    sim_logistic_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_logistic )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_logistic_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_quad) {
    sim_quad_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_quad )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_quad_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_skew) {
    sim_skew_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_skew )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_skew_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_cubic) {
    sim_cubic_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_cubic )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_cubic_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_step) {
    sim_step_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_step )
    lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
           y = sim_step_effect[seq(0, length(Gammat), by=Global_CI_spacing)], 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  ### Plot gamma(t) smaple mean 
  lines( x = Gammat[seq(0,length(INLA_model$marginals.lincomb.derived), by = mean_spacing)], 
         y = unlist(lapply(INLA_model$marginals.lincomb.derived[seq(0,length(INLA_model$marginals.lincomb.derived), 
                                                                    by = mean_spacing)], 
                           inla.emarginal, fun = exp)),
         lwd = 3, col = mean_color )
  
  ### Plot gamma(t) global envelopes
  lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
         y = exp( Global_CI$Global_LCL )[seq(0, length(Global_CI$Global_LCL), by=Global_CI_spacing)], 
         type = "l", lwd = 3, lty = "dashed", col = blue_color )
  lines( x = Gammat[seq(0, length(Gammat), by=Global_CI_spacing)], 
         y = exp( Global_CI$Global_UCL )[seq(0, length(Global_CI$Global_UCL), by=Global_CI_spacing)], 
         type = "l", lwd = 3, lty = "dashed", col = blue_color )
  
  ### Add legend
  if ( full_legend ) {
    legend( legend_position, c("Mean", "True Value", "Global Envelope", "Samples"), 
            lty = c("solid", "solid", "dashed", "solid"), lwd = 3, cex = cex_size,
            col = c(mean_color, red_color, blue_color, "grey80" ), bty = "n", ncol = 2 ) 
  }
  
  if ( legend_mean_only ) {
    legend( legend_position, c("Mean"), lty = c("solid"), lwd = 3, 
            cex = cex_size, col = mean_color, bty = "n", ncol = 2 ) 
  }
}




# Plot Fixed Effects Gamma ------------------------------------------------
calculate_year_effect_mean <- function(INLA_annual_model, lincombs) {
  
  year_effect <- vector(length = length(year_count))
  year_estimates <- vector(length = length(unique(year_count))-1)
  year_estimates_lincomb <- vector(length = length(unique(year_count)))
  
  if (lincombs) {
    for ( ll in 1:length(year_estimates_lincomb) ) { 
      year_estimates_lincomb[ll] <- inla.emarginal( INLA_annual_model$marginals.lincomb.derived[[ll]], fun = exp )
    }
    for ( ll in 1:length(year_count) ) {
      year_effect[ll] <- year_estimates_lincomb[year_count[ll]]
    }
    return(year_effect)
  }
  
}

calculate_year_effect_quantile <- function(INLA_annual_model, quantile, lincombs) {
  
  year_effect <- vector(length = length(year_count))
  year_estimates <- vector(length = length(unique(year_count))-1)
  year_estimates_lincomb <- vector(length = length(unique(year_count)))
  
  if (lincombs) {
    for (ll in 1:length(year_estimates_lincomb)) { 
      year_estimates_lincomb[ll] <- exp( inla.qmarginal( quantile, INLA_annual_model$marginals.lincomb.derived[[ll]] ) )
    }
    for (ll in 1:length(year_count)) {
      year_effect[ll] <- year_estimates_lincomb[year_count[ll]]
    }
    return(year_effect)
  }
  
}



plot_Gamma_lincomb_annual <- function( INLA_fixedeffects_model,
                                       sim_const = FALSE, sim_logistic = FALSE, sim_quad = FALSE,
                                       sim_skew = FALSE, sim_cubic = FALSE,  sim_step = FALSE, 
                                       full_legend = FALSE, legend_mean_only = FALSE,
                                       mean_color = alpha("black", 0.6), lincomb_bool ) { 
  
  ### Calculate mean and quantiles 
  gamma_mean_effect <- calculate_year_effect_mean( INLA_annual_model = INLA_fixedeffects_model,
                                                   lincombs = lincomb_bool ) 
  gamma_lower_quantile <- calculate_year_effect_quantile( INLA_annual_model = INLA_fixedeffects_model,
                                                          quantile = 0.025, lincombs = lincomb_bool ) 
  gamma_upper_quantile <- calculate_year_effect_quantile( INLA_annual_model = INLA_fixedeffects_model,
                                                          quantile = 0.975, lincombs = lincomb_bool ) 
  
  
  if ( sim_const ) { range_min <- 1.02; range_max <- 1.05 }
  if ( sim_logistic ) { range_min <- 0.96; range_max <- 1.10 }
  if ( sim_quad ) { range_min <- 1.02; range_max <- 1.07 }
  
  if ( sim_skew ) { range_min <- 1.00; range_max <- 1.09 }
  if ( sim_cubic ) { range_min <- 0.97; range_max <- 1.08 }
  if ( sim_step ) { range_min <- 1.01; range_max <- 1.11 }
  
  ### Plot 
  par( mar = c(3,2,3,0.1), xaxs = "r", yaxs = "r" )
  plot( x = sim_data$Date, y = gamma_mean_effect, xlab = "", ylab = "",
        col = mean_color, type = "l", lwd = 2, ylim = c( range_min, range_max) ) 
  lines( x = sim_data$Date, y = gamma_lower_quantile, 
         type = "l", lwd = 1, lty = "dashed", col = blue_color ) 
  lines( x = sim_data$Date, y = gamma_upper_quantile, 
         type = "l", lwd = 1, lty = "dashed", col = blue_color ) 
  
  if (sim_const) {
    sim_const_effect <- exp( fixed_effects$gamma_const + rep(0,6210) )
    lines( x = sim_data$Date, y = sim_const_effect, 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_logistic) {
    sim_logistic_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_logistic )
    lines( x = sim_data$Date, y = sim_logistic_effect, 
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_quad) {
    sim_quad_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_quad )
    lines( x = sim_data$Date, y = sim_quad_effect,
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_skew) {
    sim_skew_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_skew )
    lines( x = sim_data$Date, y = sim_skew_effect,
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_cubic) {
    sim_cubic_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_cubic )
    lines( x = sim_data$Date, y = sim_cubic_effect,
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  if (sim_step) {
    sim_step_effect <- exp( fixed_effects$gamma_const + random_effects$gamma_step )
    lines( x = sim_data$Date, y = sim_step_effect,
           type = "l", lwd = 3, lty = "solid", col = red_color )
  }
  
  ### Add legend
  if ( full_legend ) {
    legend( "top", c("Mean", "True Value", "Credible Interval"),
            lty = c("solid", "solid", "dashed"), 
            lwd = 3, cex = 1, bty = "n", ncol = 2,
            col = c(alpha("black", 0.6), red_color, blue_color) )  
  }
  
  if ( legend_mean_only ) {
    legend( "top", c("Mean"), lty = c("solid"), 
            lwd = 3, cex = 1, bty = "n", ncol = 2,  
            col = alpha("black", 0.6) ) 
  }
}



# Fit INLA models --------------------------------------------------------
sim_data_INLA <- function() {
  return( data.frame( 
    # Response variable
    "mortality_const" = mortality_const(),
    "mortality_logistic" = mortality_logistic(),
    "mortality_quad" = mortality_quad(), 
    "mortality_skew" = mortality_skew(), 
    "mortality_cubic" = mortality_cubic(), 
    "mortality_step" = mortality_step(), 
    # Regression Matrix
    "Intercept" = X[,1],
    "Tuesday" = X[,2],
    "Wednesday" = X[,3],
    "Thursday" = X[,4], 
    "Friday" = X[,5], 
    "Saturday" = X[,6],
    "Sunday" = X[,7],
    "cos12" = S[,1],
    "cos6" = S[,2],
    "sin12" = S[,3],
    "sin6" = S[,4],
    "Pollutant" = pollutant, # main effect of O3
    # Random Effects
    "daily_time" = 1:length(num_obs),
    "daily_iid" = 1:length(num_obs),
    "Temperature_RW2" = temperature - min(temperature) + 1,
    "Pollutant_RW2" = 1:length(num_obs) ) )
  # Track year
  # "year_count" = year_count ) )
}


make_lincombs_fixedeffect <- function() { 
  lc1 = inla.make.lincomb( "Pollutant" = 1 ); names(lc1) = "lc1"
  lc2 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)2" = 1 ); names(lc2) = "lc2"
  lc3 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)3" = 1 ); names(lc3) = "lc3"
  lc4 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)4" = 1 ); names(lc4) = "lc4"
  lc5 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)5" = 1 ); names(lc5) = "lc5"
  lc6 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)6" = 1 ); names(lc6) = "lc6"
  lc7 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)7" = 1 ); names(lc7) = "lc7"
  lc8 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)8" = 1 ); names(lc8) = "lc8"
  lc9 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)9" = 1 ); names(lc9) = "lc9"
  lc10 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)10" = 1 ); names(lc10) = "lc10"
  lc11 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)11" = 1 ); names(lc11) = "lc11"
  lc12 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)12" = 1 ); names(lc12) = "lc12"
  lc13 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)13" = 1 ); names(lc13) = "lc13"
  lc14 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)14" = 1 ); names(lc14) = "lc14"
  lc15 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)15" = 1 ); names(lc15) = "lc15"
  lc16 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)16" = 1 ); names(lc16) = "lc16"
  lc17 = inla.make.lincomb( "Pollutant" = 1, "Pollutant:factor(year_count)17" = 1 ); names(lc17) = "lc17"
  
  return( c( lc1, lc2, lc3, lc4, lc5, lc6, lc7, lc8, lc9, 
             lc10, lc11, lc12, lc13, lc14, lc15, lc16, lc17 ) )
}



fit_dummy_model <- function(response, city_data) {
  dummy_model = inla( eval(substitute( response ~ Tuesday ) ),
                      data = city_data, family = 'poisson',
                      silent = 2L )
  
  m = get("inla.models", INLA:::inla.get.inlaEnv())
  m$latent$rw2$min.diff = NULL; m$latent$rw1$min.diff = NULL
  assign("inla.models", m, INLA:::inla.get.inlaEnv())
  
}



fit_RW2_model <- function( response, city_data, step_factor, step_len, init_theta = 4, 
                           day_constraint, temp_constraint, h_val = 0.005 ) {
  inla( eval(substitute( response ~ Tuesday + Wednesday + Thursday + Friday 
                         + Saturday + Sunday + cos12 + cos6 + sin12 + sin6 + Pollutant 
                         + f(daily_time, 
                             model = 'rw2', 
                             values = unique(daily_time), 
                             scale.model = TRUE,
                             constr = FALSE, 
                             extraconstr = list(A = matrix(day_constraint, nrow = 1), e = 0),
                             hyper = list(theta = list(prior = 'pc.prec', initial = 4,
                                                       param = c(log(2)/2, 0.50))))
                         + f(daily_iid, 
                             model = 'iid',
                             values = unique(daily_time), 
                             hyper = list(theta = list(prior = 'pc.prec', initial = 4, 
                                                       param = c(log(2)/2, 0.50))))
                         + f(Temperature_RW2, 
                             model = 'rw2', 
                             values = temp_bins,
                             scale.model = TRUE,
                             constr = FALSE, 
                             extraconstr = list(A = matrix(temp_constraint, nrow = 1), e = 0), 
                             hyper = list(theta = list(prior = 'pc.prec', initial = 4, 
                                                       param = c(log(2)/2, 0.50))))
                         + f(Pollutant_RW2, Pollutant, 
                             model = 'rw2', 
                             values = num_obs,
                             scale.model = TRUE,
                             constr = FALSE, 
                             extraconstr = list(A = matrix(day_constraint, nrow = 1), e = 0),
                             hyper = list(theta = list(prior = 'pc.prec', initial = init_theta,
                                                       param = c(log(2)/2, 0.50)))) ) ),
        data = city_data, 
        family = 'poisson', 
        # marginal posterior density of gamma_{0} + gamma_{t}
        lincomb = inla.make.lincombs( Pollutant = rep(1,length(num_obs)), 
                                      Pollutant_RW2 = diag(1,length(num_obs)) ),
        control.predictor = list( compute = TRUE ),
        control.fixed = list( prec = 1, 
                              prec.intercept = 1 ),
        control.inla = list( strategy = "simplified.laplace",
                             improved.simplified.laplace = TRUE,
                             restart = 5,
                             stencil = 9,
                             fast = TRUE,
                             step.len = step_len,
                             step.factor = step_factor,
                             h = h_val ),
        control.compute = list( config = TRUE, 
                                dic = TRUE, 
                                waic = TRUE,
                                openmp.strategy = "huge", 
                                smtp = "taucs" ),
        silent = 2L ) # close INLA call 
} # close function definition 




fit_FE_model <- function( response, city_data, 
                          day_constraint, temp_constraint,
                          step_factor, step_len ) {
  inla( eval( substitute( response ~ Tuesday + Wednesday + Thursday + Friday
                          + Saturday + Sunday + cos12 + cos6 + sin12 + sin6 + Pollutant 
                          + Pollutant:factor(year_count) + 
                            + f(daily_time, 
                                model = 'rw2', 
                                values = num_obs, 
                                scale.model = TRUE,
                                constr = FALSE, 
                                extraconstr = list(A = rbind( matrix(day_constraint, nrow = 1),
                                                              matrix(scale(num_obs), nrow = 1) ),
                                                   e = c(0,0)),
                                hyper = list(theta = list(prior = 'pc.prec', 
                                                          param = c(log(2)/2, 0.50))))
                          + f(daily_iid, 
                              model = 'iid',
                              values = num_obs, 
                              hyper = list(theta = list(prior = 'pc.prec', 
                                                        param = c(log(2)/1, 0.50))))
                          + f(Temperature_RW2, 
                              model = 'rw2', 
                              values = temp_bins,
                              scale.model = TRUE,
                              constr = FALSE, 
                              extraconstr = list(A = rbind( matrix(temp_constraint, nrow = 1), 
                                                            matrix(scale(1:length(temp_constraint)), nrow = 1) ),
                                                 e = c(0,0)), 
                              hyper = list(theta = list(prior = 'pc.prec', 
                                                        param = c(log(2)/2, 0.50)))) ) ),
        data = city_data, 
        family = 'poisson', 
        lincomb = make_lincombs_fixedeffect(),
        control.predictor = list( compute = TRUE ),
        control.fixed = list( prec = 1/1, 
                              prec.intercept = 1/1 ),
        control.inla = list( strategy = "simplified.laplace",
                             improved.simplified.laplace = TRUE,
                             restart = 5,
                             stencil = 9,
                             step.len = step_len,
                             step.factor = step_factor,
                             h = 0.001 ),
        control.compute = list( config = TRUE, 
                                dic = TRUE, 
                                waic = TRUE,
                                smtp = "taucs" ),
        silent = 2L ) # close INLA call 
} # close function definition 
