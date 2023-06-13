# Modify estimate_R_agg to just return reconstructed incidence matrix
recon_inc_matrix <- function(incid,
                               dt = 7L,
                               dt_out = 7L,
                               iter = 10L,
                               tol = 1e-6,
                               recon_opt = "naive",
                               config = make_config(), 
                               method = c("non_parametric_si", "parametric_si"),
                               grid = list(precision = 0.001,
                                           min = -1,
                                           max = 1)){ 
  
  
  config_out <- config 
  n_dt <- length(incid) # number of aggregations
  
  if (length(dt) == 1){
    T <- n_dt * dt
    config$t_start <- seq(from = dt + 1, to = T - (dt - 1), dt)
    config$t_end <- seq(from = min(config$t_start) + (dt - 1),to = T, dt)
  } else if (length(dt) == length(incid)){
    T <- sum(dt)
    config$t_start <- cumsum(c(dt[1] + 1, dt[2:length(dt[-1])]))
    config$t_end <- cumsum(c(config$t_start[1] + dt[2] - 1, dt[3:length(dt)]))
  } else { # vector of repeating aggregations
    T <- sum(rep(dt, length.out = n_dt))
    # reorder dt as R estimation starts on second aggregation window
    reo_dt_start <- c(dt[2:length(dt)], dt[1])
    reo_dt_end <- c(reo_dt_start[2:length(reo_dt_start)], reo_dt_start[1])
    config$t_start <- cumsum(c(dt[1] + 1, 
                               rep(reo_dt_start, length.out = n_dt - 2)))
    config$t_end <- cumsum(c(config$t_start[1] + reo_dt_start[1] - 1, 
                             rep(reo_dt_end, length.out = n_dt - 2)))
  }
  
  niter <- seq(1, iter, 1) 
  sim_inc <- matrix(NA, nrow = T, ncol = iter)
  
  for (i in seq_along(niter)){
    if (niter[i] == 1){
      # Initialisation of EM. Aggregated incidence split evenly:
      if (length(dt) == 1){
        dis <- incid / dt
        dis_inc <- rep(dis, each = dt)
        full_dt <- rep(dt, n_dt)
      } else if (length(dt) == n_dt){
        dis <- incid / dt
        dis_inc <- rep(dis, times = dt)
        full_dt <- dt
      } else {
        full_dt <- rep(dt, length.out = n_dt)
        dis <- incid / full_dt
        dis_inc <- rep(dis, times = full_dt)
      }
      
      
      # Estimate R
      R <- estimate_R(dis_inc, 
                      method = method,
                      config = config)
      
      message("Estimated R for iteration: ", i)
      Mean_R <- R$R$`Mean(R)`
      
      if (anyNA(Mean_R)){
        idx_na <- which(is.na(Mean_R))
        idx_reconstruct <- seq(min(R$R$t_start[-idx_na]), length(dis_inc))
        Mean_R <- Mean_R[!is.na(Mean_R)]
      } else {
        idx_reconstruct <- seq(min(R$R$t_start), length(dis_inc))
      }
      
      # Index for aggregation windows
      idx_aggregation <- rep(seq(1:n_dt), times=full_dt) 
      
      # Two options:
      # Opt 1) To keep naive disaggregation of the incidence for the aggregation 
      # window which precedes the first window that R can be estimated for:
      if (recon_opt == "naive"){
        aggs_to_reconstruct <- seq(idx_aggregation[idx_reconstruct[1]], n_dt)
      } 
      
      # Opt 2) To reconstruct the incidence in the preceding aggregation window by 
      # assuming that the growth rate matches that of the first estimation window:
      if (recon_opt == "match"){
        aggs_with_estimate <- seq(idx_aggregation[idx_reconstruct[1]], n_dt)
        aggs_to_reconstruct <- c(min(aggs_with_estimate) - 1, aggs_with_estimate)
        add_idx <- rev((min(idx_reconstruct) - 1) : (min(idx_reconstruct) - dt))
        idx_reconstruct <- c(add_idx, idx_reconstruct)
      }
      
      ## Incidence that can't be reconstructed (e.g. if recon_opt = "naive" and 
      # can't estimate R for first agg window or if incidence is too low)
      incid_not_to_reconstruct <- dis_inc[-idx_reconstruct]
      incid_to_reconstruct <- incid[aggs_to_reconstruct]
      
      # Translate R to growth rate
      get_r_from_R <- function(R, gt_mean, gt_sd, 
                               gt_distr,
                               grid) {
        r_grid <- seq(grid$min, grid$max, grid$precision)
        if (is.null(gt_distr)) {
          gt_pars <- gamma_mucv2shapescale(mu = gt_mean, cv = gt_sd / gt_mean)
          gt_distr <- distcrete("gamma", interval = 1,
                                shape = gt_pars$shape,
                                scale = gt_pars$scale, w = 0.5)
        }
        # using a grid of r values translate that into R using r2R0
        R_grid <- epitrix::r2R0(r = r_grid, w = gt_distr)
        
        # find location of the value in the R grid which has the smallest 
        # difference to the input of R the user provided e.g. R_grid[idx_r]:
        idx_r <- vapply(R, function(e) which.min(abs(R_grid - e)), numeric(1L)) 
        
        while (any(idx_r == 1) || any(idx_r == length(r_grid))) {
          # if necessary rerun get_r_from_R with a wider r_grid
          grid_multiplier <- 5
          if (grid$max > 0) grid$max <- grid_multiplier * grid$max else grid$max <- - grid$max
          if (grid$min < 0) grid$min <- grid_multiplier * grid$min else grid$min <- - grid$min
          r_grid <- seq(grid$min, grid$max, grid$precision)
          R_grid <- r2R0(r = r_grid, w = gt_distr)
          idx_r <- vapply(R, function(e) which.min(abs(R_grid - e)), numeric(1L))
        }
        r <- vapply(idx_r, function(e) r_grid[e], numeric(1L))
      }
      
      gr <- get_r_from_R(R = Mean_R, 
                         gt_mean = config$mean_si, gt_sd = config$std_si, 
                         gt_distr = config$si_distr,
                         grid = grid)
      
      # Assume the growth rates match to reconstruct preceding aggregation window:
      if (recon_opt == "match") {
        gr <- c(gr[1], gr)
      }
      
      # Estimate incidence using growth rate
      
      # Assume that It is a constant (k) multiplied by exp(gr[for that dt]*t)
      d <- numeric(length(aggs_to_reconstruct))
      k <- numeric(length(aggs_to_reconstruct))
      dt_seq <- full_dt[aggs_to_reconstruct]
      
      for (w in seq_along(k)){
        if (dt_seq[w] > 1){
          d[w] <- sum(exp(gr[w] * seq(1, dt_seq[w] - 1, 1)))
          k[w] <- incid_to_reconstruct[w] / (exp(gr[w]) * (1 + d[w]))
        } else { # if dt is 1 no need to reconstruct
          d[w] <- 0
          k[w] <- incid_to_reconstruct[w]
          gr[w] <- 0
        }
      }
      
      recon_df <- data.frame(k = k, gr = gr, dt = dt_seq)
      k_ls <- list()
      gr_ls <- list()
      
      for (f in seq_along(incid_to_reconstruct)){
        k_ls[[f]] <- rep(recon_df$k[f], recon_df$dt[f])
        gr_ls[[f]] <- rep(recon_df$gr[f], recon_df$dt[f])
      }
      k_seq <- unlist(k_ls)
      gr_seq <- unlist(gr_ls)
      
      recon_days <- seq(1, sum(dt_seq))
      w_day <- list()
      
      for (x in seq_along(dt_seq)){
        w_day[[x]] <- seq(1,dt_seq[x])
      }
      
      w_day <- unlist(w_day)
      
      # Reconstruct daily incidence
      est_inc <- rep(NA, length(recon_days))
      
      for (t in seq_along(recon_days)){
        est_inc[t] <- k_seq[t] * exp(gr_seq[t] * w_day[t])
      }
      
      ## For the incidence that can't be reconstructed (can't estimate R for
      ## first agg window or if incidence is too low), using the initial dis_inc
      sim_inc[,i] <- c(incid_not_to_reconstruct, est_inc)
      
      message("Reconstructed incidence for iteration: ", i)
      
      
    } else {
      
      # Use new adjusted incidence as starting point
      new_inc <- sim_inc[,i-1]
      
      # Re-Estimate R
      R <- estimate_R(new_inc,
                      method = method,
                      config = config)
      
      message("Estimated R for iteration: ", i)
      
      Mean_R <- R$R$`Mean(R)`
      
      if (anyNA(Mean_R)){
        idx_na <- which(is.na(Mean_R))
        idx_reconstruct <- seq(min(R$R$t_start[-idx_na]), length(dis_inc))
        Mean_R <- Mean_R[!is.na(Mean_R)]
      } else {
        idx_reconstruct <- seq(min(R$R$t_start), length(dis_inc))
      }
      
      
      # Index for aggregation windows
      idx_aggregation <- rep(seq(1:n_dt), times=full_dt) 
      
      if (recon_opt == "naive"){
        aggs_to_reconstruct <- seq(idx_aggregation[idx_reconstruct[1]], n_dt)
      } 
      
      if (recon_opt == "match"){
        aggs_with_estimate <- seq(idx_aggregation[idx_reconstruct[1]], n_dt)
        aggs_to_reconstruct <- c(min(aggs_with_estimate) - 1, aggs_with_estimate)
        add_idx <- rev((min(idx_reconstruct) - 1) : (min(idx_reconstruct) - dt))
        idx_reconstruct <- c(add_idx, idx_reconstruct)
      }
      
      incid_not_to_reconstruct <- dis_inc[-idx_reconstruct]
      incid_to_reconstruct <- incid[aggs_to_reconstruct]
      
      # Translate R to growth rate again
      gr <- get_r_from_R(R = Mean_R, 
                         gt_mean = config$mean_si, gt_sd = config$std_si, 
                         gt_distr = config$si_distr,
                         grid = grid)
      
      if (recon_opt == "match"){
        gr <- c(gr[1], gr)
      }
      
      # Estimate incidence
      d <- numeric(length(aggs_to_reconstruct))
      k <- numeric(length(aggs_to_reconstruct))
      dt_seq <- full_dt[aggs_to_reconstruct]
      
      for (w in seq_along(k)){
        if (dt_seq[w] > 1){
          d[w] <- sum(exp(gr[w] * seq(1, dt_seq[w] - 1, 1)))
          k[w] <- incid_to_reconstruct[w] / (exp(gr[w]) * (1 + d[w]))
        } else { # if dt is 1 no need to reconstruct
          d[w] <- 0
          k[w] <- incid_to_reconstruct[w]
          gr[w] <- 0
        }
      }
      
      recon_df <- data.frame(k = k, gr = gr, dt = dt_seq)
      k_ls <- list()
      gr_ls <- list()
      
      for (f in seq_along(incid_to_reconstruct)){
        k_ls[[f]] <- rep(recon_df$k[f], recon_df$dt[f])
        gr_ls[[f]] <- rep(recon_df$gr[f], recon_df$dt[f])
      }
      k_seq <- unlist(k_ls)
      gr_seq <- unlist(gr_ls)
      
      w_day <- list()
      for (x in seq_along(dt_seq)){
        w_day[[x]] <- seq(1 , dt_seq[x])
      }
      w_day <- unlist(w_day)
      
      recon_days <- seq(1, sum(dt_seq))    
      est_inc <- rep(NA, length(recon_days))
      for (t in seq_along(recon_days)){
        est_inc[t] <- k_seq[t] * exp(gr_seq[t] * w_day[t])
      }
      
      sim_inc[,i] <- c(incid_not_to_reconstruct, est_inc)
      
      # monitor progress:
      message("Reconstructed incidence for iteration: ", i)
      
      # Final estimate R starting on the first aggregation window 
      # that incidence was able to be reconstructed over
      
      if (is.null(config_out$t_start)) {
        config_out$t_start <- seq(from = min(R$R$t_start), 
                                  to = T - (dt_out - 1), 1)
      }
      if (is.null(config_out$t_end)) {
        config_out$t_end <- config_out$t_start + (dt_out - 1)
      }
      
      if (niter[i] == max(niter)){
        if (any(abs(sim_inc[,i] - sim_inc[,i-1]) > tol)){
          message("Reconstructed incidence has not converged within the set
                  tolerance. Please run again with greater number of iterations.")
        }
        R_out <- estimate_R(sim_inc[,i],
                            method = method,
                            config = config_out)
        message("R estimation starts on day ", R_out$R$t_start[1])
      }
      
    }
  }
  sim_inc
}