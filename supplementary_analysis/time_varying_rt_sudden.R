## Simulation study: time-varying Rt with sudden change

# Packages required
packages <- c("EpiEstim", "dplyr", "tidyr", "ggplot2",
              "hrbrthemes", "cowplot", "gdata")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

####################
# Load simulations #
####################

# Simulated incidence for sudden increase/decrease in Rt (simulate_incidence.R):
sims <- list()
names <- c("1to1.25", "1.25to1.5", "1.5to1.75",
           "1.25to0.75", "1.25to1", "1.5to1.25", "1.75to1.5")

for (x in seq_along(names)) {
sims[[x]] <- readRDS(
  paste0("supplementary_analysis/simulated_incidence/sims_rt_sud_",
         names[x], ".rds"))
}

n_results <- length(sims)
names(sims) <- names

# Serial interval
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si, std_si = sd_si)

# Create 7 matrices 'all_inc' with each column = simulation, and each row = day
ndays <- 70
sim_idx <- seq(1, 100, 1)
all_inc <- lapply(1:n_results, matrix,
                  data = NA, nrow = ndays, ncol = max(sim_idx))

split_start <- (sim_idx * ndays) - (ndays - 1)
split_end <- sim_idx * ndays

for (i in seq_along(sims)){
  for (x in seq_along(sim_idx)){
    inc <- sims[[i]]$incidence[split_start[x] : split_end[x]]
    all_inc[[i]][, x] <- inc
  }
}

# Visualise simulated incidence (100 epidemic trajectories for each)
par(mfrow=c(2, 4))
layout(matrix(c(1, 4, 2, 5, 3, 6, 0, 7), ncol = 4, nrow = 2))
for (i in seq_along(sims)){
  matplot(all_inc[[i]], type = "l", ylab = "Simulated Incidence")
}
dev.off()

# Aggregate daily simulated incidence into weekly incidence
nweeks <- 10
sims_weekly <- lapply(1:n_results, matrix, data = NA,
                      nrow = nweeks, ncol = max(sim_idx))

for (i in seq_along(sims_weekly)){
  for (x in seq_along(sim_idx)){
    sims_weekly[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][, x]))
  }
}

##################################
# Estimate R for each simulation #
##################################

#############################################################################
# WARNING: n = 700, so takes a while! You can run this code or you can just
# load the result of this below.
#############################################################################

# list_vec <- setNames(vector("list", n_results), paste0("rt_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_weekly[[res]][, x],
#                                  dt = 7L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("rt_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("rt_", x)),
#           file = paste0("supplementary_analysis/tv_rt_sudden/result_rt_sud_",
#                         names(sims[x]), ".rds"))
# }
# 
# # Also estimate Rt from simulated daily data as reference
# list_vec <- setNames(vector("list", n_results), paste0("rt_d_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = all_inc[[res]][, x],
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("rt_d_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("rt_d_", x)),
#           file = paste0("supplementary_analysis/tv_rt_sudden/result_rt_sud_",
#                         names(sims[x]), "_daily.rds"))
# }

#################
## Load results #
#################

result <- list()
all_results <- n_results * 2
daily_names <- paste(names, "daily", sep = "_")
all_names <- c(names, daily_names)

for (n in seq_along(all_names)) {
  result[[n]] <- readRDS(
    paste0("supplementary_analysis/tv_rt_sudden/result_rt_sud_",
           all_names[n], ".rds"))
}

# function to select t_end, mean, lower and upper quantiles for each result
select_res <- function(result) {
  result %>% select(t_end, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`)
}

# create ests1-14
list_vec <- setNames(vector("list", all_results),
                     paste0("ests", 1:all_results))
list2env(list_vec, .GlobalEnv)
est <- list()

# apply function for each result
for (res in seq_along(all_names)){
  for (i in seq_along(sim_idx)){
    est[[i]] <- select_res(result[[res]][i][[1]]$R)
    assign(paste0("ests", res), est)
  }
}

all_ests <- list()

# function to get means
process_res <- function(result, true_r) {
  bind_rows(result, .id = "id") %>%
    group_by(t_end) %>%
    summarise(mean_r = mean(`Mean(R)`),
              upper = mean(`Quantile.0.975(R)`),
              lower = mean(`Quantile.0.025(R)`))
}

for (res in seq_along(all_names)) {
  all_ests[[res]] <- process_res(get(paste0("ests", res)))
}

###########
# R plots #
###########

r_plot <- function(data, legend_pos = "none", ymin = 0.7, ymax = 2,
                   xmin = 14, xmax = 70, legend_size = 11) {
  ggplot(data, aes(x = t_end,
                   color = source, fill = source, lty = source)) +
    geom_hline(yintercept = 1, col = "grey", lty = 3) +
    geom_line(aes(y = mean_r)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), col = NA) +
    scale_color_manual("",
                       values = c("daily" = "darkgrey",
                                  "weekly" = "darkgreen",
                                  "true" = "black"),
                       labels = c("daily" = "Daily data",
                                  "weekly" = "Reconstructed data",
                                  "true" = "True R")) +
    scale_fill_manual("",
                      values = c("daily" = alpha("grey", 0.2),
                                 "weekly" = alpha("mediumseagreen", 0.2),
                                 "true" = alpha("white", 0)),
                      labels = c("daily" = "Daily data",
                                 "weekly" = "Reconstructed data",
                                 "true" = "True R")) +
    scale_linetype_manual("",
                          values = c("daily" = "solid",
                                     "weekly" = "solid",
                                     "true" = "dashed"),
                          labels = c("daily" = "Daily data",
                                     "weekly" = "Reconstructed data",
                                     "true" = "True R")) +
    ylab("Weekly sliding R estimate") +
    xlab("Time") +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    theme_ipsum() +
    theme(
      legend.position = legend_pos,
      legend.text = element_text(size = legend_size),
      legend.background = element_rect(fill="white", colour="white"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "lightgrey", linewidth = 0.3))
}

add_true_r <- function(true) {
  data.frame(t_end = all_ests[[1]]$t_end,
             mean_r = true,
             upper = rep(NA, length = nrow(all_ests[[1]])),
             lower = rep(NA, length = nrow(all_ests[[1]])))
}

true_all <- strsplit(names, "to")

for (res in seq_along(names)) {
  true <- c(rep(true_all[[res]][1], 22), rep(true_all[[res]][2], 35))
  assign(paste0("true", res), as.numeric(true))
}

for (res in seq_along(names)) {
  datx <- gdata::combine(all_ests[[res]], all_ests[[res + n_results]],
                        add_true_r(get(paste0("true", res))),
                        names = c("weekly", "daily", "true"))
  assign(paste0("dat",res), datx)
}

# panels with legend or no legend
no_legend <- c(1, 2, 4, 5, 6)
legend <- c(3, 7)

for (l in seq_along(no_legend)) {
  plot_r <- r_plot(data = get(paste0("dat", no_legend[l])),
                   legend_pos = "none", ymin = 0.5)
  assign(paste0("plot_r_", no_legend[l]), plot_r)
}

for (l in seq_along(legend)) {
  plot_r <- r_plot(data = get(paste0("dat", legend[l])),
                   legend_pos = c(0.80, 0.37), ymin = 0.5)
  assign(paste0("plot_r_", legend[l]), plot_r)
}

##############
# Bias plots #
##############

# Mean and SD of the bias (mean posterior estimate of Rt - true Rt)

process_bias <- function(dat) {
  t_end <- dat$t_end[dat$source == "weekly"]
  bias <- dat$mean_r[dat$source == "weekly"] - dat$mean_r[dat$source == "true"]
  upper_bias <- dat$upper[dat$source == "weekly"] - dat$mean_r[dat$source == "true"]
  lower_bias <- dat$lower[dat$source == "weekly"] - dat$mean_r[dat$source == "true"]
  data.frame(t_end, bias, upper_bias, lower_bias)
}

for (res in seq_along(names)) {
  bias <- process_bias(get(paste0("dat", res)))
  assign(paste0("bias", res), bias)
}

bias_plot <- function(data) {
  ggplot(data, aes(x = t_end, y = bias)) +
    geom_line(col = "darkgreen", linewidth = 0.5) +
    geom_vline(xintercept = 35, lty = 2, col = "lightgrey") +
    geom_ribbon(aes(ymin = lower_bias, ymax = upper_bias), col = NA,
                fill = "mediumseagreen", alpha = 0.3) +
    geom_hline(yintercept = 0, lty = 2) +
    ylim(-1, 1) +
    ylab("Bias") +
    xlab("Time") +
    theme_ipsum() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "lightgrey", linewidth = 0.3),
      axis.ticks.x = element_blank())
}

for (res in seq_along(names)) {
  bias <- bias_plot(get(paste0("bias", res)))
  assign(paste0("plot_bias_", res), bias)
}

#####################
# Uncertainty plots #
#####################

# Mean +/- the standard deviation of the uncertainty

process_sd <- function(result, true_r) {
  bind_rows(result, .id = "id") %>%
    group_by(t_end) %>%
    summarise(mean_r = mean(`Mean(R)`),
              upper = mean(`Quantile.0.975(R)`),
              lower = mean(`Quantile.0.025(R)`),
              sd_upper = sd(`Quantile.0.975(R)`),
              sd_lower = sd(`Quantile.0.025(R)`),
              mean_unc = upper - lower,
              sd_unc_upper = (upper + sd_upper) - (lower + sd_lower),
              sd_unc_lower = (upper - sd_upper) - (lower - sd_lower))
}

all_sds <- list()
for (res in seq_along(names)) {
  all_sds[[res]] <- process_sd(get(paste0("ests",res)))
}

unc_plot <- function(data) {
  ggplot(data, aes(x = t_end, y = mean_unc)) +
    geom_line(col = "#061d33", linewidth = 0.5) +
    geom_vline(xintercept = 35, lty = 2, col = "lightgrey") +
    ylim(0, 0.8) +
    ylab("Uncertainty") +
    xlab("Time") +
    theme_ipsum() +
    theme(
      legend.position = c(0.50, 0.87),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "lightgrey", linewidth = 0.3),
      axis.ticks.x = element_blank())
}

for (res in seq_along(names)) {
  unc <- unc_plot(all_sds[[res]])
  assign(paste0("plot_unc_", res), unc)
}

######################
# 95% Coverage plots #
######################

# Proportion of weekly estimates (n = 100 per simulation per time point) where
# the 95% CrI contains true R or R estimated from daily data

# function to get mean R from daily data and true R for each scenario
process_ref <- function(data) {
  mean_r_d <- data$mean_r[data$source == "daily"][7:63]
  true_r <- data$mean_r[data$source == "true"]
  t_end <- data$t_end[data$source == "true"]
  data.frame(mean_r_d, true_r, t_end)
}

# apply function to each scenario
ref_mx <- lapply(1:n_results, matrix, data = NA, nrow = 57, ncol = 3)
for (res in seq_along(names)) {
  ref_mx[[res]] <- process_ref(get(paste0("dat", res)))
}

# create matrix for coverage
time <- length(ref_mx[[1]]$mean_r_d)
sims <- max(sim_idx)
cov_mx <- lapply(1:all_results, matrix, data = NA, nrow = time, ncol = sims)
cov <- matrix(NA, nrow = time, ncol = all_results)

# function to work out whether 95% CrI encompasses R from daily data or true R)
process_cov <- function(dat, ests) {
  ifelse(between(dat,
                 ests$`Quantile.0.025(R)`,
                 ests$`Quantile.0.975(R)`) == "TRUE", 1, 0)
}

# apply function to each scenario
for (res in seq_along(names)) {
  for (i in seq_along(sim_idx)) {
    # 95% CrI encompasses mean R estimated from daily data (true = 1)
    cov_mx[[res]][, i] <- process_cov(dat = ref_mx[[res]]$mean_r_d,
                                      ests = get(paste0("ests", res))[[i]])
    # 95% CrI encompasses true R (true = 1)
    cov_mx[[res + n_results]][, i] <- process_cov(dat = ref_mx[[res]]$true_r,
                                          ests = get(paste0("ests", res))[[i]])
  }
}

# get proportions and add to dataframe for plots
for (x in 1:all_results) {
  cov[, x] <- rowSums(cov_mx[[x]]) / sims
  colnames(cov) <- c(paste("cov_d", 1:n_results, sep = ""),
                     paste("cov_t", 1:n_results, sep = ""))
  cov_dat <- data.frame(cov, t_end = ref_mx[[1]]$t_end)
}

# credible interval encompassed true Rt x% of the time:
perc_time <- cov_dat %>% select(starts_with("cov_t")) %>%
  summarise(across(where(is.numeric), mean))
round(perc_time[1:3], 2) # increasing [87 - 96%]
round(perc_time[4:7], 2) # decreasing [85 - 90%]

# excluding week after step change:
perc_time_excl <- cov_dat %>% filter(!t_end %in% (36:42)) %>%
  select(starts_with("cov_t")) %>%
  summarise(across(where(is.numeric), mean))
round(perc_time_excl[1:3], 2) # increasing [95 - 97%]
round(perc_time_excl[4:7], 2) # decreasing [95 - 97%]

cov_dat <- cov_dat %>% pivot_longer(cols = starts_with("cov"),
                                    names_to = "type",
                                    values_to = "coverage")
cov_dat$result <- rep(c(1:n_results), length = nrow(cov_dat))
cov_dat$ref <- rep(rep(c("daily", "true"), each = n_results), 57)


plot_coverage <- function(dat, legend_pos = "none", legend_size = 11) {
  ggplot(dat, aes(x = t_end, colour = ref)) +
    geom_point(aes(y = coverage), size = 0.75) +
    geom_line(aes(y = coverage), linewidth = 0.5) +
    geom_vline(xintercept = 35, lty = 2, col = "lightgrey") +
    scale_color_manual("",
                       values = c("daily" = "#1b82e6",
                                  "true" = "#061d33"),
                       labels = c("daily" = "Daily data R",
                                  "true" = "True R")) +
    geom_hline(yintercept = 0.95, lty = 2) +
    ylim(0, 1) +
    ylab("95% Coverage") +
    xlab("Time") +
    theme_ipsum() +
    theme(
      legend.position = legend_pos,
      legend.text = element_text(size = legend_size),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "lightgrey", linewidth = 0.3),
      axis.ticks.x = element_blank())
}

for (res in seq_along(names)){
  datx <- cov_dat %>% filter(result == res)
  assign(paste0("cov_dat", res), datx)
}

# panels with legend or no legend
no_legend <- c(1, 2, 4, 5, 6)
legend <- c(3, 7)

for (l in seq_along(no_legend)) {
  plot_cov <- plot_coverage(dat = get(paste0("cov_dat", no_legend[l])),
                            legend_pos = "none")
  assign(paste0("cov", no_legend[l]), plot_cov)
}

for (l in seq_along(legend)) {
  plot_cov <- plot_coverage(dat = get(paste0("cov_dat", legend[l])),
                            legend_pos = c(0.90, 0.30))
  assign(paste0("cov", legend[l]), plot_cov)
}

#################
# 16 panel plot #
#################

inc_plots <- plot_grid(plot_r_1, plot_r_2, plot_r_3,
                       plot_bias_1, plot_bias_2, plot_bias_3,
                       plot_unc_1, plot_unc_2, plot_unc_3,
                       cov1, cov2, cov3,
                       ncol = 3, align = "hv",
                       labels = c("A", "B", "C",
                                  "D", "E", "F",
                                  "G", "H", "I",
                                  "J", "K", "L"),
                       label_size = 12)

dec_plots <- plot_grid(plot_r_4, plot_r_5, plot_r_6, plot_r_7,
                       plot_bias_4, plot_bias_5, plot_bias_6, plot_bias_7,
                       plot_unc_4, plot_unc_5, plot_unc_6, plot_unc_7,
                       cov4, cov5, cov6, cov7,
                       ncol = 4, align = "hv",
                       labels = c("A", "B", "C", "D",
                                  "E", "F", "G", "H",
                                  "I", "J", "K", "L",
                                  "M", "N", "O", "P"),
                       label_size = 12)

ggsave(plot = dec_plots, "supplementary_figures/tv_rt_sudden_decrease.pdf",
       device = cairo_pdf, width = 12, height = 12)

ggsave(plot = inc_plots, "supplementary_figures/tv_rt_sudden_increase.pdf",
       device = cairo_pdf, width = 12, height = 12)


###################
# Mid-point plots #
###################

# For the "influence of Rt plotting time relative to time window used" section
# Change t_end to mid-point of window instead of end-point

# re-create ests1-14 for mid-time windows
list_vec <- setNames(vector("list", all_results),
                     paste0("mid_ests", 1:all_results))
list2env(list_vec, .GlobalEnv)
est <- list()

# apply function for each result
for (res in seq_along(all_names)){
  for (i in seq_along(sim_idx)){
    est[[i]] <- select_res(result[[res]][i][[1]]$R)
    est[[i]]$t_end <- est[[i]]$t_end - 3.5
    assign(paste0("mid_ests", res), est)
  }
}

mid_all_ests <- list()

# get means
for (res in seq_along(all_names)) {
  mid_all_ests[[res]] <- process_res(get(paste0("mid_ests", res)))
}

for (res in seq_along(names)) {
  true <- c(rep(true_all[[res]][1], 25), rep(true_all[[res]][2], 32))
  assign(paste0("mid_true", res), as.numeric(true))
}

add_true_r_mid <- function(true) {
  data.frame(t_end = mid_all_ests[[1]]$t_end,
             mean_r = true,
             upper = rep(NA, length = nrow(mid_all_ests[[1]])),
             lower = rep(NA, length = nrow(mid_all_ests[[1]])))
}

for (res in seq_along(names)) {
  datx <- gdata::combine(mid_all_ests[[res]], mid_all_ests[[res + n_results]],
                         add_true_r_mid(get(paste0("mid_true", res))),
                         names = c("weekly", "daily", "true"))
  assign(paste0("mid_dat", res), datx)
}

#####################
# Mid-point R plots #
#####################

increasing <- c(1:3)
decreasing <- c(4:7)

for (l in seq_along(increasing)) {
  plot_r <- r_plot(data = get(paste0("mid_dat", increasing[l])),
                   legend_pos = "none", xmin = 10.5)
  assign(paste0("mid_plot_r_", increasing[l]), plot_r)
}

for (l in seq_along(decreasing)) {
  plot_r <- r_plot(data = get(paste0("mid_dat", decreasing[l])),
                   legend_pos = "none", ymin = 0.5, ymax = 2, xmin = 10.5)
  assign(paste0("mid_plot_r_", decreasing[l]), plot_r)
}

mid_plot_r_legend <- r_plot(data = get(paste0("mid_dat", 1)),
                            legend_pos = c(0.50, 0.60),
                            legend_size = 11)

legend <- get_legend(mid_plot_r_legend)

############################
# Mid-point coverage plots #
############################

# apply process_ref to each scenario
mid_ref_mx <- lapply(1:n_results, matrix, data = NA, nrow = 57, ncol = 3)
for (res in seq_along(names)) {
  mid_ref_mx[[res]] <- process_ref(get(paste0("mid_dat", res)))
}

# create matrix for coverage
time <- length(mid_ref_mx[[1]]$mean_r_d)
sims <- max(sim_idx)
mid_cov_mx <- lapply(1:all_results, matrix, data = NA, nrow = time, ncol = sims)
mid_cov <- matrix(NA, nrow = time, ncol = all_results)

# apply process_cov function to each scenario
for (res in seq_along(names)) {
  for (i in seq_along(sim_idx)) {
    # 95% CrI encompasses mean R estimated from daily data (true = 1)
    mid_cov_mx[[res]][, i] <- process_cov(dat = mid_ref_mx[[res]]$mean_r_d,
                                          ests = get(paste0("mid_ests", res))[[i]])
    # 95% CrI encompasses true R (true = 1)
    mid_cov_mx[[res + n_results]][, i] <- 
      process_cov(dat = mid_ref_mx[[res]]$true_r,
                  ests = get(paste0("mid_ests", res))[[i]])
  }
}

# get proportions and add to dataframe for plots
for (x in 1:all_results) {
  mid_cov[, x] <- rowSums(mid_cov_mx[[x]]) / sims
  colnames(mid_cov) <- c(paste("mid_cov_d", 1:n_results, sep = ""),
                         paste("mid_cov_t", 1:n_results, sep = ""))
  mid_cov_dat <- data.frame(mid_cov, t_end = mid_ref_mx[[1]]$t_end)
}

# credible interval encompassed true Rt x% of the time:
mid_perc_time <- mid_cov_dat %>% select(starts_with("mid_cov_t")) %>%
  summarise(across(where(is.numeric), mean))
round(mid_perc_time[1:3], 2) # increasing [88 - 96%]
round(mid_perc_time[4:7], 2) # decreasing [85 - 93%]

mid_cov_dat <- mid_cov_dat %>% pivot_longer(cols = starts_with("mid_cov"),
                                            names_to = "type",
                                            values_to = "coverage")
mid_cov_dat$result <- rep(c(1:n_results), length = nrow(mid_cov_dat))
mid_cov_dat$ref <- rep(rep(c("daily", "true"), each = n_results), 57)

# plots
for (res in seq_along(names)) {
  datx <- mid_cov_dat %>% filter(result == res)
  assign(paste0("mid_cov_dat", res), datx)
}

comb_dat <- list()
for (i in seq_along(names)) {
  comb_dat[[i]] <- gdata::combine(get(paste0("cov_dat", i)),
                                  get(paste0("mid_cov_dat", i)),
                                  names = c("end", "mid"))
  comb_dat[[i]]$cat <- paste(comb_dat[[i]]$ref, comb_dat[[i]]$source, sep = "_")
  comb_dat[[i]] <- filter(comb_dat[[i]], cat != "daily_end")
  comb_dat[[i]]$cat <- factor(comb_dat[[i]]$cat,
                              levels = rev(c("daily_mid", "true_mid", "true_end")))
}

plot_coverage_comb <- function(dat, legend_pos = "none") {
  ggplot(dat, aes(x = t_end, colour = cat)) +
    geom_vline(xintercept = 35, lty = 2, col = "lightgrey")+
    geom_point(aes(y = coverage), size = 0.75) +
    geom_line(aes(y = coverage), linewidth = 0.5) +
    scale_color_manual("",
                       values = c("daily_mid" = "#1b82e6",
                                  "true_mid" = "#061d33",
                                  "true_end" = "grey"),
                       labels = c("daily_mid" = "Daily data R",
                                  "true_mid" = "True R (mid-window)",
                                  "true_end" = "True R (end of window)"),
                       guide = guide_legend(reverse = TRUE)) +
    geom_hline(yintercept = 0.95, lty = 2) +
    ylim(0, 1) +
    ylab("95% Coverage") +
    xlab("Time") +
    theme_ipsum() +
    theme(
      legend.position = legend_pos,
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                  size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "lightgrey", linewidth = 0.3),
      axis.ticks.x = element_blank())
}


for (l in seq_along(names)) {
  plot_cov <- plot_coverage_comb(dat = comb_dat[[l]],
                                 legend_pos = "none")
  assign(paste0("mid_cov", l), plot_cov)
}

############################
# Combined mid-point plots #
############################

cov_legend <- get_legend(plot_coverage_comb(dat = comb_dat[[1]],
                                            legend_pos = c(0.5, 0.6)))

legends <- plot_grid(legend, cov_legend, nrow=2, align = 'hv')

inc_plots <- plot_grid(mid_plot_r_1, mid_plot_r_2, mid_plot_r_3,
                       mid_cov1, mid_cov2, mid_cov3,
                       ncol = 3, align = "hv",
                       labels = c("A", "B", "C",
                                  "D", "E", "F"),
                       label_size = 12)
inc_plots_with_legend <- plot_grid(inc_plots, legends, ncol=2, nrow=1, rel_widths = c(3,1))

dec_plots <- plot_grid(mid_plot_r_4, mid_plot_r_5, mid_plot_r_6, mid_plot_r_7,
                       mid_cov4, mid_cov5, mid_cov6, mid_cov7,
                       ncol = 4, align = "hv",
                       labels = c("G", "H", "I", "J",
                                  "K", "L", "M", "N"),
                       label_size = 12)

all_plots <- plot_grid(inc_plots_with_legend, dec_plots, ncol = 1, align = "hv")

ggsave(plot = all_plots, "supplementary_figures/tv_rt_sudden_midwindow.pdf",
       device = cairo_pdf, width = 12, height = 11)
