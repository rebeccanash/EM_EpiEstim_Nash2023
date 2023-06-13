## Simulation study: constant Rt

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

# Simulated incidence for Rt of 1, 1.25, 1.5, and 1.75 (simulate_incidence.R)
sims <- list()
names <- c("1", "1.25", "1.5", "1.75")

for (x in seq_along(names)) {
  sims[[x]] <- readRDS(
    paste0("supplementary_analysis/simulated_incidence/sims_rt_",
           names[x], ".rds"))
}

n_results <- length(sims)
names(sims) <- names

# Serial interval
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si, std_si = sd_si)

# Create 4 matrices 'all_inc' with each column = simulation, and each row = day
ndays <- 70
sim_idx <- seq(1, 100, 1)
all_inc <- lapply(1:4, matrix, data = NA, nrow = ndays, ncol = max(sim_idx))

split_start <- (sim_idx * ndays) - (ndays - 1)
split_end <- sim_idx * ndays

for (i in seq_along(sims)){
  for (x in seq_along(sim_idx)){
  inc <- sims[[i]]$incidence[split_start[x] : split_end[x]]
  all_inc[[i]][, x] <- inc
  }
}

# Visualise simulated incidence (100 epidemic trajectories for each)
par(mfrow=c(1, n_results))
for (i in seq_along(sims)){
  matplot(all_inc[[i]], type = "l", ylab = "Simulated Incidence")
}
dev.off()

# Aggregate daily simulated incidence into weekly incidence
nweeks <- 10
sims_weekly <- lapply(1:4, matrix, data = NA,
                        nrow = nweeks, ncol = max(sim_idx))

for (i in seq_along(sims_weekly)){
  for (x in seq_along(sim_idx)){
    sims_weekly[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][, x]))
  }
}

##################################
# Estimate R for each simulation #
##################################

list_vec <- setNames(vector("list", 4), paste0("rt_", 1:4))
list2env(list_vec, .GlobalEnv)
rt_result <- list()

#############################################################################
# WARNING: n = 400, so takes a while! You can run this code or you can just
# load the result of this below.
#############################################################################

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
#           file = paste0("supplementary_analysis/constant_rt/result_rt_",
#                         names(sims[x]), ".rds"))
# }
# 
# # Also estimate Rt from simulated daily data as reference
# list_vec <- setNames(vector("list", 4), paste0("rt_d_", 1:4))
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
#           file = paste0("supplementary_analysis/constant_rt/result_rt_",
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
    paste0("supplementary_analysis/constant_rt/result_rt_", all_names[n], ".rds"))
}

# function to select t_end, mean, lower and upper quantiles for each result
select_res <- function(result) {
  result %>% select(t_end, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`)
}

# create ests1-8
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
  all_ests[[res]] <- process_res(get(paste0("ests",res)))
}

###########
# R plots #
###########

r_plot <- function(data, legend_pos = "none", legend_size = 11) {
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
    xlim(14, 70) +
    ylim(0.7, 2) +
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
             mean_r = rep(true, length = nrow(all_ests[[1]])),
             upper = rep(NA, length = nrow(all_ests[[1]])),
             lower = rep(NA, length = nrow(all_ests[[1]])))
}

for (res in seq_along(names)) {
  datx <- gdata::combine(all_ests[[res]], all_ests[[res + n_results]],
                         add_true_r(as.numeric(names[res])),
                         names = c("weekly", "daily", "true"))
  assign(paste0("dat",res), datx)
}

plot_r_1 <- r_plot(data = dat1, legend_pos = "none")
plot_r_2 <- r_plot(data = dat2, legend_pos = "none")
plot_r_3 <- r_plot(data = dat3, legend_pos = "none")
plot_r_4 <- r_plot(data = dat4, legend_pos = c(0.75, 0.37))

##############
# Bias plots #
##############

# Mean and SD of the bias (mean posterior estimate of Rt - true Rt)

process_bias <- function(dat) {
  t_end <- dat$t_end[dat$source == "weekly"]
  bias <- dat$mean_r[dat$source == "weekly"] -
    dat$mean_r[dat$source == "true"]
  upper_bias <- dat$upper[dat$source == "weekly"] -
    dat$mean_r[dat$source == "true"]
  lower_bias <- dat$lower[dat$source == "weekly"] -
    dat$mean_r[dat$source == "true"]
  data.frame(t_end, bias, upper_bias, lower_bias)
}

for (res in seq_along(names)) {
  bias <- process_bias(get(paste0("dat", res)))
  assign(paste0("bias", res), bias)
}

bias_plot <- function(data) {
  ggplot(data, aes(x = t_end, y = bias)) +
  geom_line(col = "darkgreen", linewidth = 0.5) +
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
  all_sds[[res]] <- process_sd(get(paste0("ests", res)))
}

unc_plot <- function(data) {
  ggplot(data, aes(x = t_end, y = mean_unc)) +
    geom_line(col = "#061d33", linewidth = 0.5) +
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

cov_dat <- cov_dat %>% pivot_longer(cols = starts_with("cov"),
                                    names_to = "type",
                                    values_to = "coverage")
cov_dat$result <- rep(c(1:n_results), length = nrow(cov_dat))
cov_dat$ref <- rep(rep(c("daily", "true"), each = n_results), 57)


plot_coverage <- function(dat, legend_pos = "none") {
  ggplot(dat, aes(x = t_end, colour = ref)) +
    geom_point(aes(y = coverage), size = 0.75) +
    geom_line(aes(y = coverage), linewidth = 0.5) +
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

for (res in seq_along(names)){
  datx <- cov_dat %>% filter(result == res)
  assign(paste0("cov_dat", res), datx)
}

cov1 <- plot_coverage(cov_dat1, legend_pos = c(0.8, 0.37))
cov2 <- plot_coverage(cov_dat2)
cov3 <- plot_coverage(cov_dat3)
cov4 <- plot_coverage(cov_dat4)

#################
# 16 panel plot #
#################

all_plots <- plot_grid(plot_r_1, plot_r_2, plot_r_3, plot_r_4,
                       plot_bias_1, plot_bias_2, plot_bias_3, plot_bias_4,
                       plot_unc_1, plot_unc_2, plot_unc_3, plot_unc_4,
                       cov1, cov2, cov3, cov4,
                       ncol = 4, align = "hv",
                       labels = c("A", "B", "C", "D",
                                  "E", "F", "G", "H",
                                  "I", "J", "K", "L",
                                  "M", "N", "O", "P"),
                       label_size = 12)

ggsave(plot = all_plots, "supplementary_figures/constant_rt.pdf",
       device = cairo_pdf, width = 12, height = 11)
