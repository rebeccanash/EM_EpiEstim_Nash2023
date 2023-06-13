## Mid-aggregation window variations in transmissibility

# Packages required
packages <- c("EpiEstim", "projections", "dplyr", "ggplot2", "hrbrthemes",
              "cowplot")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

######################
# Simulate incidence #
######################

# Simulate over 70 days (10 weeks)
n_days <- 70

# Seed epidemic with 7 days of 10 daily cases
initial_i <- incidence::incidence(c(rep(1, 10), rep(2, 10), rep(3, 10),
                                    rep(4, 10), rep(5, 10), rep(6, 10),
                                    rep(7, 10)))
initial_i$counts

# Serial interval (assume the same as influenza)
mean_si <- 3.6
sd_si <- 1.6
si <- discr_si(0:30, mu = mean_si, si = sd_si)
si_zero <- si / sum(si)
si_one <- si_zero[-1]
n_sim <- 100

proj <- list()
rt <- list()

# Increase in transmissibility on weekends
rt[[1]] <- c(rep(c(1.25, 1.5), 10))

# Genuine decrease in transmissibility on weekends
rt[[2]] <- c(rep(c(1.5, 1.25), 10))

# time change last day you want old R minus 1
s1 <- seq(4, 70, 7)
s2 <- seq(6, 70, 7)
idx_tc <- sort(c(s1, s2))
tc <- idx_tc[-length(idx_tc)]

for (i in seq_along(rt)){
  set.seed(5)
  proj[[i]] <- project(
    x = initial_i,
    R = rt[[i]],
    si = si_one,
    n_days = n_days,
    n_sim = n_sim,
    time_change = tc,
    instantaneous_R = TRUE)
}

proj_data <- list()
for (i in seq_along(rt)){
  proj_data[[i]] <- as.data.frame(proj[[i]], long = TRUE)
}

name_rt <- c("1.25to1.5", "1.5to1.25")
for(x in seq_along(rt)) {
  saveRDS(proj_data[[x]], file = paste0(
    "supplementary_analysis/mid_aggregation_variations/sims_rt_",
    name_rt[x], ".rds"))
}

####################
# Load simulations #
####################

sims <- list()
names <- name_rt

for (x in seq_along(names)) {
  sims[[x]] <- readRDS(
    paste0("supplementary_analysis/mid_aggregation_variations/sims_rt_",
           names[x], ".rds"))
}

n_results <- length(sims)
names(sims) <- names

# Serial interval
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si, std_si = sd_si)

# Create 2 matrices 'all_inc' with each column = simulation, and each row = day
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

# This takes a while, can just load results below

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
#   saveRDS(get(paste0("rt_", x)), file = paste0(
#     "supplementary_analysis/mid_aggregation_variations/result_rt_",
#     names(sims[x]), ".rds"))
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
#   saveRDS(get(paste0("rt_d_", x)), file = paste0(
#     "supplementary_analysis/mid_aggregation_variations/result_rt_",
#     names(sims[x]), "_daily.rds"))
# }


## Estimates using non-sliding time windows (daily R estimates)

# list_vec <- setNames(vector("list", n_results), paste0("rt_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_weekly[[res]][, x],
#                                  dt = 7L,
#                                  dt_out = 2L,
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
#           file = paste0("supplementary_analysis/mid_aggregation_variations/result_rt_",
#                         names(sims[x]), "_nonsliding.rds"))
# }
# 
# # Also estimate non-sliding Rt from simulated daily data as reference
# list_vec <- setNames(vector("list", n_results), paste0("rt_d_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# t_start <- seq(from = 2, to = ndays - 1, 1)
# config_daily <- EpiEstim::make_config(mean_si = mean_si,
#                                       std_si = sd_si,
#                                       t_start = t_start,
#                                       t_end = t_start + 1)
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = all_inc[[res]][, x],
#                                  method = "parametric_si",
#                                  config = config_daily)
#     assign(paste0("rt_d_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("rt_d_", x)),
#           file = paste0("supplementary_analysis/mid_aggregation_variations/result_rt_",
#                         names(sims[x]), "_nonsliding_daily.rds"))
# }

#################
## Load results #
#################

result <- list()
all_results <- n_results * 4
all_names <- c(names,
               paste(names, "daily", sep = "_"),
               paste(names, "nonsliding", sep = "_"),
               paste(names, "nonsliding_daily", sep = "_"))

for (n in seq_along(all_names)) {
  result[[n]] <- readRDS(
    paste0("supplementary_analysis/mid_aggregation_variations/result_rt_",
           all_names[n], ".rds"))
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
  all_ests[[res]] <- process_res(get(paste0("ests", res)))
}

###########
# R plots #
###########

r_plot <- function(data, legend_pos = "none", ymin = 0.7, ymax = 2,
                   xmin = 14, xmax = 70, legend_size = 11,
                   y_axis_lab = "Weekly sliding R estimate") {
  ggplot(data, aes(x = t_end,
                   color = source, fill = source, lty = source)) +
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
    ylab(y_axis_lab) +
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
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
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
  true <- rep(c(rep(true_all[[res]][1], 5), rep(true_all[[res]][2], 2)), 8)
  true <- c(true[1], true)
  assign(paste0("true", res), as.numeric(true))
}

dat1 <- gdata::combine(all_ests[[1]], all_ests[[3]],
                       add_true_r(true1),
                       names = c("weekly", "daily", "true"))
dat2 <- gdata::combine(all_ests[[2]], all_ests[[4]],
                       add_true_r(true2),
                       names = c("weekly", "daily", "true"))
dat3 <- gdata::combine(all_ests[[5]], all_ests[[7]],
                       add_true_r(true1),
                       names = c("weekly", "daily", "true"))
dat4 <- gdata::combine(all_ests[[6]], all_ests[[8]],
                       add_true_r(true2),
                       names = c("weekly", "daily", "true"))

# panels with legend or different y axis labels
weekly <- c(1, 2)

for (l in seq_along(weekly)) {
  plot_r <- r_plot(data = get(paste0("dat", weekly[l])),
                   legend_pos = "none", ymin = 1, ymax = 1.7)
  assign(paste0("plot_r_", weekly[l]), plot_r)
}

plot_r_3 <- r_plot(data = dat3,
                   legend_pos = c(0.8, 0.9), ymin = 0.7, ymax = 2.1,
                   y_axis_lab = "Daily R estimate")

plot_r_4 <- r_plot(data = dat4,
                   legend_pos = "none", ymin = 0.7, ymax = 2.1,
                   y_axis_lab = "Daily R estimate")

# Plot both
all <- plot_grid(plot_r_3, plot_r_4,
                 plot_r_1, plot_r_2,
                 ncol = 2,
                 labels = c("A", "B", "C", "D"),
                 label_size = 12)

ggsave(plot = all, "supplementary_figures/mid_window_variability.pdf",
       device = cairo_pdf, width = 12, height = 8)
