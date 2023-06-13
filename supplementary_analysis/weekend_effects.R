## Weekend effects

# Packages required
packages <- c("EpiEstim", "dplyr", "ggplot2", "hrbrthemes",
              "cowplot")

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

# One example of simulated incidence from each scenario: constant rt = 1.5,
# stepwise decrease from 1.5 to 1.25, and gradual increase from 1.25 to 1.5

sims <- list()
sims[[1]] <- readRDS(
  "supplementary_analysis/simulated_incidence/sims_rt_1.5.rds")
sims[[2]] <- readRDS(
  "supplementary_analysis/simulated_incidence/sims_rt_sud_1.5to1.25.rds")
sims[[3]] <- readRDS(
  "supplementary_analysis/simulated_incidence/sims_rt_grad_1.25to1.5.rds")

names <- c("const", "sud", "grad")
n_results <- length(sims)
names(sims) <- names

# Create 3 matrices 'all_inc' with each column = simulation, and each row = day
ndays <- 70
sim_idx <- seq(1, 100, 1)
all_inc <- lapply(1:n_results, matrix,
                  data = NA, nrow = ndays, ncol = max(sim_idx))

split_start <- (sim_idx * ndays) - (ndays - 1)
split_end <- sim_idx * ndays

for (i in seq_along(sims)) {
  for (x in seq_along(sim_idx)){
    inc <- sims[[i]]$incidence[split_start[x] : split_end[x]]
    all_inc[[i]][, x] <- inc
  }
}

# Artificially move 80% of cases from Sat/Sun to Mon/Tue
mon <- seq(1, 70, by = 7)
tue <- seq(2, 70, by = 7)
sat <- seq(6, 70, by = 7)
sun <- seq(7, 70, by = 7)
mt_index <- sort(c(mon, tue), decreasing = FALSE)
ss_index <- sort(c(sat, sun), decreasing = FALSE)

for (i in seq_along(sims)) {
# Find 80% cases that occur on Sat/Sun
prop_weekend <- all_inc[[i]][ss_index,] * 0.8
# Remove 80% from Sat/Sun
all_inc[[i]][ss_index,] <- all_inc[[i]][ss_index,] - prop_weekend
# Re-allocate to Mon/Tue
all_inc[[i]][mt_index,] <- all_inc[[i]][mt_index,] + prop_weekend
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

#############################################################################
# WARNING: n = 300 for each, so takes a while! You can run this code or you 
# can just load the result of this below.
#############################################################################

# # Serial interval
# mean_si <- 3.6
# sd_si <- 1.6
# 
# ## Estimates using non-sliding time windows (daily R estimates)
# 
# config <- EpiEstim::make_config(mean_si = mean_si, std_si = sd_si)
# 
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
#                         names(sims[x]), "_nonsliding_daily.rds"))
# }
# 
# ## Estimates using weekly sliding time windows
# 
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
#                         names(sims[x]), "_sliding.rds"))
# }
# 
# # Also estimate weekly sliding Rt from simulated daily data as reference
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
#                         names(sims[x]), "_sliding_daily.rds"))
# }
#
## Estimates using weekly sliding time windows
# list_vec <- setNames(vector("list", n_results), paste0("rt_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_weekly[[res]][, x],
#                                  dt = 7L,
#                                  dt_out = 14L,
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
#                         names(sims[x]), "_twsliding.rds"))
# }
#
# # Also estimate non-sliding Rt from simulated daily data as reference
# list_vec <- setNames(vector("list", n_results), paste0("rt_d_", 1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# t_start <- seq(from = 2, to = ndays - 13, 1)
# config_daily <- EpiEstim::make_config(mean_si = mean_si,
#                                       std_si = sd_si,
#                                       t_start = t_start,
#                                       t_end = t_start + 13)
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
#           file = paste0("supplementary_analysis/weekend_effects/result_rt_",
#                         names(sims[x]), "_twsliding_daily.rds"))
# }

#################
## Load results #
#################

files <- list.files("supplementary_analysis/weekend_effects", pattern = "rds")

# name cleaning
names <- tools::file_path_sans_ext(files)
names <- gsub("result_rt_", "", names)
all_names <- names

result <- list()
all_results <- (n_results * 2) * 3 # non-sliding, sliding, two-weekly sliding ests

for (n in seq_along(all_names)) {
  result[[n]] <- readRDS(
    paste0("supplementary_analysis/weekend_effects/result_rt_",
           all_names[n], ".rds"))
}

names(result) <- all_names

# function to select t_end, mean, lower and upper quantiles for each result
select_res <- function(result) {
  result %>% select(t_end, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`)
}

# create ests1-18
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
process_res <- function(result) {
  bind_rows(result, .id = "id") %>%
    group_by(t_end) %>%
    summarise(mean_r = mean(`Mean(R)`),
              upper = mean(`Quantile.0.975(R)`),
              lower = mean(`Quantile.0.025(R)`))
}

for (res in seq_along(all_names)) {
  all_ests[[res]] <- process_res(get(paste0("ests", res)))
}

names(all_ests) <- all_names

###########
# R plots #
###########

r_plot <- function(data,
                   ymin = 0.5, ymax = 2,
                   xmin = 14, xmax = 70,
                   legend_pos = "none", legend_size = 11,
                   axis_title = "Weekly sliding R estimate",
                   hline = NULL) {
  ggplot(data, aes(x = t_end,
                   color = source, fill = source, lty = source)) +
    geom_hline(yintercept = hline, col = "grey", lty = 3) +
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
    ylab(axis_title) +
    xlab("Time") +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
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
      axis.line = element_line(color = "lightgrey", linewidth = 0.3))
}

add_true_r <- function(true) {
  data.frame(t_end = 1:70,
             mean_r = true,
             upper = rep(NA, length = 70),
             lower = rep(NA, length = 70))
}

# true R const
true_const <- rep(1.5, length = 70)

# true R sudden
true_sud <- c(rep(1.5, 35), rep(1.25, 35))

# true R grad
true_grad <- c(rep(1.25, 18), seq(from = 1.25, to = 1.5, length = 32), rep(1.5, 20))

scenario <- rep(c("const", "grad", "sud"), each = 6)
comb_idx <- seq(2, 18, 2)

for (res in comb_idx) {
  datx <- gdata::combine(all_ests[[res - 1]], all_ests[[res]],
                         add_true_r(get(paste0("true_", scenario[res]))),
                         names = c("daily", "weekly", "true"))
  datx <- datx[!(datx$source == "daily" &
                   datx$t_end < min(datx$t_end[datx$source == "weekly"])),]
  assign(paste0("dat", res/2), datx)
}

# panels for non-sliding, weekly, two weekly estimates
ns_plots <- c(1, 4, 7)
sliding_plots <- c(2, 5, 8)
tw_sliding_plots <- c(3, 6, 9)
legend <- 6

for (l in seq_along(sliding_plots)) {
  plot_r <- r_plot(data = get(paste0("dat", sliding_plots[l])),
                   ymin = 0.9, ymax = 2,
                   xmin = 0, xmax = 70,
                   legend_pos = "none")
  assign(paste0("plot_r_", sliding_plots[l]), plot_r)
}

for (l in seq_along(ns_plots)) {
  plot_r <- r_plot(data = get(paste0("dat", ns_plots[l])),
                   ymin = 0, ymax = 5,
                   xmin = 0, xmax = 70,
                   legend_pos = "none",
                   axis_title = "Daily R estimate")
  assign(paste0("plot_r_", ns_plots[l]), plot_r)
}

for (l in seq_along(tw_sliding_plots)) {
  plot_r <- r_plot(data = get(paste0("dat", tw_sliding_plots[l])),
                   ymin = 1, ymax = 2,
                   xmin = 0, xmax = 70,
                   legend_pos = "none",
                   axis_title = "Two week sliding R estimate")
  assign(paste0("plot_r_", tw_sliding_plots[l]), plot_r)
}

for (l in seq_along(legend)) {
  plot_r <- r_plot(data = get(paste0("dat", legend[l])),
                   ymin = 1, ymax = 2,
                   xmin = 0, xmax = 70,
                   legend_pos = c(0.35, 0.90),
                   axis_title = "Two week sliding R estimate")
  assign(paste0("plot_r_", legend[l]), plot_r)
}

###################
# Incidence plots #
###################

inc_ex <- matrix(data = NA, nrow = ndays, ncol = 6)

# constant example
inc_ex[, 1] <- result[[4]][[1]]$I
inc_ex[, 2] <- all_inc[[1]][, 1]

# sudden example
inc_ex[, 3] <- result[[16]][[1]]$I
inc_ex[, 4] <- all_inc[[2]][, 1]

# gradual example
inc_ex[, 5] <- result[[10]][[1]]$I
inc_ex[, 6] <- all_inc[[3]][, 1]
time <- seq(1,70,1)

inc_plots_dat <- data.frame(inc_ex, time)

inc_plot <- function(dat, sim_inc, recon_inc,
                     legend_pos = "none", legend_size = 11) {
  ggplot(dat, aes(x = time))+
  geom_line(aes(y = sim_inc, colour = "Simulated"), alpha = 0.8) +
  geom_line(aes(y = recon_inc, colour = "Reconstructed")) +
  scale_colour_manual("",
                      breaks = c("Simulated", "Reconstructed"),
                      values = c("darkgrey", "#297a4d")) +
  ylab("Daily incidence") + 
  xlab("Time") +
  theme_ipsum() +
  theme(
    legend.position = legend_pos,
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.ticks.x = element_blank())+
  coord_cartesian(clip = "off")+
  guides(colour = guide_legend(override.aes = list(alpha = c(0.8, 1))))
}

inc1 <- inc_plot(inc_plots_dat,
         sim_inc = inc_plots_dat$X2, recon_inc = inc_plots_dat$X1)
inc2 <- inc_plot(inc_plots_dat,
         sim_inc = inc_plots_dat$X4, recon_inc = inc_plots_dat$X3)
inc3 <- inc_plot(inc_plots_dat,
         sim_inc = inc_plots_dat$X6, recon_inc = inc_plots_dat$X5,
         legend_pos = c(0.30, 0.87))

############
# Plot all #
############

all <- plot_grid(inc1, inc2, inc3,
                 plot_r_1, plot_r_7, plot_r_4,
                 plot_r_2, plot_r_8, plot_r_5,
                 plot_r_3, plot_r_9, plot_r_6,
                 ncol=3, align='hv', labels=c("A","B","C",
                                              "D","E","F",
                                              "G","H","I",
                                              "J","K","L"), 
                 label_size = 12)

ggsave(plot = all,"supplementary_figures/weekend_effects.pdf",
       device = cairo_pdf, width = 12, height = 12)
