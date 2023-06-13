## Different temporal aggreations - time windows misalign with step change

# Packages required
packages <- c("EpiEstim", "dplyr", "ggplot2", "hrbrthemes", "cowplot")

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

##################
# Aggregate data #
##################

# t_start and t_end for each aggregation to misalign with step change
time <- list()
time[[1]] <- c(1, 69) # 3 day
time[[2]] <- c(3, 65) # 7 day
time[[3]] <- c(1, 70) # 10 day (already misalign)
time[[4]] <- c(1, 70) # 14 day (already misalign)

# 3-day aggregations misaligned
start <- time[[1]][1]
end <- time[[1]][2]
ndays <- length(start:end)

agg <- 3L
naggs <- ndays / agg

sims_agg_3 <- lapply(1:n_results, matrix, data = NA,
                     nrow = naggs, ncol = max(sim_idx))

for (i in seq_along(sims_agg_3)){
  for (x in seq_along(sim_idx)){
    sims_agg_3[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][start:end, x],
                                                 dt = agg))
  }
}
names(sims_agg_3) <- c("const", "sud", "grad")

# 7-day aggregations misaligned
start <- time[[2]][1]
end <- time[[2]][2]
ndays <- length(start:end)

agg <- 7L
naggs <- ndays / agg

sims_agg_7 <- lapply(1:n_results, matrix, data = NA,
                     nrow = naggs, ncol = max(sim_idx))

for (i in seq_along(sims_agg_7)){
  for (x in seq_along(sim_idx)){
    sims_agg_7[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][start:end, x],
                                                 dt = agg))
  }
}
names(sims_agg_7) <- c("const", "sud", "grad")

# 10-day aggregations misaligned
start <- time[[3]][1]
end <- time[[3]][2]
ndays <- length(start:end)

agg <- 10L
naggs <- ndays / agg

sims_agg_10 <- lapply(1:n_results, matrix, data = NA,
                      nrow = naggs, ncol = max(sim_idx))

for (i in seq_along(sims_agg_10)){
  for (x in seq_along(sim_idx)){
    sims_agg_10[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][start:end, x],
                                                  dt = agg))
  }
}
names(sims_agg_10) <- c("const", "sud", "grad")

# 14-day aggregations misaligned
start <- time[[4]][1]
end <- time[[4]][2]
ndays <- length(start:end)

agg <- 14L
naggs <- ndays / agg

sims_agg_14 <- lapply(1:n_results, matrix, data = NA,
                      nrow = naggs, ncol = max(sim_idx))

for (i in seq_along(sims_agg_14)){
  for (x in seq_along(sim_idx)){
    sims_agg_14[[i]][, x] <- matrix(aggregate_inc(all_inc[[i]][start:end, x],
                                                  dt = agg))
  }
}
names(sims_agg_14) <- c("const", "sud", "grad")

##############
# Estimate R #
##############

## This takes a long time, load results below

# mean_si <- 3.6
# sd_si <- 1.6
# config <- EpiEstim::make_config(mean_si = mean_si, std_si = sd_si)
# 
# # Weekly sliding estimates for 3-day
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_3day_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_3[[res]][, x],
#                                  dt = 3L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_3day_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_3day_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_3[x]), "_3day_misaligned.rds"))
# }
# 
# # Weekly sliding estimates for 7-day
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_7day_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_7[[res]][, x],
#                                  dt = 7L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_7day_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_7day_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_3[x]), "_7day_misaligned.rds"))
# }
# 
# # Weekly sliding estimates for 10-day
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_10day_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_10[[res]][, x],
#                                  dt = 10L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_10day_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_10day_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_10[x]), "_10day_misaligned.rds"))
# }
# 
# # Weekly sliding estimates for 14-day
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_14day_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_14[[res]][, x],
#                                  dt = 14L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_14day_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_14day_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_14[x]), "_14day_misaligned.rds"))
# }
# 
# 
# # Match sliding estimates time window to 10-day aggregations
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_10day_match_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_10[[res]][, x],
#                                  dt = 10L,
#                                  dt_out = 10L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_10day_match_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_10day_match_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_10[x]), "_10day_misaligned_match.rds"))
# }
# 
# # Match sliding estimates time window to 14-day aggregations
# 
# list_vec <- setNames(vector("list", n_results), paste0("misaligned_rt_14day_match_",
#                                                        1:n_results))
# list2env(list_vec, .GlobalEnv)
# rt_result <- list()
# 
# for (res in seq_along(names)){
#   for (x in seq_along(sim_idx)){
#     rt_result[[x]] <- estimate_R(incid = sims_agg_14[[res]][, x],
#                                  dt = 14L,
#                                  dt_out = 14L,
#                                  recon_opt = "match",
#                                  method = "parametric_si",
#                                  config = config)
#     assign(paste0("misaligned_rt_14day_match_", res), rt_result)
#   }
# }
# 
# # Save result
# for(x in seq_along(names)) {
#   saveRDS(get(paste0("misaligned_rt_14day_match_", x)),
#           file = paste0("supplementary_analysis/different_aggregations/result_rt_",
#                         names(sims_agg_14[x]), "_14day_misaligned_match.rds"))
# }

################
# Load results #
################

all_aggs <- c(3, 7, 10, 14)
all_aggs_match <- c(all_aggs, "10_match", "14_match")

for (a in seq_along(all_aggs_match)) {
  assign(paste0("result", all_aggs_match[a]), list()) 
}

# 3 day
for (n in seq_along(names)) {
  result3[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_3day_misaligned.rds"))
}

# 7 day
for (n in seq_along(names)) {
  result7[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_7day_misaligned.rds"))
}

# 10 day
for (n in seq_along(names)) {
  result10[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_10day_misaligned.rds"))
}

for (n in seq_along(names)) {
  result10_match[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_10day_misaligned_match.rds"))
}

# 14 day
for (n in seq_along(names)) {
  result14[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_14day_misaligned.rds"))
}

for (n in seq_along(names)) {
  result14_match[[n]] <- readRDS(
    paste0("supplementary_analysis/different_aggregations/result_rt_",
           names[n], "_14day_misaligned_match.rds"))
}

# function to select t_end, mean, lower and upper quantiles for each result
select_res <- function(result) {
  result %>% select(t_end, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`)
}

# apply function for each result
est <- list()
for (agg in seq_along(all_aggs_match)) {
  for (res in seq_along(names)) {
    for (i in seq_along(sim_idx)) {
      est[[i]] <- select_res(get(
        paste0("result", all_aggs_match[agg]))[[res]][i][[1]]$R)
      assign(paste0("ests", all_aggs_match[agg], "_", res), est)
    }
  }
}

# function to get means
process_res <- function(result) {
  bind_rows(result, .id = "id") %>%
    group_by(t_end) %>%
    summarise(mean_r = mean(`Mean(R)`),
              upper = mean(`Quantile.0.975(R)`),
              lower = mean(`Quantile.0.025(R)`))
}

n_ests <- seq(1, length(all_aggs_match) * 3)
n_res <- rep(seq(1:3), 6)
seq_aggs <- rep(all_aggs_match, each = 3)
all_ests <- list()

for(n in seq_along(n_ests)) {
  all_ests[[n]] <- process_res(get(
    paste0("ests", seq_aggs[n], "_", n_res[n])))
}


#########
# Plots #
#########

r_plot <- function(data,
                   ymin = 1, ymax = 2,
                   xmin = 1, xmax = 70,
                   legend_pos = "none", legend_size = 11,
                   agg_guide = NULL,
                   axis_title = "Weekly sliding R estimate",
                   hline = NULL) {
  ggplot(data, aes(x = t_end,
                   color = source, fill = source, lty = source)) +
    geom_hline(yintercept = hline, col = "grey", lty = 3) +
    geom_line(aes(y = mean_r)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), col = NA) +
    geom_vline(xintercept = agg_guide, lty = 3, col = "darkgrey")+
    scale_color_manual("",
                       values = c("weekly" = "darkgreen",
                                  "true" = "black"),
                       labels = c("weekly" = "Reconstructed data",
                                  "true" = "True R")) +
    scale_fill_manual("",
                      values = c("weekly" = alpha("mediumseagreen", 0.2),
                                 "true" = alpha("white", 0)),
                      labels = c("weekly" = "Reconstructed data",
                                 "true" = "True R")) +
    scale_linetype_manual("",
                          values = c("weekly" = "solid",
                                     "true" = "dashed"),
                          labels = c("weekly" = "Reconstructed data",
                                     "true" = "True R")) +
    ylab(axis_title) +
    xlab("Time") +
    xlim(xmin, xmax) +
    coord_cartesian(ylim = c(ymin, ymax)) +
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

r_plot_match <- function(data,
                         ymin = 1, ymax = 2,
                         xmin = 1, xmax = 70,
                         legend_pos = "none", legend_size = 11,
                         agg_guide = NULL,
                         axis_title = "Sliding R estimate",
                         legend_label = "10 day window",
                         hline = NULL) {
  ggplot(data, aes(x = t_end,
                   color = source, fill = source, lty = source)) +
    geom_hline(yintercept = hline, col = "grey", lty = 3) +
    geom_line(aes(y = mean_r)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), col = NA) +
    geom_vline(xintercept = agg_guide, lty = 3, col = "darkgrey")+
    scale_color_manual("",
                       values = c("weekly" = "darkgreen",
                                  "weekly_match" = "navy",
                                  "true" = "black"),
                       labels = c("weekly" = "7 day window",
                                  "weekly_match" = legend_label,
                                  "true" = "True R")) +
    scale_fill_manual("",
                      values = c("weekly" = alpha("mediumseagreen", 0.2),
                                 "weekly_match" = alpha("dodgerblue", 0.2),
                                 "true" = alpha("white", 0)),
                      labels = c("weekly" = "7 day window",
                                 "weekly_match" = legend_label,
                                 "true" = "True R")) +
    scale_linetype_manual("",
                          values = c("weekly" = "solid",
                                     "weekly_match" = "solid",
                                     "true" = "dashed"),
                          labels = c("weekly" = "7 day window",
                                     "weekly_match" = legend_label,
                                     "true" = "True R")) +
    ylab(axis_title) +
    xlab("Time") +
    xlim(xmin, xmax) +
    coord_cartesian(ylim = c(ymin, ymax)) +
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


# true R const
true_const <- rep(1.5, length = 70)

# true R sudden
true_sud <- c(rep(1.5, 35), rep(1.25, 35))

# true R grad
true_grad <- c(rep(1.25, 18), seq(from = 1.25, to = 1.5, length = 32), rep(1.5, 20))

add_true_r <- function(true) {
  data.frame(t_end = 1:70,
             mean_r = true,
             upper = rep(NA, length = length(true)),
             lower = rep(NA, length = length(true)))
}

seq_res <- rep(names, 6)
no_match <- seq_aggs[1:6]
match <- seq_aggs[7:12]

for (res in seq_along(no_match)) {
  datx <- gdata::combine(all_ests[[res]],
                         add_true_r(get(paste0("true_", seq_res[res]))),
                         names = c("weekly", "true"))
  assign(paste0("dat", res), datx)
}

for (res in seq_along(match)) {
  datx <- gdata::combine(all_ests[[res + 6]], all_ests[[res + 12]],
                         add_true_r(get(paste0("true_", seq_res[res]))),
                         names = c("weekly", "weekly_match", "true"))
  assign(paste0("dat", res + 6), datx)
}

# adjust t_end so match misaligning
# 7-day (removed first two days from inc)
for (i in 4:6) {
  adj <- get(paste0("dat", i)) %>% 
    mutate(t_end = ifelse(source == "weekly", t_end + 2, t_end))
  assign(paste0("dat", i), adj)
}

# grey polygons to show which times were not included to misalign time windows
excl_b <- 0
excl_e <- c(69:70, rev(69:70))
excl_yb <- 0
excl_ye <- c(c(rep(0.95, 2)), c(rep(2.25, 2)))
excl_3 <- data.frame(excl_b, excl_e, excl_yb, excl_ye)

excl_b <- c(1:3, rev(1:3))
excl_e <- c(65:70, rev(65:70))
excl_yb <- c(c(rep(0.95, 3)), c(rep(2.25, 3)))
excl_ye <- c(c(rep(0.95, 6)), c(rep(2.25, 6)))
excl_7 <- data.frame(excl_b, excl_e, excl_yb, excl_ye)

# for dashed aggregration lines
guide_3 <- seq(3, 69, 3)
guide_7 <- seq(9, 65, 7)
guide_10 <- seq(10, 70, 10)
guide_14 <- seq(14, 70, 14)

rep_agg <- rep(c(3, 7), each = 3)

for (p in seq_along(no_match)) {
  plot <- r_plot(data = get(paste0("dat", p)),
                 ymin = 1, ymax = 2,
                 legend_pos = "none",
                 agg_guide = get(paste0("guide_", rep_agg[p])),
                 axis_title = "Weekly sliding R estimate",
                 hline = NULL)
    plot <- plot +
      geom_polygon(data = get(paste0("excl_", rep_agg[p])), aes(y = excl_yb, x = excl_b),
                   fill = "grey96", alpha = 1, inherit.aes = FALSE) +
      geom_polygon(data = get(paste0("excl_", rep_agg[p])), aes(y = excl_ye, x = excl_e),
                   fill = "grey96", alpha = 1, inherit.aes = FALSE)
  assign(paste0("plot_r", p), plot)
}

rep_agg_match <- rep(c(10, 14), each = 3)

for (p in seq_along(match)) {
  plot <- r_plot_match(data = get(paste0("dat", p + 6)),
                       ymin = 1, ymax = 2,
                       legend_pos = "none",
                       agg_guide = get(paste0("guide_", rep_agg_match[p])),
                       axis_title = "Sliding R estimate",
                       legend_label = paste(rep_agg_match[p], "day window"),
                       hline = NULL)
  assign(paste0("plot_r", p + 6), plot)
}

# add legends
plot_r7 <- r_plot_match(data = dat7,
                        ymin = 1, ymax = 2,
                        legend_pos = c(0.75, 0.87),
                        agg_guide = guide_10,
                        axis_title = "Sliding R estimate",
                        legend_label = "10 day window",
                        hline = NULL)

plot_r10 <- r_plot_match(data = dat10,
                         ymin = 1, ymax = 2,
                         legend_pos = c(0.75, 0.87),
                         agg_guide = guide_14,
                         axis_title = "Sliding R estimate",
                         legend_label = "14 day window",
                         hline = NULL)

# all together

all <- plot_grid(plot_r1, plot_r2, plot_r3,
                 plot_r4, plot_r5, plot_r6,
                 plot_r7, plot_r8, plot_r9,
                 plot_r10, plot_r11, plot_r12,
                 ncol = 3,
                 labels=c("A","B","C",
                          "D","E","F",
                          "G","H","I",
                          "J","K","L"), 
                 label_size = 12)

ggsave(plot = all,"supplementary_figures/different_aggregations_misaligned.pdf",
       device = cairo_pdf, width = 12, height = 12)
