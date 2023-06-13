## Alternative method: LOESS

# packages
library(stats)
library(EpiEstim)
library(cowplot)
library(gdata)

# load data
daily_inc <- readRDS("real_data/daily_flu.rds")
weekly_inc <- aggregate_inc(daily_inc$Incidence, dt = 7L)

# results produced in main_analysis/flu.R
rep_result <- readRDS("supplementary_analysis/flu/r_rep_sliding_result.rds")
em_result <- readRDS("supplementary_analysis/flu/R_agg_sliding_result.rds")

recon_inc <- em_result$I
rep_inc <- daily_inc$Incidence
naive_inc <- weekly_inc / 7

# plot naive inc at end of time window or middle
opt1 <- seq(7, by = 7, to = length(weekly_inc) * 7)
opt2 <- seq(3, by = 7, to = length(weekly_inc) * 7)
day_opt1 <- seq(7, by = 1, to = max(opt1))
day_opt2 <- seq(3, by = 1, to = max(opt2))

x0_opt1 <- loess(naive_inc ~ opt1, degree = 0)
x1_opt1 <- loess(naive_inc ~ opt1, degree = 1)
x2_opt1 <- loess(naive_inc ~ opt1, degree = 2)
x0_opt2 <- loess(naive_inc ~ opt2, degree = 0)
x1_opt2 <- loess(naive_inc ~ opt2, degree = 1)
x2_opt2 <- loess(naive_inc ~ opt2, degree = 2)

y0_opt1 <- predict(x0_opt1, day_opt1)
y1_opt1 <- predict(x1_opt1, day_opt1)
y2_opt1 <- predict(x2_opt1, day_opt1)
y0_opt2 <- predict(x0_opt2, day_opt2)
y1_opt2 <- predict(x1_opt2, day_opt2)
y2_opt2 <- predict(x2_opt2, day_opt2)

dat <- data.frame(recon_inc, rep_inc)
dat_naive <- data.frame(opt1, opt2, naive_inc)
dat_opt1 <- data.frame(day_opt1, y0_opt1, y1_opt1, y2_opt1)
dat_opt2 <- data.frame(day_opt2, y0_opt2, y1_opt2, y2_opt2)


plot_opt1 <- ggplot(dat, aes(x = 1:35)) +
  geom_line(aes(y = rep_inc, col = "Reported daily")) +
  geom_line(aes(y = recon_inc, col = "EM Reconstruction")) +
  geom_point(data = dat_naive, aes(x = opt1, y = naive_inc), pch = 1) +
  geom_line(data = dat_opt1, aes(x = day_opt1, y = y0_opt1,
                                 col = "LOESS (d = 0)")) +
  geom_line(data = dat_opt1, aes(x = day_opt1, y = y1_opt1,
                                 col = "LOESS (d = 1)")) +
  geom_line(data = dat_opt1, aes(x = day_opt1, y = y2_opt1,
                                 col = "LOESS (d = 2)")) +
  scale_colour_manual("",
                      values = c("Reported daily" = "grey",
                                 "EM Reconstruction" = "black",
                                 "LOESS (d = 0)" = "red",
                                 "LOESS (d = 1)" = "dodgerblue",
                                 "LOESS (d = 2)" = "magenta")) +
  ylab("Incidence") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = c(0.42, 0.85),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3))

plot_opt2 <- ggplot(dat, aes(x = 1:35)) +
  geom_line(aes(y = rep_inc, col = "Reported daily")) +
  geom_line(aes(y = recon_inc, col = "EM Reconstruction")) +
  geom_point(data = dat_naive, aes(x = opt2, y = naive_inc), pch = 1) +
  geom_line(data = dat_opt2, aes(x = day_opt2, y = y0_opt2,
                                 col = "LOESS (d = 0)")) +
  geom_line(data = dat_opt2, aes(x = day_opt2, y = y1_opt2,
                                 col = "LOESS (d = 1)")) +
  geom_line(data = dat_opt2, aes(x = day_opt2, y = y2_opt2,
                                 col = "LOESS (d = 2)")) +
  scale_colour_manual("",
                      values = c("Reported daily" = "grey",
                                 "EM Reconstruction" = "black",
                                 "LOESS (d = 0)" = "red",
                                 "LOESS (d = 1)" = "dodgerblue",
                                 "LOESS (d = 2)" = "magenta")) +
  ylab("Incidence") +
  xlab("Day") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 20, r = 0, b = 0, l = 0)),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3))

loess_comp <- plot_grid(plot_opt1, plot_opt2, ncol = 1,
                        labels = c("A", "B"), label_size = 12)
ggsave(plot = loess_comp, "supplementary_figures/loess_reconstructions.pdf",
       device = cairo_pdf, width = 10, height = 7)

# Values for Table (aggregate full weeks of daily data)
y0_agg <- aggregate_inc(y0_opt1[2:length(y0_opt1)], dt = 7L)
y1_agg <- aggregate_inc(y1_opt1[2:length(y0_opt1)], dt = 7L)
y2_agg <- aggregate_inc(y2_opt1[2:length(y0_opt1)], dt = 7L)

y0_agg_2 <- aggregate_inc(y0_opt2[6:26], dt = 7L)
y1_agg_2 <- aggregate_inc(y1_opt2[6:26], dt = 7L)
y2_agg_2 <- aggregate_inc(y2_opt2[6:26], dt = 7L)

##############
# Estimate R #
##############

# Serial interval - Cowling 2009
method <- "parametric_si"
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si,
                                std_si = sd_si)

result_loess_y0_opt1 <- EpiEstim::estimate_R(incid = y0_opt1,
                                             method = method,
                                             config = config)
result_loess_y1_opt1 <- EpiEstim::estimate_R(incid = y1_opt1,
                                             method = method,
                                             config = config)
result_loess_y2_opt1 <- EpiEstim::estimate_R(incid = y2_opt1,
                                             method = method,
                                             config = config)

result_loess_y0_opt2 <- EpiEstim::estimate_R(incid = y0_opt2,
                                             method = method,
                                             config = config)
result_loess_y1_opt2 <- EpiEstim::estimate_R(incid = y1_opt2,
                                             method = method,
                                             config = config)
result_loess_y2_opt2 <- EpiEstim::estimate_R(incid = y2_opt2,
                                             method = method,
                                             config = config)

# R estimates
r_rep <- data.frame(R = rep_result$R$`Mean(R)`,
                    lower = rep_result$R$`Quantile.0.025(R)`,
                    upper = rep_result$R$`Quantile.0.975(R)`,
                    t_end = rep_result$R$t_end)
r_em <- data.frame(R = em_result$R$`Mean(R)`,
                   lower = em_result$R$`Quantile.0.025(R)`,
                   upper = em_result$R$`Quantile.0.975(R)`,
                   t_end = em_result$R$t_end)
# Option 1: Recon starts on day 7, first t_start is day 8, first t_end day 14
r_loess_y0_opt1 <- data.frame(R = result_loess_y0_opt1$R$`Mean(R)`,
                              lower =
                                result_loess_y0_opt1$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y0_opt1$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y0_opt1$R$t_end + 6)
r_loess_y1_opt1 <- data.frame(R = result_loess_y1_opt1$R$`Mean(R)`,
                              lower =
                                result_loess_y1_opt1$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y1_opt1$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y1_opt1$R$t_end + 6)
r_loess_y2_opt1 <- data.frame(R = result_loess_y2_opt1$R$`Mean(R)`,
                              lower =
                                result_loess_y2_opt1$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y2_opt1$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y2_opt1$R$t_end + 6)
# Option 2: Recon starts on day 3, first t_start is day 4, first t_end day 10
r_loess_y0_opt2 <- data.frame(R = result_loess_y0_opt2$R$`Mean(R)`,
                              lower =
                                result_loess_y0_opt2$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y0_opt2$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y0_opt2$R$t_end + 2)
r_loess_y1_opt2 <- data.frame(R = result_loess_y1_opt2$R$`Mean(R)`,
                              lower =
                                result_loess_y1_opt2$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y1_opt2$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y1_opt2$R$t_end + 2)
r_loess_y2_opt2 <- data.frame(R = result_loess_y2_opt2$R$`Mean(R)`,
                              lower =
                                result_loess_y2_opt2$R$`Quantile.0.025(R)`,
                              upper =
                                result_loess_y2_opt2$R$`Quantile.0.975(R)`,
                              t_end =
                                result_loess_y2_opt2$R$t_end + 2)
comb_dat1 <- combine(r_rep, r_em,
                    r_loess_y0_opt1, r_loess_y1_opt1, r_loess_y2_opt1)
comb_dat2 <- combine(r_rep, r_em,
                     r_loess_y0_opt2, r_loess_y1_opt2, r_loess_y2_opt2)


plot_r_opt1 <- ggplot(comb_dat1, aes(x = t_end, color = source,
                                     fill = source)) +
  geom_hline(yintercept = 1, lty = 2, col = "lightgrey") +
  geom_line(aes(y = R)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, col = NA) +
  scale_color_manual("",
                    values = c("r_rep" = "grey",
                               "r_em" = "black",
                               "r_loess_y0_opt1" = "red",
                               "r_loess_y1_opt1" = "dodgerblue",
                               "r_loess_y2_opt1" = "magenta"),
                    labels = c("Reported daily", "EM reconstruction",
                               "LOESS (d = 0)",
                               "LOESS (d = 1)",
                               "LOESS (d = 2)")) +
  scale_fill_manual("",
                    values = c("r_rep" = "grey",
                               "r_em" = "black",
                               "r_loess_y0_opt1" = "red",
                               "r_loess_y1_opt1" = "dodgerblue",
                               "r_loess_y2_opt1" = "magenta"),
                    labels = c("Reported daily", "EM reconstruction",
                               "LOESS (d = 0)",
                               "LOESS (d = 1)",
                               "LOESS (d = 2)")) +
  ylab("R estimate") +
  xlab("") +
  ylim(0.5, 2) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 20, r = 0, b = 0, l = 0)))
plot_r_opt1



# Panel B

plot_r_opt2 <- ggplot(comb_dat2, aes(x = t_end, color = source,
                                     fill = source)) +
  geom_hline(yintercept = 1, lty = 2, col = "lightgrey") +
  geom_line(aes(y = R)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, col = NA) +
  scale_color_manual("",
                     values = c("r_rep" = "grey",
                                "r_em" = "black",
                                "r_loess_y0_opt2" = "red",
                                "r_loess_y1_opt2" = "dodgerblue",
                                "r_loess_y2_opt2" = "magenta"),
                     labels = c("Reported daily", "EM reconstruction",
                                "LOESS (d = 0)",
                                "LOESS (d = 1)",
                                "LOESS (d = 2)")) +
  scale_fill_manual("",
                    values = c("r_rep" = "grey",
                               "r_em" = "black",
                               "r_loess_y0_opt2" = "red",
                               "r_loess_y1_opt2" = "dodgerblue",
                               "r_loess_y2_opt2" = "magenta"),
                    labels = c("Reported daily", "EM reconstruction",
                               "LOESS (d = 0)",
                               "LOESS (d = 1)",
                               "LOESS (d = 2)")) +
  ylab("R estimate") +
  xlab("Day") +
  ylim(0.5, 2) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 20, r = 0, b = 0, l = 0)))
plot_r_opt2


r_loess_comp <- plot_grid(plot_r_opt1, plot_r_opt2, ncol = 1,
                          labels = c("A", "B"), label_size = 12)
ggsave(plot = r_loess_comp, "supplementary_figures/loess_rt.pdf",
       device = cairo_pdf, width = 10, height = 7)


# Bias

bias_dat1 <- comb_dat1 %>% filter(t_end >= 15)
bias_dat1$t_end <- as.factor(bias_dat1$t_end)

bias_dat1 <- bias_dat1 %>%
  dplyr::group_by(t_end) %>%
  dplyr::mutate(diff = R - R[1],
                diff = replace(diff, row_number() == 1, NA)) %>%
  ungroup(t_end) %>%
  filter(source != "r_rep")

bias_dat1$t_end <- as.numeric(paste(bias_dat1$t_end))

bias_opt1 <- ggplot(bias_dat1, aes(x = t_end, color = source)) +
  geom_hline(yintercept = 0, lty = 2, col = "lightgrey") +
  geom_line(aes(y = diff)) +
  ylim(-0.6, 0.6) +
  theme_minimal() +
  ylab("R (reported data) - \n R (reconstructed data)") +
  xlab("") +
  scale_color_manual("",
                     values = c("r_em" = "black",
                                "r_loess_y0_opt1" = "red",
                                "r_loess_y1_opt1" = "dodgerblue",
                                "r_loess_y2_opt1" = "magenta"),
                     labels = c("EM reconstruction",
                                "LOESS (d = 0)",
                                "LOESS (d = 1)",
                                "LOESS (d = 2)")) +
  theme(legend.position = c(0.7, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 20, r = 0, b = 0, l = 0)))


bias_dat2 <- comb_dat2 %>% filter(t_end >= 14 & t_end <= 31)
bias_dat2$t_end <- as.factor(bias_dat2$t_end)

bias_dat2 <- bias_dat2 %>%
  dplyr::group_by(t_end) %>%
  dplyr::mutate(diff = R - R[1],
                diff = replace(diff, row_number() == 1, NA)) %>%
  ungroup(t_end) %>%
  filter(source != "r_rep")

bias_dat2$t_end <- as.numeric(paste(bias_dat2$t_end))

bias_opt2 <- ggplot(bias_dat2, aes(x = t_end, color = source)) +
  geom_hline(yintercept = 0, lty = 2, col = "lightgrey") +
  geom_line(aes(y = diff)) +
  ylim(-0.6, 0.6) +
  theme_minimal() +
  ylab("R (reported data) - \n R (reconstructed data)") +
  xlab("Day") +
  scale_color_manual("",
                     values = c("r_em" = "black",
                                "r_loess_y0_opt2" = "red",
                                "r_loess_y1_opt2" = "dodgerblue",
                                "r_loess_y2_opt2" = "magenta"),
                     labels = c("EM reconstruction",
                                "LOESS (d = 0)",
                                "LOESS (d = 1)",
                                "LOESS (d = 2)")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "lightgrey", linewidth = 0.3),
        axis.title.y = element_text(margin =
                                      margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =
                                      margin(t = 20, r = 0, b = 0, l = 0)))

p1 <- bias_opt1 + coord_cartesian(xlim = c(14, 35))
p2 <- bias_opt2 + coord_cartesian(xlim = c(14, 35))

r_loess_bias <- plot_grid(p1, p2, ncol = 1,
                          labels = c("A", "B"), label_size = 12)

ggsave(plot = r_loess_bias, "supplementary_figures/loess_rt_bias.pdf",
       device = cairo_pdf, width = 10, height = 7)
