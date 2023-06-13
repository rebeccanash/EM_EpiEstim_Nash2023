## US Influenza analysis

# Packages required
packages <- c("EpiEstim", "scales", "dplyr", "ggplot2",
              "hrbrthemes", "cowplot", "grDevices")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

#############
# Load data #
#############

# Load daily influenza data
daily_inc <- readRDS("real_data/daily_flu.rds")
ndays <- length(daily_inc$Incidence)
index <- seq(1, ndays)

# Aggregate to weekly
weekly_inc <- EpiEstim::aggregate_inc(daily_inc$Incidence, dt = 7L)
saveRDS(weekly_inc, "supplementary_analysis/flu/weekly_inc_flu.rds")

# Serial interval - Cowling 2009
method <- "parametric_si"
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si,
                                std_si = sd_si)
t_start <- seq(from = 2, to = ndays - 1, 1)
config_daily <- EpiEstim::make_config(mean_si = mean_si,
                                      std_si = sd_si,
                                      t_start = t_start,
                                      t_end = t_start + 1)

# Estimate R from reported daily incidence (weekly sliding windows)
r_rep_sliding <- EpiEstim::estimate_R(incid = daily_inc$Incidence,
                                        method = method,
                                        config = config)
# save to use in supplementary analyses:
saveRDS(r_rep_sliding, "supplementary_analysis/flu/r_rep_sliding_result.rds")

# Estimate R from aggregated incidence (weekly sliding windows)
r_agg_sliding <- EpiEstim::estimate_R(incid = weekly_inc,
                                 dt = 7L,
                                 dt_out = 7L,
                                 recon_opt = "match",
                                 method = method,
                                 config = config)
saveRDS(r_agg_sliding, "supplementary_analysis/flu/r_agg_sliding_result.rds")

# Estimate R from reported daily incidence (non sliding - daily R estimate)
r_rep <- EpiEstim::estimate_R(incid = daily_inc$Incidence,
                                        method = method,
                                        config = config_daily)
saveRDS(r_rep, "supplementary_analysis/flu/r_rep_result.rds")

# Estimate R from aggregated incidence (no sliding windows - daily R estimate)
r_agg <- EpiEstim::estimate_R(incid = weekly_inc,
                              dt = 7L,
                              dt_out = 2L,
                              recon_opt = "match",
                              method = method,
                              config = config)
saveRDS(r_agg, "supplementary_analysis/flu/r_agg_result.rds")

recon_inc <- r_agg_sliding$I
inc_dat <- data.frame(index, daily_inc, recon_inc)
day1 <- min(r_agg$R$t_end) # day 9
day2 <- min(r_agg_sliding$R$t_end) # day 14
int <- length(seq(day2, ndays)) / 3

inc_dat$interval <- c(rep("no_est", length = day1 - 1),
                          rep("early_daily", length = day2 - day1),
                          rep(c("early", "mid", "late"), each = int), "late")

########################
## PANEL A: Incidence ##
########################

# Strip to show epidemic phase
early_days <- which(inc_dat$interval == "early")
mid_days <- which(inc_dat$interval == "mid")
late_days <- which(inc_dat$interval == "late")

early_y <- c(rep(-20, round(int) + 1), rep(50, round(int) + 1))
early_x <- c(c(early_days, max(early_days) + 1),
             rev(c(early_days, max(early_days) + 1)))
early_dat <- data.frame(early_y, early_x)

mid_y <- c(rep(-20, round(int) + 1), rep(50, round(int) + 1))
mid_x <- c(c(mid_days, max(mid_days) + 1), rev(c(mid_days, max(mid_days) + 1)))
mid_dat <- data.frame(mid_y, mid_x)

late_y <- c(rep(-20, round(int) + 1), rep(50, round(int) + 1))
late_x <- c(late_days, rev(late_days))
late_dat <- data.frame(late_y, late_x)

x_dates <- seq(as.Date("2009-12-11"), as.Date("2010-01-14"), by = "1 day")
date_1 <- format(x_dates[10], format = "%d %b %y")
date_2 <- format(x_dates[20], format = "%d %b %y")
date_3 <- format(x_dates[30], format = "%d %b %y")

# Plot

plot_inc <- ggplot(inc_dat, aes(x = index)) +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = c(10, 20, 30),
                     labels = c(date_1, date_2, date_3)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  geom_polygon(data = early_dat, aes(y = early_y, x = early_x),
               fill = "#B3E5FC", alpha = 0.8) +
  geom_polygon(data = mid_dat, aes(y = mid_y, x = mid_x),
               fill = "dodgerblue", alpha = 0.7) +
  geom_polygon(data = late_dat, aes(y = late_y, x = late_x),
               fill = "#01579B", alpha = 0.7) +
  geom_line(aes(y = Incidence, colour = "Reported"), alpha = 0.8) +
  geom_line(aes(y = recon_inc, colour = "Reconstructed")) +
  scale_colour_manual("",
                      values = c("Reported" = "darkgrey",
                                 "Reconstructed" = "#297a4d")) +
  ylab("Daily incidence by presentation at clinic") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.45, 0.87),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.05)) +
  coord_cartesian(clip = "off") +
  guides(colour = guide_legend(override.aes = list(alpha = c(0.8, 1))))
plot_inc

###################################
##  PANEL B: Squared error plot  ##
###################################

mod1 <- lm(r_agg_sliding$R$`Mean(R)` ~ r_rep_sliding$R$`Mean(R)`[7:28])
res1 <- signif(residuals(mod1), 5)
sq_err_sliding <- as.numeric(res1^2)
sliding_index <- r_agg_sliding$R$t_end

mod2 <- lm(r_agg$R$`Mean(R)` ~ r_rep$R$`Mean(R)`[6:32])
res2 <- signif(residuals(mod2), 5)
sq_err_daily <- as.numeric(res2^2)
daily_index <- r_agg$R$t_end

plot_diff_sliding <- data.frame(sliding_index, sq_err_sliding) %>%
  dplyr::rename(day = sliding_index)
plot_diff_daily <- data.frame(daily_index, sq_err_daily) %>%
  dplyr::rename(day = daily_index)
plot_diff <- plot_diff_sliding %>%
  dplyr::right_join(plot_diff_daily, by = c("day"))

# Plot

squared_err <- ggplot(plot_diff, aes(x = day)) +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = c(10, 20, 30),
                     labels = c(date_1, date_2, date_3)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  geom_line(aes(y = sq_err_sliding, colour = "Weekly sliding R")) +
  geom_line(aes(y = sq_err_daily, colour = "Daily R"), alpha = 0.3) +
  scale_colour_manual("",
                      breaks = c("Daily R", "Weekly sliding R"),
                      values = c(alpha("#ff6600", 0.3), "#e35874")) +
  ylab(expression(atop("Squared error of R estimates",
                       paste("using reported and reconstructed data")))) +
  xlab("Date") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.4, 0.87),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 10, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.x = element_line(size = 0.05),
    axis.line = element_line(color = "lightgrey", size = 0.3))
squared_err

######################################################
##  PANEL C: Difference in mean R (weekly sliding)  ##
######################################################

scatter_dat1 <- data.frame(agg = r_agg_sliding$R$`Mean(R)`,
                          rep = r_rep_sliding$R$`Mean(R)`[7:28],
                          agg_upper =
                            r_agg_sliding$R$`Quantile.0.975(R)`,
                          agg_lower =
                            r_agg_sliding$R$`Quantile.0.025(R)`,
                          rep_upper =
                            r_rep_sliding$R$`Quantile.0.975(R)`[7:28],
                          rep_lower =
                            r_rep_sliding$R$`Quantile.0.025(R)`[7:28],
                          col = inc_dat$interval[14:35])

corr_coeff <- cor(scatter_dat1$agg, scatter_dat1$rep)
r2 <- corr_coeff^2
round(r2, 2) # 0.91

corr1 <- ggplot(scatter_dat1, aes(x = rep, y = agg)) +
  geom_errorbar(aes(ymin = agg_lower, ymax = agg_upper, col = col),
                alpha = 1, size = 0.1) +
  geom_errorbarh(aes(xmin = rep_lower, xmax = rep_upper, col = col),
                 alpha = 1, size = 0.1) +
  geom_point(aes(col = col), size = 1.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  geom_hline(yintercept = 1, lty = 3) +
  geom_vline(xintercept = 1, lty = 3) +
  scale_colour_manual("",
                      breaks = c("early", "mid", "late"),
                      values = c("#B3E5FC", "dodgerblue", "#01579B")) +
  annotate("text", x = 1.3, y = 0.8, label = "R^2==0.91", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ylab("Mean weekly R (reconstructed data)") +
  xlab("Mean weekly R (reported data)") +
  ylim(0.6, 1.5) +
  xlim(0.6, 1.5) +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.20, 0.87),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3))
corr1

##########################
## PANEL D: Uncertainty ##
##########################

# Work out width of CrI's
scatter_dat1$width_weekly <- scatter_dat1$agg_upper - scatter_dat1$agg_lower
scatter_dat1$width_daily <- scatter_dat1$rep_upper - scatter_dat1$rep_lower

corr_coeff <- cor(scatter_dat1$width_weekly,
                  scatter_dat1$width_daily)
r2 <- corr_coeff^2
round(r2, 2) # 0.93

uncertainty1 <- ggplot(scatter_dat1, aes(x = width_daily, y = width_weekly)) +
  geom_point(aes(col = col), size = 1.5) +
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  scale_colour_manual("Time",
                      breaks = c("early", "mid", "late"),
                      values = c("#B3E5FC", "dodgerblue", "#01579B"),
                      guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(n = 3), limits = c(0.028, 0.049)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3), limits = c(0.028, 0.049)) +
  annotate("text", x = 0.044, y = 0.0315, label  = "R^2==0.93", size = 3.5,
           parse = TRUE) +
  ylab("95% CrI width for weekly R (reconstructed data)") +
  xlab("95% CrI width for weekly R (reported data)") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.77, 0.87),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 10, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3))
uncertainty1

###############################################
##  PANEL E: Difference in mean R (daily R)  ##
###############################################

scatter_dat2 <- data.frame(agg = r_agg$R$`Mean(R)`,
                           rep = r_rep$R$`Mean(R)`[6:32],
                           agg_upper = r_agg$R$`Quantile.0.975(R)`,
                           agg_lower = r_agg$R$`Quantile.0.025(R)`,
                           rep_upper = r_rep$R$`Quantile.0.975(R)`[6:32],
                           rep_lower = r_rep$R$`Quantile.0.025(R)`[6:32],
                           col = as.factor(inc_dat$interval[9:35]))


scatter_dat2$col <- dplyr::recode_factor(scatter_dat2$col,
                                         "early_daily" = "early")

corr_coeff <- cor(scatter_dat2$agg, scatter_dat2$rep)
r2 <- corr_coeff^2
round(r2, 2) # 0.13

corr2 <- ggplot(scatter_dat2, aes(x = rep, y = agg)) +
  geom_errorbar(aes(ymin = agg_lower, ymax = agg_upper, col = col),
                alpha = 1, size = 0.1) +
  geom_errorbarh(aes(xmin = rep_lower, xmax = rep_upper, col = col),
                 alpha = 1, size = 0.1) +
  geom_point(aes(col = col), size = 1.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  geom_hline(yintercept = 1, lty = 3) +
  geom_vline(xintercept = 1, lty = 3) +
  scale_colour_manual("",
                      breaks = c("early", "mid", "late"),
                      values = c("#B3E5FC", "dodgerblue", "#01579B")) +
  annotate("text", x = 3.1, y = 0.5, label  = "R^2==0.13", size = 3.5,
           parse = TRUE) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ylab("Mean daily R (reconstructed data)") +
  xlab("Mean daily R (reported data)") +
  ylim(0, 3.6) +
  xlim(0, 3.6) +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3))
corr2

######################################
##  PANEL F: Uncertainty (daily R)  ##
######################################

# Work out width of CrI's
scatter_dat2$width_weekly <- scatter_dat2$agg_upper - scatter_dat2$agg_lower
scatter_dat2$width_daily <- scatter_dat2$rep_upper - scatter_dat2$rep_lower

corr_coeff <- cor(scatter_dat2$width_weekly,
                  scatter_dat2$width_daily)
r2 <- corr_coeff^2
round(r2, 2) # 0.14

uncertainty2 <- ggplot(scatter_dat2, aes(x = width_daily, y = width_weekly)) +
  geom_point(aes(col = col), size = 1.5) +
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  scale_colour_manual("Time",
                      breaks = c("early", "mid", "late"),
                      values = c("#B3E5FC", "dodgerblue", "#01579B"),
                      guide = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  xlim(0, 0.18) +
  ylim(0, 0.18) +
  annotate("text", x = 0.15, y = 0.027, label  = "R^2==0.14", size = 3.5,
           parse = TRUE) +
  ylab("95% CrI width for daily R (reconstructed data)") +
  xlab("95% CrI width for daily R (reported data)") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.77, 0.87),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 10, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3))
uncertainty2

##################
## 6 PANEL PLOT ##
##################

plot_grid(plot_inc, squared_err, ncol = 1, align = "v")
p1 <- plot_inc + coord_cartesian(xlim = c(0, 35))
p2 <- squared_err + coord_cartesian(xlim = c(0, 35))

flu_plot <- plot_grid(p1, corr1, corr2, p2, uncertainty1, uncertainty2,
          ncol = 3, align = "hv",
          labels = c("A", "C", "E", "B", "D", "F"), label_size = 12)

ggsave(plot = flu_plot, "figures/Figure_2.pdf", device = cairo_pdf,
       width = 11, height = 7)
