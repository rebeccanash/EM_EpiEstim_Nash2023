## UK COVID by date of death: Supplementary analysis

library(hrbrthemes)
library(scales)
library(cowplot)

# data and results produced in main_analysis/flu.R
weekly_inc <- readRDS("supplementary_analysis/covid_deaths/weekly_inc_covid_deaths.rds")
r_rep <- readRDS("supplementary_analysis/covid_deaths/r_rep_result.rds")
r_agg <- readRDS("supplementary_analysis/covid_deaths/r_agg_result.rds")
r_rep_sliding <- readRDS("supplementary_analysis/covid_deaths/r_rep_sliding_result.rds")
r_agg_sliding <- readRDS("supplementary_analysis/covid_deaths/r_agg_sliding_result.rds")

daily_inc <- readRDS("real_data/daily_covid_deaths.rds")
recon_inc <- r_agg$I
index <- seq_along(daily_inc$Incidence)
inc_dat <- data.frame(index, daily_inc, recon_inc)

# non-sliding R
rep_dat <- data.frame(R = r_rep$R$`Mean(R)`[3:nrow(r_rep$R)],
                      upper = r_rep$R$`Quantile.0.975(R)`[3:nrow(r_rep$R)],
                      lower = r_rep$R$`Quantile.0.025(R)`[3:nrow(r_rep$R)],
                      day = r_rep$R$t_end[3:nrow(r_rep$R)])
agg_dat <- data.frame(R = r_agg$R$`Mean(R)`,
                      upper = r_agg$R$`Quantile.0.975(R)`,
                      lower = r_agg$R$`Quantile.0.025(R)`,
                      day = r_agg$R$t_end)
dat1 <- combine(rep_dat, agg_dat)

# sliding R
adj_t <- nrow(r_rep_sliding$R)
rep_dat_s <- data.frame(R = r_rep_sliding$R$`Mean(R)`[7:adj_t],
                        upper = r_rep_sliding$R$`Quantile.0.975(R)`[7:adj_t],
                        lower = r_rep_sliding$R$`Quantile.0.025(R)`[7:adj_t],
                        day = r_rep_sliding$R$t_end[7:adj_t])
agg_dat_s <- data.frame(R = r_agg_sliding$R$`Mean(R)`,
                        upper = r_agg_sliding$R$`Quantile.0.975(R)`,
                        lower = r_agg_sliding$R$`Quantile.0.025(R)`,
                        day = r_agg_sliding$R$t_end)

dat2 <- combine(rep_dat_s, agg_dat_s)


x_dates <- seq(as.Date('2020-03-02'), as.Date('2022-01-02'), by = "1 day")
date_1 <- format(x_dates[1], format = "%d %b %y")
date_2 <- format(x_dates[200], format = "%d %b %y")
date_3 <- format(x_dates[400], format = "%d %b %y")
date_4 <- format(x_dates[600], format = "%d %b %y")

####################
# R estimate plots #
####################

r_plot1 <- ggplot(dat1, aes(x = day, color = source,
                            fill = source)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  geom_hline(yintercept = 1, lty = 3) +
  geom_line(aes(y = R)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, col = NA) +
  scale_color_manual("",
                     values = c("rep_dat" = "darkgrey",
                                "agg_dat" = "darkgreen"),
                     labels = c("Reported data",
                                "Reconstructed data")) +
  scale_fill_manual("",
                    values = c("rep_dat" = alpha("darkgrey", 0.7),
                               "agg_dat" = alpha("seagreen4", 0.5)),
                    labels = c("Reported data",
                               "Reconstructed data")) +
  ylab("Daily R estimate") +
  xlab("") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.85, 0.75),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.ticks.x = element_line(linewidth = 0.1))
r_plot1

r_plot2 <- ggplot(dat2, aes(x = day, color = source,
                            fill = source)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  geom_hline(yintercept = 1, lty = 3) +
  geom_line(aes(y = R)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, col = NA) +
  scale_color_manual("",
                     values = c("rep_dat_s" = "darkgrey",
                                "agg_dat_s" = "darkgreen"),
                     labels = c("Reported data",
                                "Reconstructed data")) +
  scale_fill_manual("",
                    values = c("rep_dat_s" = alpha("darkgrey", 0.7),
                               "agg_dat_s" = alpha("seagreen4", 0.5)),
                    labels = c("Reported data",
                               "Reconstructed data")) +
  ylab("Weekly sliding R estimate") +
  xlab("Date") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.ticks.x = element_line(linewidth = 0.1))
r_plot2

p1 <- r_plot1 + coord_cartesian(xlim = c(5, 679), ylim = c(0, 4))
p2 <- r_plot2 + coord_cartesian(xlim = c(5, 679), ylim = c(0, 4))
r_plots_covid <- plot_grid(p1, p2, ncol = 1, align = "v",
                           labels = c("A", "B"), label_size = 12)

ggsave(plot = r_plots_covid, "supplementary_figures/covid_deaths_rt.pdf",
       device = cairo_pdf, width = 9, height = 7.5)

##############
# Bias plots #
##############

## Incidence

plot_inc <- ggplot(inc_dat, aes(x = index)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  scale_y_continuous(labels = comma) +
  geom_line(aes(y = Incidence, colour = "Reported")) +
  geom_line(aes(y = recon_inc, colour = "Reconstructed")) +
  scale_colour_manual("",
                      breaks = c("Reported", "Reconstructed"),
                      values = c("Reported" = alpha("darkgrey", 0.8),
                                 "Reconstructed" = "#297a4d")) +
  ylab("Daily incidence \n by date of death") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.73, 0.77),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.05)) +
  guides(colour = guide_legend(override.aes = list(alpha = c(0.8, 1))))
plot_inc

plot_log_inc <- ggplot(inc_dat, aes(x = index)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  scale_y_log10(breaks = scales::log_breaks(), labels = comma) +
  geom_line(aes(y = Incidence, colour = "Reported"), alpha = 0.8) +
  geom_line(aes(y = recon_inc, colour = "Reconstructed")) +
  scale_colour_manual("",
                      values = c("Reported" = alpha("darkgrey", 0.8),
                                 "Reconstructed" = "#297a4d")) +
  ylab("Daily incidence by \n date of death (log scale)") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "lightgrey", linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.05)) +
  guides(colour = guide_legend(override.aes = list(alpha = c(0.8, 1))))
plot_log_inc

## Absolute difference

diff_daily <- rep_dat$R - agg_dat$R
day_daily <- rep_dat$day 
plot_diff_daily <- data.frame(diff_daily, day_daily)

diff_plot1 <- ggplot(plot_diff_daily,
                     aes(x = day_daily, y = diff_daily)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Daily R (reported data) \n - Daily R (reconstructed data)") +
  xlab("") +
  hrbrthemes::theme_ipsum() +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),                       
    panel.grid.major.y = element_blank(),                        
    panel.grid.major.x = element_line(linetype = 2),                   
    axis.line = element_line(color="lightgrey", linewidth = 0.3))
diff_plot1

diff_weekly <- rep_dat_s$R - agg_dat_s$R
day_weekly <- rep_dat_s$day 
plot_diff_weekly <- data.frame(diff_weekly, diff_weekly)

diff_plot2 <- ggplot(plot_diff_weekly,
                     aes(x = day_weekly, y = diff_weekly)) +
  scale_x_continuous(expand = c(0.05, 0.01), breaks = c(1, 200, 400, 600),
                     labels = c(date_1, date_2, date_3, date_4)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Weekly R (reported data) \n - Weekly R (reconstructed data)") +
  xlab("Date") +
  hrbrthemes::theme_ipsum() +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.tag.position = c(0, 1.05),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                                size = 12, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),
                                size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),                       
    panel.grid.major.y = element_blank(),                        
    panel.grid.major.x = element_line(linetype = 2),                   
    axis.line = element_line(color="lightgrey", linewidth = 0.3))
diff_plot2

## Plot together

panel_a <- plot_inc + coord_cartesian(xlim = c(1, 672))
panel_b <- plot_log_inc + coord_cartesian(xlim = c(1, 672))
panel_c <- diff_plot1 + coord_cartesian(xlim = c(1, 672))
panel_d <- diff_plot2 + coord_cartesian(xlim = c(1, 672))
cov_bias_d <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 1,
                        align = "v", labels=c("A", "B", "C", "D"), 
                        label_size = 12)

ggsave(plot = cov_bias_d, "supplementary_figures/covid_deaths_rt_bias.pdf",
       device = cairo_pdf, width = 7.5, height = 9)
