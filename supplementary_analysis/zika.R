## Zika virus disease - Compare typical workaround method to new method

library(EpiEstim)
library(tidyr)
library(hrbrthemes)

# Taken data for zika in Rio de Janeiro from the supplements of:
# https://www.science.org/doi/10.1126/science.aag0219

# load data (22 weeks)
dat <- read.csv("supplementary_analysis/zika/Brazil_RioJaneiro.csv")

# Serial interval taken from Ferguson 2016
# Workaround involved supplying SI at same time scale
method <- "parametric_si"
old_mean_si <- 20.0 / 7
old_sd_si <- 7.4 / 7

nweeks <- nrow(dat)
t_start <- seq(from = 2, to = nweeks - 1, 1)
old_config <- EpiEstim::make_config(mean_si = old_mean_si,
                                    std_si = old_sd_si,
                                    t_start = t_start,
                                    t_end = t_start + 1)

old_result <- estimate_R(incid = dat$incidence,
                         method = method,
                         config = old_config)

# EM method allows you to supply SI as normal (on a daily timescale)
mean_si <- 20.0
sd_si <- 7.4
config <- EpiEstim::make_config(mean_si = mean_si,
                                std_si = sd_si)
  
new_result <- estimate_R(incid = dat$incidence,
                         dt = 7L,
                         dt_out = 7L,
                         recon_opt = "match",
                         method = method,
                         config = config)
ndays <- length(t_start) * 7
overall_time <- 1:(nweeks * 7)

old_r <- data.frame(R = c(rep(NA, 14), rep(old_result$R$`Mean(R)`, each = 7)),
           upper = c(rep(NA, 14), rep(old_result$R$`Quantile.0.975(R)`, each = 7)),
           lower = c(rep(NA, 14), rep(old_result$R$`Quantile.0.025(R)`, each = 7)),
           time = overall_time,
           source = "old")

new_r <- data.frame(R = c(rep(NA, 20), new_result$R$`Mean(R)`),
                    upper = c(rep(NA, 20), new_result$R$`Quantile.0.975(R)`),
                    lower = c(rep(NA, 20), new_result$R$`Quantile.0.025(R)`),
                    time = overall_time,
                    source = "new")

plot_dat <- combine(old_r, new_r)
plot_dat <- plot_dat[plot_dat$time > 35,]

x_dates <- seq(as.Date("2015-11-14"), as.Date("2016-04-15"), by = "1 day")
date_1 <- format(x_dates[60], format = "%d %b %y")
date_2 <- format(x_dates[90], format = "%d %b %y")
date_3 <- format(x_dates[120], format = "%d %b %y")
date_4 <- format(x_dates[150], format = "%d %b %y")

# plot
zika_plot <- ggplot(plot_dat, aes(x = time, color = source, fill = source)) +
  scale_x_continuous(breaks = c(60, 90, 120, 150),
                     labels = c(date_1, date_2, date_3, date_4)) +
  geom_hline(yintercept = 1, lty = 3) +
  geom_line(aes(y = R)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, col = NA) +
  scale_color_manual("",
                     values = c("old" = "grey",
                                "new" = "darkgreen"),
                     labels = c("old" = "Workaround method",
                                "new" = "EM algorithm")) +
  scale_fill_manual("",
                    values = c("old" = alpha("grey", 0.6),
                               "new" = alpha("seagreen4", 0.5)),
                    labels = c("old" = "Workaround method",
                               "new" = "EM algorithm")) +
  ylab("R estimate") +
  xlab("Time") +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position = c(0.8, 0.85),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 0, r = 20, b = 20, l = 20),
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
    axis.ticks.x = element_line(linewidth = 0.1)) +
  coord_cartesian(ylim = c(0,4), clip = "off")

# save
ggsave(plot = zika_plot, "supplementary_figures/zika_rt.pdf",
       device = cairo_pdf, width = 9, height = 4)
