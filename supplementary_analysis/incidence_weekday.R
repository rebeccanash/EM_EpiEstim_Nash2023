## Relative difference in reporting by weekday

# Packages required
packages <- c("dplyr", "ggplot2", "hrbrthemes", "cowplot")

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

flu <- readRDS("real_data/daily_flu.rds")
cov_c <- readRDS("real_data/daily_covid_cases.rds")
cov_d <- readRDS("real_data/daily_covid_deaths.rds")

weekdays <- c("Monday", "Tuesday", "Wednesday", "Thursday",
              "Friday", "Saturday", "Sunday")

# flu
flu$day <- factor(weekdays(as.Date(flu$Date, format="%d/%m/%Y")),
                  levels = weekdays)
flu$week <- as.numeric(format(as.Date(flu$Date, format="%d/%m/%Y"), "%V"))
comp_weeks <- aggregate(data.frame(count = flu$week),
                        list(nday = flu$week), length)
rm_weeks <- comp_weeks$nday[comp_weeks$count != 7]
flu <- flu[!flu$week %in% rm_weeks, ]
flu$weekly_total <- with(flu, ave(Incidence, week, FUN = sum))
flu$weekly_mean <- with(flu, ave(Incidence, week, FUN = mean))

for (row in 1:nrow(flu)) { 
  flu$rel_diff[row] <-
    (flu$Incidence[row] - flu$weekly_mean[row]) / flu$weekly_mean[row]
}

# covid cases
cov_c$day <- factor(weekdays(as.Date(cov_c$Date, format="%d/%m/%Y")),
                    levels = weekdays)
cov_c$week <- as.numeric(format(as.Date(cov_c$Date, format="%d/%m/%Y"), "%V"))
cov_c$week[319:679] <- cov_c$week[319:679] + 53
comp_weeks <- aggregate(data.frame(count = cov_c$week),
                        list(nday = cov_c$week), length)
rm_weeks <- c(comp_weeks$nday[comp_weeks$count != 7], 9)
cov_c <- cov_c[!cov_c$week %in% rm_weeks, ]
cov_c$weekly_total <- with(cov_c, ave(Incidence, week, FUN = sum))
cov_c$weekly_mean <- with(cov_c, ave(Incidence, week, FUN = mean))

for (row in 1:nrow(cov_c)) { 
  cov_c$rel_diff[row] <-
    (cov_c$Incidence[row] - cov_c$weekly_mean[row]) / cov_c$weekly_mean[row]
}

# covid deaths
cov_d$day <- factor(weekdays(as.Date(cov_d$Date, format="%d/%m/%Y")),
                    levels = weekdays)
cov_d$week <- as.numeric(format(as.Date(cov_d$Date, format="%d/%m/%Y"), "%V"))
cov_d$week[372:672] <- cov_d$week[372:672] + 53
comp_weeks <- aggregate(data.frame(count = cov_d$week),
                        list(nday = cov_d$week), length)
rm_weeks <- comp_weeks$nday[comp_weeks$count != 7]
cov_d <- cov_d[!cov_d$week %in% rm_weeks, ]
cov_d$weekly_total <- with(cov_d, ave(Incidence, week, FUN = sum))
cov_d$weekly_mean <- with(cov_d, ave(Incidence, week, FUN = mean))

for (row in 1:nrow(cov_d)) { 
  cov_d$rel_diff[row] <-
    (cov_d$Incidence[row] - cov_d$weekly_mean[row]) / cov_d$weekly_mean[row]
}

# plot
weekday_plot <- function(dat, x_axis_title = element_blank()) {
  ggplot(dat, aes(x = day, y = rel_diff)) + 
    geom_hline(yintercept = 0, linetype = 2, colour = "grey")+
    theme_ipsum() +
    xlab("Weekday") +
    ylab("Relative difference from weekly mean") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.75) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.tag.position = c(0, 1.05),
      axis.title.x = x_axis_title,
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                  size=12, hjust=0.5),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),                    
      axis.line = element_line(color="lightgrey", linewidth = 0.3),
      axis.ticks.x = element_blank()) + 
    geom_jitter(aes(fill = day), shape = 21,
                position = position_jitter(0.2, seed = 5)) 
} 

plot1 <- weekday_plot(flu)
plot2 <- weekday_plot(cov_c)
plot3 <- weekday_plot(cov_d, x_axis_title = element_text(
  margin = margin(t = 20, r = 0, b = 0, l = 0), size=12, hjust=0.5))

# plot together
weekday_3panel <- plot_grid(plot1, plot2, plot3, ncol = 1, align = "v",
                            labels=c("A", "B", "C"), label_size = 12)

# save
ggsave(plot = weekday_3panel,"supplementary_figures/incidence_weekday.pdf",
       device = cairo_pdf, width = 9,height = 8.5)
