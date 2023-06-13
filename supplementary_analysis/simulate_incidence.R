## Simulate incidence

library(projections)
library(distcrete)

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

###############
# Constant Rt #
###############

proj <- list()
rt <- c(1, 1.25, 1.5, 1.75)

for (i in seq_along(rt)){
set.seed(6)
proj[[i]] <- project(
  x = initial_i,
  R = rt[i],
  si = si_one,
  n_days = n_days,
  n_sim = n_sim,
  instantaneous_R = TRUE
)
}

proj_data <- list()
for (i in seq_along(rt)){
proj_data[[i]] <- as.data.frame(proj[[i]], long = TRUE)
}

name_rt <- as.character(rt)
for(x in seq_along(rt)) {
  saveRDS(proj_data[[x]],
          file = paste0("supplementary_analysis/simulated_incidence/sims_rt_",
                        name_rt[x], ".rds"))
}

###########################
# Time-varying Rt: Sudden #
###########################

proj_step <- list()
rt <- list()

# Increasing
rt[[1]] <- c(1, 1.25)
rt[[2]] <- c(1.25, 1.5)
rt[[3]] <- c(1.5, 1.75)

# Decreasing
rt[[4]] <- c(1.25, 0.75)
rt[[5]] <- c(1.25, 1)
rt[[6]] <- c(1.5, 1.25)
rt[[7]] <- c(1.75, 1.5)

# Time change day 34 (last day you want old R minus 1)
# days 1 to 35 = step_from, days 36 to 70 = step_to

for (i in seq_along(rt)) {
set.seed(5)
proj_step[[i]] <- project(
  x = initial_i,
  R = rt[[i]],
  si = si_one,
  n_days = n_days,
  n_sim = n_sim,
  time_change = 34,
  instantaneous_R = TRUE
)
}

proj_step_data <- list()
for (i in seq_along(rt)){
  proj_step_data[[i]] <- as.data.frame(proj_step[[i]], long = TRUE)
}


name_rt <- c("1to1.25", "1.25to1.5", "1.5to1.75",
             "1.25to0.75", "1.25to1", "1.5to1.25", "1.75to1.5")

for(x in seq_along(rt)) {
  saveRDS(proj_step_data[[x]], file =
            paste0("supplementary_analysis/simulated_incidence/sims_rt_sud_",
                   name_rt[x], ".rds"))
}


############################
# Time-varying Rt: Gradual #
############################

proj_grad <- list()
rt <- list()

# Increasing
rt[[1]] <- seq(from = 1, to = 1.25, length = 32)
rt[[2]] <- seq(from = 1.25, to = 1.5, length = 32)
rt[[3]] <- seq(from = 1.5, to = 1.75, length = 32)

# Decreasing
rt[[4]] <- seq(from = 1.25, to = 0.75, length = 32)
rt[[5]] <- seq(from = 1.25, to = 1, length = 32)
rt[[6]] <- seq(from = 1.5, to = 1.25, length = 32)
rt[[7]] <- seq(from = 1.75, to = 1.5, length = 32)


# time change last day you want old R minus 1 
# step_from = days 1-20, gradual change = days 21-49, step_to = days 50-70

for (i in seq_along(rt)){
  set.seed(5)
  proj_grad[[i]] <- project(
    x = initial_i,
    R = rt[[i]],
    si = si_one,
    n_days = n_days,
    n_sim = n_sim,
    time_change = c(seq(19,49)),
    instantaneous_R = TRUE)
}

proj_grad_data <- list()
for (i in seq_along(rt)){
  proj_grad_data[[i]] <- as.data.frame(proj_grad[[i]], long = TRUE)
}

name_rt <- c("1to1.25", "1.25to1.5", "1.5to1.75",
             "1.25to0.75", "1.25to1", "1.5to1.25", "1.75to1.5")

for(x in seq_along(rt)) {
  saveRDS(proj_grad_data[[x]], file =
            paste0("supplementary_analysis/simulated_incidence/sims_rt_grad_",
                   name_rt[x], ".rds"))
}
