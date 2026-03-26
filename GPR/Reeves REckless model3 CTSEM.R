setwd('/home/moosehunter/R/Fujita Lab/GPR/')
library(openesm)
library(rlang)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ctsem)

ls = list_datasets()
View(ls)
d1 = get_dataset('0034')
d2 = d1$data
View(d2)
head(d2)
# all id
d_all_id_clean <- d2 %>%
  dplyr::select(completion_date, reckless, id) %>%
  filter(!is.na(completion_date),
         !is.na(reckless)) %>%
  group_by(id) %>%
  mutate(
    time_day = as.numeric(
      difftime(
        completion_date,
        as.POSIXct(
          as.Date(min(completion_date)),
          tz = attr(completion_date, "tzone")
        ),
        units = "days"
      )
    )
  ) %>%
  ungroup()

# uncentered reckless for each id 
ggplot(d_all_id_clean,
       aes(x = time_day, y = reckless)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ id, scales = "free_x") +
  theme_minimal()


# C002 making the midnight of the 1st day 0, and making 1 day 1
d_C002 <- d2 %>%
  filter(id == 'C002') %>%
  dplyr::select(reckless, completion_date)%>%
  filter(!is.na(completion_date)) %>%
  mutate(
    time_day = as.numeric(
      difftime(
        completion_date,
        as.POSIXct(
          as.Date(min(completion_date)),
          tz = attr(completion_date, "tzone")
        ),
        units = "days"
      )
    )
  ) %>%
  mutate(centered_reckless = reckless - mean(reckless))


# uncentered C002 
ggplot(d_C002, aes(x = time_day, y = reckless)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "day",
    y = "Reckless",
    title = "Time Series of Reckless of C002"
  )

# centered C002 
ggplot(d_C002, aes(x = time_day, y = centered_reckless)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "day",
    y = "Reckless",
    title = "Time Series of Reckless of C002"
  )


# model 1  with train = 84, test = 21 20%
# prior specified by ChatGPT 
# matern 3/2 for trend, periodic kernel  p = 7 and Gaussian noise 
head(d_C002)
d_C002_train = d_C002 %>%
  slice(1:84)
d_C002_test = d_C002 %>%
  slice(85:104)

N_train <- 84

x_train <- d_C002$time_day[1:N_train]
y_train <- d_C002$centered_reckless[1:N_train]

x_new <- d_C002_test$time_day

d_ct <- d_C002 %>%
  dplyr::select(
    time = time_day,
    reckless = centered_reckless
  ) %>%
  mutate(id = 1) %>%        # ctsem needs an id column
  dplyr::arrange(time)

d_ct <- d_ct %>%
  dplyr::select(id, time, reckless)



ctmodel <- ctModel(
  type = "stanct",
  n.latent = 1,
  n.manifest = 1,
  manifestNames = "reckless",
  latentNames = "eta",
  LAMBDA = matrix(1),
  DRIFT = matrix("drift"),
  DIFFUSION = matrix("diffusion"),
  MANIFESTVAR = matrix("measerr")
)

fit <- ctStanFit(
  datalong = d_ct,
  ctstanmodel = ctmodel,
  iter = 2000,
  chains = 4
)
pred <- ctsem::ctStanGenerate(
  fit,
  nsamples = 100,
  burnin = 0,
  fullposterior = TRUE
)

# Extract simulated reckless values
sim_y <- pred$y
sim_time <- pred$t
