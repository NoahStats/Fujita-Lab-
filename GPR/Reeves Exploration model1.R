setwd('/home/moosehunter/R/Fujita Lab/GPR/')
library(openesm)
library(rlang)
library(dplyr)
library(ggplot2)
library(tidyr)

ls = list_datasets()
View(ls)
d1 = get_dataset('0034')
d2 = d1$data
View(d2)
head(d2)

# 1. Prepare the data for id 'C069'
d_c069_long <- d2 %>%
  # Filter for the specific individual and valid dates
  filter(id == 'C069', !is.na(completion_date)) %>%
  # Calculate time in days relative to their first entry
  mutate(
    time_day = as.numeric(
      difftime(
        completion_date,
        as.POSIXct(as.Date(min(completion_date)), tz = attr(completion_date, "tzone")),
        units = "days"
      )
    )
  ) %>%
  # Select the time variable and all numeric symptom variables
  # (Adjusting based on your column names; assuming columns from 'felt_enthusiastic' to the end are symptoms)
  select(time_day, felt_enthusiastic:distant, reckless) %>% 
  # Reshape to long format for plotting
  pivot_longer(
    cols = -time_day, 
    names_to = "variable", 
    values_to = "value"
  ) %>%
  # Remove NAs in the values so the lines connect properly
  filter(!is.na(value))

# 2. Plot all variables
ggplot(d_c069_long, aes(x = time_day, y = value)) +
  geom_line(color = "steelblue") +
  geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) + # Scales free_y as ranges vary
  theme_minimal() +
  labs(
    title = "Symptom Profiles for ID: C069",
    x = "Days since start",
    y = "Score"
  ) +
  theme(strip.text = element_text(size = 7)) # Shrink text to fit many facets
