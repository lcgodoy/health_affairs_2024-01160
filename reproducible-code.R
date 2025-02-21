library(rstan)
library(dplyr)
library(ggplot2)
library(lubridate)

rstan_options(auto_write = TRUE)

## the dataset is not publicly available
my_dt <- readr::read_csv("data/processed.csv",
                         show_col_types = FALSE)

## reference coding
ref_county <- "NEW HAVEN"

my_dt <- my_dt |>
  mutate(c1 = ifelse(county == "FAIRFIELD", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c2 = ifelse(county == "HARTFORD", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c3 = ifelse(county == "LITCHFIELD", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c4 = ifelse(county == "MIDDLESEX", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c5 = -1,
         c6 = ifelse(county == "NEW LONDON", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c7 = ifelse(county == "TOLLAND", 1,
              ifelse(county == ref_county,
                     -1, 0)),
         c8 = ifelse(county == "WINDHAM", 1,
              ifelse(county == ref_county,
                     -1, 0)))

my_dt <- my_dt |>
  mutate(id_county = as.numeric(factor(county, levels = unique(county))))

my_dt <- my_dt |>
  mutate(time_date = ym(ym)) |>
  mutate(time = interval(min(time_date), time_date) %/%
           months(1) +
           1L, .before = ym)

my_dt <- my_dt |>
  mutate(fips = case_when(
             county == "FAIRFIELD" ~ "09001",
             county == "HARTFORD" ~ "09003",
             county == "LITCHFIELD" ~ "09005",
             county == "MIDDLESEX" ~ "09007",
             county == "NEW HAVEN" ~ "09009",
             county == "NEW LONDON" ~ "09011",
             county == "TOLLAND" ~ "09013",
             county == "WINDHAM" ~ "09015"
         ))

sd_deaths <- my_dt |>
  filter(county == "NEW HAVEN") |>
  pull(covid_deaths) |>
  sd(na.rm = TRUE)

my_dt <- my_dt |>
  mutate(covid_deaths = ifelse(is.na(covid),
                               - sd_deaths,
                               as.numeric(covid_deaths)))

## using quarters instead of months
my_dt <- my_dt |>
  mutate(time10 = time / 3,
         .after = time)

time_center <- my_dt |>
  filter(time_date == "2019-04-01") |>
  pull(time) |>
  unique()

## centering time (for interpretation of intercept and capta)
my_dt <- my_dt |>
  mutate(time1 = (time - time_center) / 3)

## variables to be centered
cont_covars <- c("nas",
                 "black_b",
                 "late_nopnc",
                 "wic_part",
                 "no_sm_pnc",
                 "m_married")

my_dt <- my_dt |>
  mutate(across(all_of(cont_covars),
                ~ 100 * as.numeric(scale(., scale = FALSE)),
                .names = "sc_{.col}"))

my_dt <- my_dt |>
  mutate(capta = if_else(time_date >= "2019-04-01", 1, 0))

my_dt <- my_dt |>
  mutate(capta_f = factor(capta, levels = c(0, 1),
                          labels = c("Before", "After")))

my_dt <- my_dt |>
  mutate(sc_covid_deaths = as.numeric(scale(covid_deaths)))

##--- Model SEI ----

poisson_reg <- stan_model(file = "stan/final-ar.stan",
                          allow_optimizations = TRUE)

## design matrix
make_x <-
  model.matrix(~ 0 +
                 capta * time1 +
                 sc_nas +
                 sc_black_b +
                 sc_late_nopnc +
                 sc_wic_part +
                 sc_no_sm_pnc +
                 sc_m_married +
                 ## covid_f +
                 sc_covid_deaths +
                 c1 + c2 + c3 + c4 +
                 c6 + c7 + c8,
               data = my_dt)

## make_x <- make_x[, - which(colnames(make_x) == "covid_fBefore")]

##--- correlations between predictors ----

## make_x[, - c(1, seq(NCOL(make_x) - 11,
cor_mat <-
  make_x[, - c(1, seq(NCOL(make_x) - 7,
                    NCOL(make_x)))] |>
  cor()

cor_mat

summary(abs(cor_mat[upper.tri(cor_mat)]))

my_dt <-
  mutate(my_dt, time_f = as.factor(time))

## data for marginal predictions
new_data <-
  data.frame(sc_nas = 0,
             sc_black_b = 0,
             sc_late_nopnc = 0,
             sc_wic_part = 0,
             sc_no_sm_pnc = 0,
             sc_m_married = 0,
             sc_covid_deaths = 0,
             ## f_covid = factor("[-0.00492,1.64)",
             ##                  levels = levels(my_dt$covid_f)),
             c1 = 0, c2 = 0, c3 = 0,
             c4 = 0, ## c5 = 0,
             c6 = 0,
             c7 = 0, c8 = 0,
             capta = 0,
             county = "HARTFORD",
             births = 10000,
             ## time10 =
             ##   c(seq(from = min(my_dt$time10),
             ##         to   = max(my_dt$time10),
             ##         length.out =
             ##           floor(length(unique(my_dt$time10)) / 3)),
             ##     max(my_dt$time10)))
             time1 = unique(my_dt$time1),
             time = unique(my_dt$time))

new_data <-
  new_data |>
  mutate(capta = ifelse(time >= 25, 1, capta)) |>
  mutate(time_capta = time1 * capta)

## design matrix for predictions
new_x <-
  model.matrix(~ 0 +
                 capta * time1 +
                 sc_nas +
                 sc_black_b +
                 sc_late_nopnc +
                 sc_wic_part +
                 sc_no_sm_pnc +
                 sc_m_married +
                 sc_covid_deaths +
                 ## f_covid +
                 c1 + c2 + c3 + c4 + c6 + c7 + c8,
               data = new_data)

## new_x <- new_x[, - which(colnames(new_x) == "f_covidBefore")]

new_data <-
  mutate(new_data, time_f = as.factor(time))

## data to feed stan
model_data <-
  list(
      N = nrow(my_dt),
      y = my_dt$sei_report,
      logoff = log(my_dt$births),
      K = ncol(make_x),
      toggle_reg = 1,
      X = make_x,
      toggle_un = 0,
      toggle_qr = 1,
      toggle_pred = 1,
      N_pred = array(NROW(new_x), dim = 1),
      lo_pred = array(log(10000), dim = 1),
      X_pred = new_x,
      time_pred = new_data$time,
      N_i = 0,
      group = integer(),
      toggle_ar = 1,
      N_t = max(my_dt$time),
      time = my_dt$time,
      p_eta = .1,
      eta_0 = .9
  )

sei_fit <- sampling(poisson_reg,
                    data = model_data,
                    seed = 2024,
                    cores = 4,
                    chains = 4)

to_from <-
  data.frame(parameter =
               c("alpha", sprintf("beta[%d]", seq_len(model_data$K)),
                 "beta_ct"),
             coefficient =
               c("Intercept", colnames(make_x),
                 "rel_increase_after"))

##--- viz simpler model ----

last_beta <-
  sprintf("beta[%d]", model_data$K)

out_spt <-
  sei_fit |>
  tidybayes::spread_draws(alpha,
                          `beta\\[.*`,
                          regex = TRUE) |>
  select(- .chain, - .iteration, - .draw) |>
  mutate(beta_ct = `beta[17]` + `beta[2]`) |>
  mutate(across(1:(model_data$K + 2),
                .fns = exp)) |>
  tidyr::pivot_longer(cols = 1:(model_data$K + 2),
                      names_to = "variable") |>
  group_by(variable) |>
  summarise(med = median(value),
            q1 = quantile(value, probs = .1),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            q9 = quantile(value, probs = .95)) |>
  ungroup() |>
  mutate(median = med,
         l_ci80 = q1,
         u_ci80 = q9,
         l_ci95 = q025,
         u_ci95 = q975) |>
  select(variable, median, l_ci80, u_ci80,
         l_ci95, u_ci95)

out_spt <-
  left_join(out_spt, to_from,
            by = c("variable" = "parameter")) |>
  relocate(coefficient, .before = variable) |>
  select(- variable)

## sei model parameter estimates
out_spt |>
  print(n = Inf)

out_hyper <-
  sei_fit |>
  tidybayes::spread_draws(`eta[1]`, `sigma_t[1]`) |>
  select(- .chain, - .iteration, - .draw) |>
  tidyr::pivot_longer(cols = `eta[1]`:`sigma_t[1]`,
                      names_to = "variable") |>
  group_by(variable) |>
  summarise(med = median(value),
            q1 = quantile(value, probs = .1),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            q9 = quantile(value, probs = .95)) |>
  ungroup() |>
  mutate(median = med,
         l_ci80 = q1,
         u_ci80 = q9,
         l_ci95 = q025,
         u_ci95 = q975) |>
  select(variable, median, l_ci80, u_ci80,
         l_ci95, u_ci95)


## sei model hyperparameters estimates
out_hyper

##--- predictions ----

preds <- extract(sei_fit, pars = "mu_pred",
                 permuted = TRUE)[[1]]

preds <- tibble(median = apply(preds, 2, median),
                l_ci80 = apply(preds, 2, quantile, probs = .1),
                u_ci80 = apply(preds, 2, quantile, probs = .9),
                l_ci95 = apply(preds, 2, quantile, probs = .025),
                u_ci95 = apply(preds, 2, quantile, probs = .975))

preds <- bind_cols(new_data, preds)

preds <- preds |>
  mutate(time10 = time / 3,
         .after = time)

## marginal predictions (all the covariates held fixed at reference level)

preds <-
  left_join(preds,
            distinct(select(my_dt, time_date, time1)),
            by = "time1")

ggplot(data = preds,
       aes(x = time_date,
           y = median)) +
  geom_line(aes(color = as.factor(capta))) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci95,
                     ymax = u_ci95)) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci80,
                     ymax = u_ci80),
                 lwd = 1.2) +
  geom_point(aes(fill = as.factor(capta)),
             pch = 21,
             color = 1, size = 3) +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Reports per 10k births",
       x = NULL,
       color = NULL)

readr::write_csv(select(preds, time_date, median, l_ci80:u_ci95),
                 file = "data/data-exhibit-4.csv")

## EXHIBIT 4
ggsave(filename = "figs/exhibit-4.pdf",
       width = 12,
       height = 7,
       dpi = 300)

##--- additional table (supplement) ----


capta_timepoints <-
  sei_fit |>
  tidybayes::spread_draws(`beta[1]`, `beta[2]`, `beta[17]`) |>
  select(- .chain, - .iteration, - .draw) |>
  mutate(capta_12 = `beta[1]` + 12 * (`beta[2]` + `beta[17]`)) |>
  mutate(capta_24 = `beta[1]` + 16 * (`beta[2]` + `beta[17]`)) |>
  mutate(capta_36 = `beta[1]` + 20 * (`beta[2]` + `beta[17]`)) |>
  select(starts_with("capta")) |>
  tidyr::pivot_longer(1:3) |>
  mutate(value = exp(value)) |>
  group_by(name) |>
  summarise(median = median(value),
            l_ci80 = quantile(value, probs = .1),
            u_ci80 = quantile(value, probs = .9),
            l_ci95 = quantile(value,  probs = .025),
            u_ci95 = quantile(value,  probs = .975)) |>
  transmute(coefficient = name,
            median = round(median, 2),
            ci_80 = sprintf("(%.2f - %.2f)", l_ci80, u_ci80),
            ci_95 = sprintf("(%.2f - %.2f)", l_ci95, u_ci95))

preds_timepoints <-
  preds |>
  filter(time_date %in% c("2019-02-01", sprintf("%s-03-01",
                                                c(2020, 2021, 2022)))) |>
      select(time_date, capta, median, l_ci80:u_ci95) |>
      transmute(coefficient = sprintf("Predictions at:  %s",
                                      as.character(time_date)),
                median = round(median / 100, 2),
                ci_80 = sprintf("(%.2f - %.2f)", l_ci80 / 100, u_ci80 / 100),
                ci_95 = sprintf("(%.2f - %.2f)", l_ci95 / 100, u_ci95 / 100))

## SUPPLEMENTAL EXHIBIT 8
bind_rows(
    preds_timepoints,
    capta_timepoints
)

##--- model fostercare ----

my_dt2 <- my_dt |>
  filter(my_dt$sei_report > 0)

make_x2 <-
  model.matrix(~ 0 +
                 capta * time1 +
                 sc_nas +
                 sc_black_b +
                 sc_late_nopnc +
                 sc_wic_part +
                 sc_no_sm_pnc +
                 sc_m_married +
                 sc_covid_deaths +
                 c1 + c2 + c3 + c4 + c6 + c7 + c8,
               data = my_dt2)

## data for predictions
new_data2 <-
  data.frame(sc_nas = 0,
             sc_black_b = 0,
             sc_late_nopnc = 0,
             sc_wic_part = 0,
             sc_no_sm_pnc = 0,
             sc_m_married = 0,
             sc_covid_deaths = 0,
             c1 = 0, c2 = 0, c3 = 0,
             c4 = 0, c6 = 0,
             c7 = 0, c8 = 0,
             capta = 0,
             county = "NEW HAVEN",
             time1 = unique(my_dt2$time1),
             time = unique(my_dt2$time))

new_data2 <-
  new_data2 |>
  mutate(capta = ifelse(time >= 25, 1, capta)) |>
  mutate(time_capta = time1 * capta)

new_x2 <-
  model.matrix(~ 0 +
                 capta * time1 +
                 sc_nas +
                 sc_black_b +
                 sc_late_nopnc +
                 sc_wic_part +
                 sc_no_sm_pnc +
                 sc_m_married +
                 sc_covid_deaths +
                 c1 + c2 + c3 + c4 + c6 + c7 + c8,
               data = new_data2)

model_data2 <-
  list(
      N = nrow(my_dt2),
      y = my_dt2$fostercare,
      logoff = log(my_dt2$sei_report),
      K = ncol(make_x2),
      toggle_reg = 1,
      X = make_x2,
      toggle_un = 0,
      toggle_qr = 1,
      toggle_pred = 1,
      N_pred = array(NROW(new_x2), dim = 1),
      lo_pred = array(log(100), dim = 1),
      X_pred = new_x2,
      time_pred = new_data2$time,
      N_i = 0,
      group = integer(),
      toggle_ar = 1,
      N_t = max(my_dt2$time),
      time = my_dt2$time,
      p_eta = .1,
      eta_0 = .9
  )

foster_fit <- sampling(poisson_reg,
                       data = model_data2,
                       seed = 2024,
                       cores = 4,
                       chains = 4)

to_from <-
  data.frame(parameter =
               c("alpha", sprintf("beta[%d]", seq_len(model_data2$K)), "beta_ct"),
             coefficient =
               c("Intercept", colnames(make_x2), "rel_increase_after"))


##--- viz simpler model ----

last_beta2 <-
  sprintf("beta[%d]", model_data2$K)

out_spt2 <-
  foster_fit |>
  tidybayes::spread_draws(alpha,
                          `beta\\[.*`,
                          regex = TRUE) |>
  select(- .chain, - .iteration, - .draw) |>
  mutate(beta_ct = `beta[17]` + `beta[2]`) |>
  mutate(across(1:(model_data$K + 2),
                .fns = exp)) |>
  tidyr::pivot_longer(cols = 1:(model_data$K + 2),
                      names_to = "variable") |>
  group_by(variable) |>
  summarise(med = median(value),
            q1 = quantile(value, probs = .1),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            q9 = quantile(value, probs = .95)) |>
  ungroup() |>
  mutate(median = med,
         l_ci80 = q1,
         u_ci80 = q9,
         l_ci95 = q025,
         u_ci95 = q975) |>
  select(variable, median, l_ci80, u_ci80,
         l_ci95, u_ci95)

out_spt2 <-
  left_join(out_spt2, to_from,
            by = c("variable" = "parameter")) |>
  relocate(coefficient, .before = variable) |>
  select(- variable)

## foster model parameter estimates
out_spt2 |>
  print(n = Inf)

out_hyper2 <-
  foster_fit |>
  tidybayes::spread_draws(`eta[1]`, `sigma_t[1]`) |>
  select(- .chain, - .iteration, - .draw) |>
  tidyr::pivot_longer(cols = `eta[1]`:`sigma_t[1]`,
                      names_to = "variable") |>
  group_by(variable) |>
  summarise(med = median(value),
            q1 = quantile(value, probs = .1),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            q9 = quantile(value, probs = .95)) |>
  ungroup() |>
  mutate(median = med,
         l_ci80 = q1,
         u_ci80 = q9,
         l_ci95 = q025,
         u_ci95 = q975) |>
  select(variable, median, l_ci80, u_ci80,
         l_ci95, u_ci95)

## foster model parameter estimates
out_hyper2 |>
  print(n = Inf)

## EXHIBIT 3
out_spt[c(10, 11, 9, 19), ]

## exhibit 5
out_spt2[c(10, 11, 9, 19), ]

##--- predictions ----

preds2 <- extract(foster_fit, pars = "mu_pred",
                  permuted = TRUE)[[1]]

preds2 <- tibble(median = apply(preds2, 2, median),
                 l_ci80 = apply(preds2, 2, quantile, probs = .1),
                 u_ci80 = apply(preds2, 2, quantile, probs = .9),
                 l_ci95 = apply(preds2, 2, quantile, probs = .025),
                 u_ci95 = apply(preds2, 2, quantile, probs = .975))

preds2 <- bind_cols(new_data2, preds2)

preds2 <-
  left_join(preds2,
            distinct(select(my_dt, time_date, time1)),
            by = "time1")

ggplot(data = preds2,
       aes(x = time_date,
           y = median)) +
  geom_line(aes(color = as.factor(capta))) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci95,
                     ymax = u_ci95)) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci80,
                     ymax = u_ci80),
                 lwd = 1.2) +
  geom_point(aes(fill = as.factor(capta)),
             pch = 21,
             color = 1, size = 3) +
  ## geom_vline(xintercept = as_date("2019-04-01"),
  ##            lty = 2, color = 2, lwd = 1.4) +
  theme_bw() +
  guides(color = FALSE,
         fill = FALSE) +
  labs(y = "Placements per 100 reports",
       x = NULL,
       color = NULL)

readr::write_csv(select(preds2, time_date, median, l_ci80:u_ci95),
                 file = "data/data-exhibit-6.csv")

## EXHIBIT 6
ggsave(filename = "figs/exhibit-6.pdf",
       width = 12,
       height = 7,
       dpi = 300)

##--- Negative control analysis ----

model_data3 <- model_data

## changing response variable
model_data3$y <- my_dt$non_sei_nb

nc_fit <- sampling(poisson_reg,
                   data = model_data3,
                   seed = 2024,
                   cores = 4,
                   chains = 4)

to_from <-
  data.frame(parameter =
               c("alpha", sprintf("beta[%d]", seq_len(model_data$K)),
                 "beta_ct"),
             coefficient =
               c("Intercept", colnames(make_x),
                 "rel_increase_after"))

nc_pars <-
  nc_fit |>
  tidybayes::spread_draws(alpha,
                          `beta\\[.*`,
                          regex = TRUE) |>
  select(- .chain, - .iteration, - .draw) |>
  mutate(beta_ct = `beta[17]` + `beta[2]`) |>
  mutate(across(1:(model_data$K + 2),
                .fns = exp)) |>
  tidyr::pivot_longer(cols = 1:(model_data$K + 2),
                      names_to = "variable") |>
  group_by(variable) |>
  summarise(med = median(value),
            q1 = quantile(value, probs = .1),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            q9 = quantile(value, probs = .95)) |>
  ungroup() |>
  mutate(median = med,
         l_ci80 = q1,
         u_ci80 = q9,
         l_ci95 = q025,
         u_ci95 = q975) |>
  select(variable, median, l_ci80, u_ci80,
         l_ci95, u_ci95)

nc_pars <-
  left_join(nc_pars, to_from,
            by = c("variable" = "parameter")) |>
  relocate(coefficient, .before = variable) |>
  select(- variable)

## Estimates negative control
nc_pars[c(10, 11, 9, 19), ]

##--- graph for negative controls ----

preds_nc <- extract(nc_fit, pars = "mu_pred",
                    permuted = TRUE)[[1]]

preds_nc <- tibble(median = apply(preds_nc, 2, median),
                   l_ci80 = apply(preds_nc, 2, quantile, probs = .1),
                   u_ci80 = apply(preds_nc, 2, quantile, probs = .9),
                   l_ci95 = apply(preds_nc, 2, quantile, probs = .025),
                   u_ci95 = apply(preds_nc, 2, quantile, probs = .975))

preds_nc <- bind_cols(new_data, preds_nc)

preds_nc <- preds_nc |>
  mutate(time10 = time / 3,
         .after = time)

## marginal predictions (all the covariates held fixed at reference level)

preds_nc <-
  left_join(preds_nc,
            distinct(select(my_dt, time_date, time1)),
            by = "time1")

ggplot(data = preds_nc,
       aes(x = time_date,
           y = median)) +
  geom_line(aes(color = as.factor(capta))) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci95,
                     ymax = u_ci95)) +
  geom_linerange(aes(color = as.factor(capta),
                     ymin = l_ci80,
                     ymax = u_ci80),
                 lwd = 1.2) +
  geom_point(aes(fill = as.factor(capta)),
             pch = 21,
             color = 1, size = 3) +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Reports per 10k births",
       x = NULL,
       color = NULL)

readr::write_csv(select(preds_nc, time_date, median, l_ci80:u_ci95),
                 file = "data/data-exhibit-negative-control.csv")

## EXHIBIT 4
ggsave(filename = "figs/exhibit-negative-control.pdf",
       width = 12,
       height = 7,
       dpi = 300)
