library(ggplot2)

##--- Exhibit 4 ----

dt_4 <-
  readr::read_csv("data/data-exhibit-4.csv")

ggplot(data = dt_4,
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

##--- Exhibit 6 ----

dt_6 <-
  readr::read_csv("data/data-exhibit-6.csv")

ggplot(data = dt_6,
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
  labs(y = "Placements per 100 reports",
       x = NULL,
       color = NULL)

##--- Negative Control ----

dt_nc <-
  readr::read_csv("data/data-exhibit-negative-control.csv")

ggplot(data = dt_nc,
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
