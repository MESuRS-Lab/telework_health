
library(ggplot2)
library(scales)
library(openxlsx)
library(here)
library(dplyr)
library(lubridate)

df <- read.xlsx(here("data", "extract_results_Laura.xlsx"), rows = c(1:6189)) %>%
  mutate(Time = as_date(Time, origin = "1899-12-30 UTC"))


df %>%
  group_by(Time) %>%
  summarise(Exposed = sum(Exposed)) %>%
  # filter(Time > 44050 & Time < 44150) %>%
  ggplot() +
  geom_line(aes(Time, Exposed)) +
  theme_bw()

ggsave("all_covid.png")

df = df %>%
  select(Time, AgeGroup, Exposed, Susceptible) %>%
  filter(AgeGroup %in% c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                         "50-54", "55-59", "60-64")) %>%
  group_by(Time) %>%
  summarise(Exposed = sum(Exposed), Susceptible = sum(Susceptible)) %>%
  mutate(Proba = Exposed/Susceptible) %>%
  mutate(Taux = -log(1-Proba)) %>%
  ungroup %>%
  filter(!is.nan(Taux)) %>%
  # filter(Time <= as_date("2020-05-30"))
  filter(Time >= as_date("2020-09-01") & Time <= as_date("2020-12-01"))


ggplot(df) +
  geom_col(aes(Time, Exposed), alpha = 0.6) +
  geom_line(aes(Time, Taux*10^(7.6)), linewidth = 1) +
  scale_y_continuous(breaks = seq(0, 100000, 20000), labels = label_number(),
                     sec.axis = sec_axis(transform=~./10^(7.6),
                                         name = bquote(lambda[C]),
                                         breaks = seq(0, 0.003, 0.0005),
                                         labels = label_number())) +
  theme_bw() +
  labs(x = "Time (days)", y = "Incidence") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

ggsave(here("figures","community_FOI.png"))
