
library(dplyr)
library(reshape2)
library(ggplot2)

source(here::here("model", "NCD_function.R"))

# LINEAR DECREASING ####

NCD_rate_LD = NCD_function("LD", TRUE, 0.3)

# LINEAR INCREASING ####

NCD_rate_LI = NCD_function("LI", TRUE, 1.7)

# INVERTED U SHAPED ####

NCD_rate_IU = NCD_function("IU", TRUE, 1.7)

# U SHAPED ####

NCD_rate_US = NCD_function("US", TRUE, 0.3)

# L SHAPED ####

NCD_rate_LS = NCD_function("LS", TRUE, 0.3)


# Plot ####

all_curves = data.frame(alpha = seq(0,1,0.02)) %>%
  mutate(LD = NCD_rate_LD(alpha)) %>%
  mutate(LI = NCD_rate_LI(alpha)) %>%
  mutate(IU = NCD_rate_IU(alpha)) %>%
  mutate(US = NCD_rate_US(alpha)) %>%
  mutate(LS = NCD_rate_LS(alpha)) %>%
  melt(id.vars = "alpha") %>%
  mutate(variable = as.character(variable))

ggplot(all_curves, aes(alpha, value, colour = variable)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0,2,0.1)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(x = "Teleworking frequency", y = "Relative risk of non-communicable disease", colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11))

ggsave(here::here("figures", "suppfig2.png"))
