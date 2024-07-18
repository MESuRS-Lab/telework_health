
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# palette set to match supplementary figure 2
pal = brewer.pal(5, "Set1")

# LINEAR DECREASING ####

#assuming 20 days of work per months, 8h per day
NCD_rate_LD = data.frame(freqs = c(0, 0.25/5, 1/5, 5/5),
                         vals = c((1/(1+exp(-(-3.33)))) / (1/(1+exp(-(-3.33)))),
                                  (1/(1+exp(-(-3.33-0.571)))) / (1/(1+exp(-(-3.33)))),
                                  (1/(1+exp(-(-3.33-0.997)))) / (1/(1+exp(-(-3.33)))),
                                  (1/(1+exp(-(-3.33-1.233)))) / (1/(1+exp(-(-3.33))))),
                         DRF = "LS")

# INVERTED U SHAPED ####

NCD_rate_IU = data.frame(freqs = c(0,0.2,0.5,1),
                         vals = c(49.9/49.9, 52.8/49.9, 53.2/49.9, 47.5/49.9),
                         DRF = "IU")


# U SHAPED ####

NCD_rate_US = data.frame(freqs = c(0,0.5,1),
                      vals = c(2.33/2.33, 1.71/2.33, 2.08/2.33), 
                      DRF = "US")

rbind(NCD_rate_LD, NCD_rate_IU, NCD_rate_US) %>%
  ggplot(aes(freqs, vals, colour = DRF)) +
  geom_point(size = 4) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0,2,0.1)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_color_discrete(type = pal[c(1,4,5)]) +
  theme_bw() +
  labs(x = "Teleworking frequency", y = "Relative risk of non-communicable disease", colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11))

ggsave(here::here("figures", "suppfig1.png"))

