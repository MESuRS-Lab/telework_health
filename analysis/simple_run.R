
# this script allows you to quickly run the model and plot the output

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(openxlsx)
library(here)
library(lubridate)

source(here::here("Model", "model.R"))
source(here("Model", "NCD_function.R"))

R0 = 2.66              # Basic reproduction number
alpha = 0.1           # proportion of teleworking
t_alpha = -1          # activation time for teleworking (default -1 = always on)
nu = 0.35             # relative force of infection of asymptomatic cases
epsilon = 0.21         # relative force of infection during teleworking
sigma = 1/6.57         # progression rate from exposed to infectious
rho = 1/1.5           # progression rate from pre-symptomatic to symptomatic
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics
baseline_NCD = 0.00014   # baseline NCD rate

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha, t_alpha = t_alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s,
          baseline_NCD = baseline_NCD)

# community force of infection
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using sin function shifted by 300 to align with start of simulation
# commu_FOI = approxfun(c(0:1000), 0.001/2*(sin((2*pi/200)*c(0:1000)+300)+1))
df <- read.xlsx(here("data", "extract_results_Laura.xlsx"), rows = c(1:6189)) %>%
  mutate(Time = as_date(Time, origin = "1899-12-30 UTC"))

df <- df %>%
  select(Time, AgeGroup, Exposed, Susceptible) %>%
  filter(AgeGroup %in% c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                         "50-54", "55-59", "60-64")) %>%
  group_by(Time) %>%
  summarise(Exposed = sum(Exposed), Susceptible = sum(Susceptible)) %>%
  mutate(Proba = Exposed/Susceptible) %>%
  mutate(Taux = -log(1-Proba)) %>%
  ungroup %>%
  filter(!is.nan(Taux)) %>%
  filter(Time >= as_date("2020-09-01") & Time <= as_date("2020-12-01"))

commu_FOI = approxfun(c(0:(nrow(df)-1)), df$Taux)

# rate of chronic disease linked to telework
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using a constant max rate multiplied by alpha
NCD_rate = NCD_function("IU")

result = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                             commu_FOI = commu_FOI, NCD_rate = NCD_rate))

p1 = result %>%
  mutate(
    S = S + S_c, 
    E = E + E_c, 
    Ia = Ia + Ia_c,
    P = P + P_c,
    Is = Is + Is_c,
    R = R + R_c
    ) %>%
  select(!contains("_c")) %>%
  melt(id.vars="time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time (days)", y = "", col = "", title = "Number of individuals")

p2 = result %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c) %>%
  ggplot() +
  geom_line(aes(time, Tot_c)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time (days)", y = "", col = "", title = "Number of individuals that will\ndevelop a NCD")

ggarrange(p1, p2, common.legend = T, legend = "bottom", align = "hv")

#ggsave(here::here("figures", "example_simple_run.png"))
