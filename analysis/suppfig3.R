
# this script runs example scenarios of teleworking proportion
# you can define the scenarios to try by changing this vector:
alpha_to_try = c(0, 0.2, 0.4, 0.6, 0.8, 1)

# the community FOI is displayed in the final plot for convenience

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(openxlsx)
library(here)
library(lubridate)

source(here("Model", "model.R"))
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
baseline_NCD = 0.000089   # baseline NCD rate

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha, t_alpha=t_alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s, baseline_NCD=baseline_NCD)

# community force of infection
# uses approxfun() to generate an interpolating function, passed to the model function
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

#with data will look something like:
# commu_FOI = approxfun(data_time_points, data_values)

# rate of chronic disease linked to telework
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using a constant max rate multiplied by alpha

result_NC_all = data.frame()
result_RR_NCD = data.frame(alpha = alpha_to_try)

# LINEAR DECREASING ####

#assuming 20 days of work per months, 8h per day
NCD_rate = NCD_function("LS", TRUE, 0.3)
DRF = "LS"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# INVERTED U SHAPED ####

NCD_rate = NCD_function("IU", TRUE, 1.7)
DRF = "IU"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# U SHAPED ####

NCD_rate = NCD_function("US", TRUE, 0.3)
DRF = "US"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# LINEAR DECREASING ####

NCD_rate = NCD_function("LD", TRUE, 0.3)
DRF = "LD"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# LINEAR INCREASING ####

NCD_rate = NCD_function("LI", TRUE, 1.7)
DRF = "LI"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}


result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))

result_NC_all = result_NC_all %>%
  mutate(DRF = factor(DRF, levels = c("LS", "US", "IU", "LD", "LI")))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# PLOT ####

result_NC_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c) %>%
  ggplot() +
  geom_point(aes(DRF, Tot_c, shape = DRF), size = 3) +
  facet_grid(cols = vars(alpha)) +
  theme_bw() +
  guides(shape = "none") +
  labs(x = "Exposure-response curve", y = "Cumulative number of individuals\neventually developing a NCD", col = "") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(size=11))

ggsave(here::here("figures", "suppfig3.png"), width = 8, height = 4)

