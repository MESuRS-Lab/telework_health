
library(lubridate)
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggh4x)
library(here)
library(openxlsx)

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
baseline_NCD = 0.000089   # baseline NCD rate

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s, baseline_NCD=baseline_NCD)

# community force of infection
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

result_all = data.frame()


for(alpha in seq(0, 1, 0.2)){
  
  param["alpha"] = alpha
  
  for(t_alpha in c(0:90)){
    
    param["t_alpha"] = t_alpha
    
    # L SHAPED ####
    NCD_rate = NCD_function("LS", TRUE, 0.3)
    DRF = "L-shaped"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # U SHAPED ####
    NCD_rate = NCD_function("US", TRUE, 0.3)
    DRF = "U-shaped"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # INVERTED U SHAPED ####
    NCD_rate = NCD_function("IU", TRUE, 1.7)
    DRF = "Inverted U-shaped"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # LINEAR INCREASING SHAPED ####
    NCD_rate = NCD_function("LI", TRUE, 1.7)
    DRF = "Linear increasing shaped"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # LINEAR DECREASING SHAPED ####
    NCD_rate = NCD_function("LD", TRUE, 0.3)
    DRF = "Linear decreasing shaped"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
  }
  
}

DALY_ID = 0.04327715*16
DALY_NCD = 0.11447948*365*2

result_all2 = result_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
         Tot_r = R + R_c) %>%
  mutate(Tot_c_DALY = Tot_c*DALY_NCD,
         Tot_r_DALY = Tot_r*DALY_ID) %>%
  # group_by(DRF) %>%
  # mutate(Tot_c = Tot_c/first(Tot_c)*100,
  #        Tot_r = Tot_r/first(Tot_r)*100) %>%
  # ungroup %>%
  select(DRF, alpha, t_alpha, Tot_c, Tot_r, Tot_c_DALY, Tot_r_DALY) %>% 
  mutate(Tot_tot = Tot_c+Tot_r,
         Tot_tot_DALY = Tot_c_DALY+Tot_r_DALY)
  # mutate(DRF = factor(DRF, levels = c("L-shaped",
  #                                     "U-shaped", "Inverted U-shaped")))

# stacked barplot of the change in NCD and ID DALY with increase telework frequency, for each DRF

result_all2 %>%
  filter(t_alpha == 0) %>%
  melt(id.vars = c("DRF", "alpha"), measure.vars = c("Tot_c_DALY", "Tot_r_DALY")) %>%
  mutate(variable = recode(variable, "Tot_c_DALY" = "NCD", "Tot_r_DALY" = "ID")) %>%
  ggplot() +
  facet_wrap(~DRF) +
  geom_bar(aes(x = alpha, y = value, fill = variable), stat = "identity", position = "stack") +
  theme_bw() +
  labs(x = "Proportion of teleworking", y = "Total DALY (stacked)", fill = "Type of DALY") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        strip.text = element_text(size = 12))

ggsave(here::here("figures","DALY_plot.png"), height = 5, width = 9)
