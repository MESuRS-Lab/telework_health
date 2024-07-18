
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

# PLOT ####

result_all2 = result_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
         Tot_r = R + R_c) %>%
  group_by(DRF) %>%
  mutate(Tot_c = Tot_c/first(Tot_c)*100,
         Tot_r = Tot_r/first(Tot_r)*100) %>%
  ungroup %>%
  select(DRF, alpha, t_alpha, Tot_c, Tot_r) %>% 
  mutate(Tot_tot = Tot_c+Tot_r) %>%
  mutate(DRF = factor(DRF, levels = c("L-shaped",
                                      "U-shaped", "Inverted U-shaped",
                                      "Linear increasing shaped",
                                      "Linear decreasing shaped")))


low_thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == min(t_alpha))
high_thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha))
thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha) | t_alpha == min(t_alpha)) %>%
  ungroup %>%
  group_by(DRF, alpha) %>%
  summarise(min_t = min(t_alpha), max_t = max(t_alpha))


ggplot(result_all2) +
  geom_rect(data = thresholds, aes(xmin = min_t, xmax = max_t, ymin=0, ymax=180), alpha = 0.2) +
  geom_hline(yintercept = 100, linetype="dashed", colour="grey10", alpha = 0.3) +
  geom_line(aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 0.8, alpha = 0.4) +
  geom_line(aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 0.8, alpha = 0.4) +
  geom_line(data = result_all2 %>% filter(Tot_c <= 90), aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 1) +
  geom_line(data = result_all2 %>% filter(Tot_r <= 50), aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 1) +
  facet_nested_wrap(~DRF+alpha, ncol = 6) +
  scale_x_continuous(breaks = seq(0,90,20)) +
  scale_y_continuous(breaks = seq(0,180,30), limits = c(0,180)) +
  scale_colour_discrete(type = c("darkorange3","royalblue3")) +
  theme_bw() +
  labs(x = "Day of teleworking implementation start",
       y = "Relative cumulative incidence compared to no teleworking baseline(%)", colour = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

ggsave(here::here("figures", "fig4.png"), width = 9, height = 10)

