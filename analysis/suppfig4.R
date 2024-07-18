
#add prcc checking impact of parameters on relative reduction in cumulative incidence?
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(openxlsx)
library(here)
library(epiR)
library(lubridate)

source(here::here("Model", "model.R"))
source(here("Model", "NCD_function.R"))

R0 = 2.66              # Basic reproduction number
alpha = 0.5           # proportion of teleworking
t_alpha = -1          # activation time for teleworking (default -1 = always on)
nu = 0.35             # relative force of infection of asymptomatic cases
epsilon = 0.5         # relative force of infection during teleworking
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


## sensitivity ###############

n_runs = 400

all_results_LS = data.frame(R0=runif(n_runs, 2, 4),
                            alpha=0.5,
                            # alpha=runif(n_runs, 0.5, 1),
                            t_alpha=t_alpha,
                            nu=runif(n_runs, 0.1, 1.27),
                            epsilon=runif(n_runs, epsilon*0.8, epsilon*1.2),
                            sigma=runif(n_runs, 1/18.87, 1/1.80),
                            rho=rho,
                            prop_a=runif(n_runs, 0.17, 0.25),
                            gamma_a=gamma_a,
                            gamma_s=gamma_s,
                            baseline_NCD = baseline_NCD,
                            thresh=runif(n_runs, 0.05, 0.7),
                            cumul_ID = 0,
                            cumul_NCD = 0)

all_results_US = all_results_IUS = all_results_LI = all_results_LD = all_results_LS

all_results_IUS$shape = "IU"
all_results_US$shape = "US"
all_results_LS$shape = "LS"
all_results_LI$shape = "LI"
all_results_LD$shape = "LD"

for(i in 1:nrow(all_results_LS)){
  
  if(i %% round(nrow(all_results_LS)/10) == 0) cat(i/round(nrow(all_results_LS))*100, "% done\n")
  
  param = as.vector(all_results_LS[i,-c((ncol(all_results_LS)-1) , ncol(all_results_LS))])
  
  #L SHAPED
  NCD_rate = NCD_function("LS", TRUE, 1-param$thresh)
  results = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                                commu_FOI = commu_FOI, NCD_rate = NCD_rate))
  results = results %>%
    filter(time == max(time)) %>%
    mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
           Tot_r = R + R_c)
  
  all_results_LS$cumul_ID[i] = results$Tot_r
  all_results_LS$cumul_NCD[i] = results$Tot_c
  
  #U SHAPED
  NCD_rate = NCD_function("US", TRUE, 1-param$thresh)
  results = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                                commu_FOI = commu_FOI, NCD_rate = NCD_rate))
  results = results %>%
    filter(time == max(time)) %>%
    mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
           Tot_r = R + R_c)
  
  all_results_US$cumul_ID[i] = results$Tot_r
  all_results_US$cumul_NCD[i] = results$Tot_c
  
  #INVERTED U SHAPED
  NCD_rate = NCD_function("IU", TRUE, 1+param$thresh)
  results = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                                commu_FOI = commu_FOI, NCD_rate = NCD_rate))
  results = results %>%
    filter(time == max(time)) %>%
    mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
           Tot_r = R + R_c)
  
  all_results_IUS$cumul_ID[i] = results$Tot_r
  all_results_IUS$cumul_NCD[i] = results$Tot_c
  
  #LINEAR INCREASE
  NCD_rate = NCD_function("LI", TRUE, 1+param$thresh)
  results = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                                commu_FOI = commu_FOI, NCD_rate = NCD_rate))
  results = results %>%
    filter(time == max(time)) %>%
    mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
           Tot_r = R + R_c)
  
  all_results_LI$cumul_ID[i] = results$Tot_r
  all_results_LI$cumul_NCD[i] = results$Tot_c
  
  #LINEAR DECREASE
  NCD_rate = NCD_function("LD", TRUE, 1-param$thresh)
  results = as.data.frame(lsoda(Init.cond, Time, model_function, param,
                                commu_FOI = commu_FOI, NCD_rate = NCD_rate))
  results = results %>%
    filter(time == max(time)) %>%
    mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
           Tot_r = R + R_c)
  
  all_results_LD$cumul_ID[i] = results$Tot_r
  all_results_LD$cumul_NCD[i] = results$Tot_c
  
}

all_results_LS_m = all_results_LS %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_ID) %>%
  epi.prcc() %>%
  mutate(shape = "LS") %>%
  rename(param = var) %>%
  mutate(cor = "ID")
all_results_IUS_m = all_results_IUS %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_ID) %>%
  epi.prcc() %>%
  mutate(shape = "IU") %>%
  rename(param = var) %>%
  mutate(cor = "ID")
all_results_US_m = all_results_US %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_ID) %>%
  epi.prcc() %>%
  mutate(shape = "US") %>%
  rename(param = var) %>%
  mutate(cor = "ID")
all_results_LI_m = all_results_LI %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_ID) %>%
  epi.prcc() %>%
  mutate(shape = "LI") %>%
  rename(param = var) %>%
  mutate(cor = "ID")
all_results_LD_m = all_results_LD %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_ID) %>%
  epi.prcc() %>%
  mutate(shape = "LD") %>%
  rename(param = var) %>%
  mutate(cor = "ID")

all_results_LS_m2 = all_results_LS %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_NCD) %>%
  epi.prcc() %>%
  mutate(shape = "LS") %>%
  rename(param = var) %>%
  mutate(cor = "NCD")
all_results_IUS_m2 = all_results_IUS %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_NCD) %>%
  epi.prcc() %>%
  mutate(shape = "IU") %>%
  rename(param = var) %>%
  mutate(cor = "NCD")
all_results_US_m2 = all_results_US %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_NCD) %>%
  epi.prcc() %>%
  mutate(shape = "US") %>%
  rename(param = var) %>%
  mutate(cor = "NCD")
all_results_LI_m2 = all_results_LI %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_NCD) %>%
  epi.prcc() %>%
  mutate(shape = "LI") %>%
  rename(param = var) %>%
  mutate(cor = "NCD")
all_results_LD_m2 = all_results_LD %>%
  select(R0, nu, epsilon, sigma, prop_a, thresh, cumul_NCD) %>%
  epi.prcc() %>%
  mutate(shape = "LD") %>%
  rename(param = var) %>%
  mutate(cor = "NCD")


all_results = rbind(all_results_LS_m, all_results_US_m, all_results_IUS_m, all_results_LI_m, all_results_LD_m,
                    all_results_LS_m2, all_results_US_m2, all_results_IUS_m2, all_results_LI_m2, all_results_LD_m2)

ggplot(all_results) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, group = cor, 
                      ymin = lower, ymax = upper, colour = cor),
                  position=position_dodge(width=0.2), size = 0.2, linewidth = 0.8) +
  facet_grid(rows = vars(shape)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient") +
  scale_x_discrete(labels = c(bquote(epsilon),
                              bquote(nu),
                              bquote(p[a]),
                              bquote(R0),
                              bquote(sigma),
                              bquote(omega))) +
  scale_colour_discrete(type = c("darkorange3","royalblue3")) +
  coord_cartesian(clip = "off", ylim = c(-1,1))

ggsave(here::here("figures", "suppfig4.png"))
