# Long-term impacts of teleworking policies on occupational health following a pandemic from a private company perspective 
Code repository for the MESuRS group project

## Model folder

Contains the script for the model.

### Equations of the model

We modelled SARS-CoV-2 transmission in a company of $N$ employees using a compartmental model. In this model, employees can be susceptible to the respiratory disease $S$, exposed to the disease but not yet infectious $E$, infectious and asymptomatic $I_A$, exposed to the disease and infectious but presymptomatic $P$, infectious and symptomatic $I_S$, or recovered $R$. 


$$\frac{dS}{dt} = - \lambda S$$

$$\frac{dE}{dt} = \lambda S - \sigma E$$

$$\frac{dI_A}{dt} = p_A \sigma E - \gamma_A I_A$$

$$\frac{dP}{dt} = (1 - p_A) \sigma E - \rho P $$

$$\frac{dI_S}{dt} = \rho P - \gamma_S I_S$$

$$\frac{dR}{dt} = \gamma_A I_A + \gamma_S I_S$$


In this compartmental model, infected individuals can develop symptoms with probability $1-p_A$ but their incubation period (here, time from infection to time of the onset of infectiousness) is the same whether they develop symptoms or not and is equal to $\frac{1}{\sigma}$. We also consider that for all individuals there is a pre-symptomatic infectious phase, with a duration equal to $\frac{1}{\rho}$. Infectious and symptomatic individuals are infectious for $\frac{1}{\gamma_S}$ days and are assumed to be on medical leave from their symptom/infectiousness onset to their recovery. Thereby, these individuals do not contribute to the propagation of the epidemic within the company. Infectious and asymptomatic individuals are infectious for $\frac{1}{\gamma_A}$ days and are responsible for the spread of the disease within the company.

The transmission rate $\lambda$ is divided into three terms as follows:

$$\lambda = \frac{5}{7} (1-\alpha) \frac{\beta (\nu I_A + P)}{N-I_S} + \frac{5}{7} \alpha \epsilon \lambda_v + \frac{2}{7} \lambda_v$$

Where $\beta$ is expressed using SARS-CoV-2 $R_0$ that we derived from the next-generation matrix for a frequency-dependent model:

$$\beta = \frac{R_0 \rho \gamma_A}{(1-\alpha) [(1 - p_A) \gamma_A + \rho \nu p_A]}$$

$\alpha$ is the proportion of employees teleworking, $\nu$ is the coefficient of relative infectivity of asymptomatic cases ($I_A$) compared to symptomatic cases ($P$ and $I_S$), $\lambda_v$ is the transmission rate from the community, and $\epsilon$ is a coefficient reducing the transmission from the community on teleworking days.

In addition to the transmission process of the infectious disease, we modelled the number of individuals who will ultimately develop a chronic disease following exposure to teleworking. To do so, we stratified the compartmental model into two populations, one population that will not develop a chronic disease, and a second that will develop a chronic disease. Given the different time scales of occurrence of the chronic disease and the infectious disease, we assume that the two populations mix homogeneously between them.

### Parameterization of the model

According to [Wu et al., 2022](https://doi.org/10.1001/jamanetworkopen.2022.28008), the pooled incubation period for all SARS-CoV-2 variants is 6.57 days. In this case, $\sigma = \frac{1}{6.57}$. 

| SARS-CoV-2 variant | Mean incubation period (in days) |
| :----------------: | :------------------------------: |
| Alpha              | 5.00                             |
| Beta               | 4.50                             |
| Delta              | 4.41                             |
| Omicron            | 3.42                             |

[Koelle et al., 2022](https://doi.org/10.1126/science.abm4915) have reviewed how our understanding of the transmission of SARS-CoV-2 has changed over the pandemic mostly using modelling studies. The first studies in Wuhan estimated that the $R_0$ lied between 2 and 4 before the lockdown. In their review, [Dhungel et al., 2022]( https://doi.org/10.3390/ijerph191811613) have estimated a pooled $R_0$ of 2.66 from studies published in the early months of the pandemic. 

[Buitrago-Garcia et al., 2020](https://doi.org/10.1371/journal.pmed.1003346) estimated that the proportion of asymptomatic cases is on average equal to 20%. This would lead to $p_A = 0.20$. They also estimated that asymptomatic cases are less infectious than symptomatic cases by 65% leading to $\nu = 0.35$.

Need to find estimates for $\gamma_A$ and $\gamma_S$.

## Analysis folder

Contains scripts used to run the model and generate figures.

## Figures folder

Contains figures produced by analysis scripts.

