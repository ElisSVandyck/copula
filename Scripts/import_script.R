#Imports script
options(warn=-1)
library(tidyverse)
library(copula)
library(VineCopula)
library(HMMcopula)
options(warn=0)
bac_raw <- read.csv('bac.csv')
c_raw <- read.csv('c.csv')
gs_raw <- read.csv('gs.csv')
ms_raw <- read.csv('ms.csv')
jpm_raw <- read.csv('jpm.csv')
returns <- rbind(bac_raw %>% select(Date, Close) %>% mutate(Name = 'bac'),
      c_raw %>% select(Date, Close) %>% mutate(Name = 'c'),
      jpm_raw %>% select(Date, Close) %>% mutate(Name = 'jpm')) %>%
    tibble() %>% 
    pivot_wider(id_cols = Date, names_from = Name, values_from = Close) %>%
    arrange(Date) %>% 
    na.omit() %>%
    mutate(bac_returns = (bac-lag(bac))/lag(bac)) %>%
    mutate(c_returns = (c-lag(c))/lag(c)) %>%
    mutate(jpm_returns = (jpm-lag(jpm))/lag(jpm)) %>%
    na.omit() 
X <- returns %>% select(x1 = bac_returns,
                        x2 = c_returns,
                        x3 = jpm_returns)
U <- X %>% filter(x1 != 0 & x2 != 0 & x3 != 0) %>% #Filtration optional ~ 10% of sample omitted
    pobs() %>% 
    data.frame() %>% tibble() %>% 
    select(u1 = x1, u2 = x2, u3 = x3)



