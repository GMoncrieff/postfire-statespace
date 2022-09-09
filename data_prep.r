library(dplyr)
library(readr)

pyro_test <- read_csv("pyro_test.csv")

nogaps<- pyro_test %>% 
  group_by(id) %>% 
  summarise(n = n()) %>%
  dplyr::filter(n==380)

fullts <- pyro_test %>%
  dplyr::filter(id %in% nogaps$id) %>% 
  group_by(id) %>% 
  mutate(firelag = DA - lag(DA)) %>%
  mutate(firelead = lead(firelag)) %>%
  mutate(fire = as.numeric((firelag<0))) %>%
  ungroup() %>%
  dplyr::select(c(id=...1,pid=id,y=ND,map,tmin07,fire)) %>%
  group_by(pid) %>%
  mutate(minid=min(id)) %>%
  mutate(tid=id-minid+1) %>%
  mutate(y = replace(y, y==0, NA)) %>%
  mutate(fire = replace_na(fire, 0)) %>%
  dplyr::select(-minid) %>%
  ungroup()

pix <- unique(fullts$pid) %>% sample(20)

full_sub  <- fullts %>%
  dplyr::filter(pid %in% pix) %>%
  mutate(rid = row_number()) %>%
  arrange(pid,tid)

y_obs_i <- which(!is.na(full_sub$y))
y_obs <- full_sub$y[y_obs_i]
y_na_i <- which(is.na(full_sub$y))
y_na <- full_sub$y[y_na_i]
y_fill <- full_sub$y
y_fill[y_na_i] = 0.3

firemat <- pivot_wider(dplyr::select(full_sub,c(pid,tid,fire)), names_from = tid, values_from = fire) %>%
  arrange(pid)

standard = function(x){
  z <- (x - mean(x)) / sd(x)
  return( z)
}

env <- full_sub %>%
  group_by(pid) %>%
  summarise(env1 = first(map),env2 = first(tmin07)) %>%
  mutate(env1 = standard(env1),env2 = standard(env2)) %>%
  arrange(pid)



############STAN
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyr)
library(stringr)
library(knitr)
library(ggplot2)
color_scheme_set("brightblue")
set.seed(83)

postfire_data <- list(N = nrow(postfire),
                      J= nrow(postfire_grp),
                      nd=postfire$nd,
                      age= postfire$age,
                      pid=postfire$pid,
                      envg1 = postfire_grp$envg1,
                      envg2 = postfire_grp$envg2)

model <- cmdstan_model('models/postfire_ss_nomiss.stan', compile = TRUE)
model$print()






  
  
