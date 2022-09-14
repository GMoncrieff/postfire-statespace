library(dplyr)
library(readr)
library(tidyr)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(stringr)
library(ggplot2)

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

obsmat <- pivot_wider(dplyr::select(full_sub,c(pid,tid,y)), names_from = tid, values_from = y) %>%
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

firenew <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


############STAN


set.seed(84)

postfire_data <- list(N = nrow(full_sub),
                      M = length(firenew),
                      TT = ncol(firemat)-1,
                      J=nrow(firemat),
                      y = y_fill,
                      fire=as.matrix(firemat[,2:ncol(firemat)]),
                      firenew=firenew,
                     # envg1 = env$env1,
                    #  envg2 = env$env2,
                      z0=pmin(pmax(rnorm(nrow(firemat),0.4,0.1),0.1),0.8))

model <- cmdstan_model('models/postfire_ss_nomiss_noenv_gen.stan', compile = TRUE)
#model$print()

fit_mcmc <- model$sample(
  data = postfire_data,
  seed = 12345,
  refresh=50,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 10000,
  iter_sampling = 10000
)

mcmc_trace(fit_mcmc$draws(c("gamma")))
sumfit <- fit_mcmc$summary(.args = list(na.rm=TRUE))

sumfit <- sumfit %>%
  dplyr::select(c(variable,median,q5,q95))

zfit <- sumfit %>%
  dplyr::filter(grepl('z\\[', variable)) %>%
  dplyr::mutate(znum=str_sub(variable, 3, -2)) %>%
  tidyr::separate(znum,into=c('pix','time'),sep=',') %>%
  dplyr::mutate(pix=as.numeric(pix),time=as.numeric(time)) %>%
  dplyr::filter(pix==16) %>%
  dplyr::mutate(flag="fit")

zsim <- sumfit %>%
  dplyr::filter(grepl('znew', variable)) %>%
  dplyr::mutate(time=str_sub(variable, 6, -2)) %>%
  dplyr::mutate(time=as.numeric(time)) %>%
  dplyr::filter(time!=1) %>%
  dplyr::mutate(time=time+max(zfit$time)-1) %>%
  dplyr::mutate(pix=16) %>%
  dplyr::relocate(pix, .before = time) %>%
  dplyr::mutate(flag="pred")

zdat <- rbind(zfit,zsim)


obsplot <- data.frame(y=t(obsmat[16,2:ncol(obsmat)]),time=seq(1:(ncol(obsmat)-1)))

  
ggplot() + 
  geom_ribbon(data=zdat,aes(ymin = q5, ymax = q95,x=time), fill = "grey70") +
  geom_line(data=zdat,aes(y = median,x=time)) +
  geom_point(data=obsplot,aes(y=y,x=time),col='red',alpha=0.2) +
  annotate("rect", xmin = 380, xmax = 419, ymin = 0, ymax = 1,
           alpha = .2, fill = "grey") +
  theme_bw() +
  labs(x='time step',y='ndvi')
