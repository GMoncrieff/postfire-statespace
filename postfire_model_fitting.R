library(dplyr)
library(readr)
library(tidyr)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(stringr)
library(ggplot2)

set.seed(1234)


# read raw data -----------------------------------------------------------

#data come from dimensions plots
pyro_test <- read_csv("pyro_test.csv")


# data prep ---------------------------------------------------------------

#selection only those plots that have a full time series
nogaps<- pyro_test %>% 
  group_by(id) %>% 
  summarise(n = n()) %>%
  dplyr::filter(n==380)

#create a fire flag -when age/DA decreases that means we had a fire
#create a date variable - I didn't include data in this dataset. but we know that there are 16 days between obs, so we just use this
#then we convert this date variable to radians
fullts <- pyro_test %>%
  dplyr::filter(id %in% nogaps$id) %>% 
  group_by(id) %>% 
  mutate(firelag = DA - lag(DA)) %>%
  mutate(firelead = lead(firelag)) %>%
  mutate(fire = as.numeric((firelag<0))) %>%
  ungroup() %>%
  mutate(ND= pmax(0,ND)) %>%
  dplyr::select(c(id=...1,pid=id,y=ND,map,tmin07,fire)) %>%
  group_by(pid) %>%
  mutate(minid=min(id)) %>%
  mutate(tid=id-minid+1) %>%
  mutate(y = replace(y, y==0, NA)) %>%
  mutate(fire = replace_na(fire, 0)) %>%
  dplyr::select(-minid) %>%
  ungroup() %>%
  dplyr::mutate(doy_rad=(((((tid-1)*16)+1)%%365.25)/365.25)*2*pi)

#subsample dataset to 20 pixels
pix <- unique(fullts$pid) %>% sample(20)
full_sub  <- fullts %>%
  dplyr::filter(pid %in% pix) %>%
  mutate(rid = row_number()) %>%
  arrange(pid,tid)

#get locations of observed and missing data
y_obs_i <- which(!is.na(full_sub$y))
y_obs <- full_sub$y[y_obs_i]
y_na_i <- which(is.na(full_sub$y))
y_na <- full_sub$y[y_na_i]
#y_fill <- full_sub$y
#y_fill[y_na_i] = 0.3

#fires as one long vector
firevec <- full_sub$fire

#fires as matrix - col = date, row = pixel
firemat <- pivot_wider(dplyr::select(full_sub,c(pid,tid,fire)), names_from = tid, values_from = fire) %>%
  arrange(pid)

#ndvi obs as matrix - col = date, row = pixel
obsmat <- pivot_wider(dplyr::select(full_sub,c(pid,tid,y)), names_from = tid, values_from = y) %>%
  arrange(pid)

#day of year as matrix - col = date, row = pixel
doymat <- pivot_wider(dplyr::select(full_sub,c(pid,tid,doy_rad)), names_from = tid, values_from = doy_rad) %>%
  arrange(pid)

#standardize env variables
standard = function(x){
  z <- (x - mean(x)) / sd(x)
  return( z)
}
#get env varaible per pixel
env <- full_sub %>%
  group_by(pid) %>%
  summarise(env1 = first(map),env2 = first(tmin07)) %>%
  mutate(env1 = standard(env1),env2 = standard(env2)) %>%
  arrange(pid)

#made-up sequence of fires to use for forecasting
firenew <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#sequence of doy to use for forecasting
doy_rad_new <- seq(from=max(full_sub$tid),by=1,length.out=length(firenew))
doy_rad_new <- (((((doy_rad_new-1)*16)+1)%%365.25)/365.25)*2*pi


# stan modelling ----------------------------------------------------------

set.seed(84)
#example pixel id
eg_id = 4
#create stan data obj
postfire_data <- list(N = nrow(full_sub),
                      M = length(firenew),
                      TT = ncol(firemat)-1,
                      J=nrow(firemat),
                      N_obs=length(y_obs_i),
                      N_mis=length(y_na_i),
                      ii_obs = y_obs_i,
                      ii_mis = y_na_i,
                      y_obs=y_obs,
                      #y = y_fill,
                      fire=firevec,
                      firenew=firenew,
                      doy_rad=as.matrix(doymat[,2:ncol(doymat)]),
                      doy_rad_new=doy_rad_new,
                      envg1 = env$env1,
                      envg2 = env$env2,
                      z0=pmin(pmax(rnorm(nrow(firemat),0.4,0.1),0.1),0.8),
                      eg_id = eg_id)
#compil model
model <- cmdstan_model('models/postfire_ss_obserr.stan', compile = TRUE)

#fit model
fit_mcmc <- model$sample(
  data = postfire_data,
  seed = 12345,
  refresh=50,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 200,
  iter_sampling = 200,
  max_treedepth = 15
)


# model evaluation --------------------------------------------------------

#examine trace for parameter
#mcmc_trace(fit_mcmc$draws(c("gamma")))

#summarize fit
sumfit <- fit_mcmc$summary(.args = list(na.rm=TRUE))
sumfit <- sumfit %>%
  dplyr::select(c(variable,median,q5,q95))

#fitted states
zfit <- sumfit %>%
  dplyr::filter(grepl('z\\[', variable)) %>%
  dplyr::mutate(znum=str_sub(variable, 3, -2)) %>%
  tidyr::separate(znum,into=c('pix','time'),sep=',') %>%
  dplyr::mutate(pix=as.numeric(pix),time=as.numeric(time)) %>%
  dplyr::filter(pix==eg_id) %>%
  dplyr::mutate(flag="fit")

#forecast states
zsim <- sumfit %>%
  dplyr::filter(grepl('znew', variable)) %>%
  dplyr::mutate(time=str_sub(variable, 6, -2)) %>%
  dplyr::mutate(time=as.numeric(time)) %>%
  dplyr::filter(time!=1) %>%
  dplyr::mutate(time=time+max(zfit$time)-1) %>%
  dplyr::mutate(pix=eg_id) %>%
  dplyr::relocate(pix, .before = time) %>%
  dplyr::mutate(flag="pred")

zdat <- rbind(zfit,zsim)
#observations
obsplot <- data.frame(y=t(obsmat[eg_id,2:ncol(obsmat)]),time=seq(1:(ncol(obsmat)-1)))

#plot!
ggplot() + 
  geom_ribbon(data=zdat,aes(ymin = q5, ymax = q95,x=time), fill = "grey70") +
  geom_line(data=zdat,aes(y = median,x=time)) +
  geom_point(data=obsplot,aes(y=y,x=time),col='red',alpha=0.2) +
  annotate("rect", xmin = 380, xmax = 419, ymin = 0, ymax = 1,
           alpha = .2, fill = "grey") +
  theme_bw() +
  labs(x='time step',y='ndvi')
