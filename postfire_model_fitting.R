library(dplyr)
library(readr)
library(tidyr)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(stringr)
library(ggplot2)
library(lubridate)

set.seed(1234)


# read raw data -----------------------------------------------------------

#data come from dimensions plots
fullts <- read_csv("example_postfire_data.csv")

#id = obs id
#pid = unique pixel id
#y = MODIS NDVI
#map = mean annual precipitation in mm
#tmin07 = july min temp (deg C)
#fire = did a fire occur?
#tid = temporal id (date)
#doy rad = day of year in radians

# data prep ---------------------------------------------------------------

#subsample dataset to 20 pixels
pix <- unique(fullts$pid) %>% sample(40)
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
model <- cmdstan_model('models/postfire_ss.stan', compile = TRUE)

#fit model
fit_mcmc <- model$sample(
  data = postfire_data,
  seed = 12345,
  refresh=50,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 1000
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

#dates of fires
fireplot <- data.frame(y=t(firemat[eg_id,2:ncol(firemat)]),time=seq(1:(ncol(firemat)-1)))
ffire <- nrow(fireplot)+which(firenew==1)
fireplot <- fireplot %>% filter(y==1)
fireplot <- rbind(fireplot,c(1,ffire))

#dates for labels
dates = dmy("1/1/2000") + days((zdat$time*16)-16)
dateslab = dmy("1/1/2000") + days((c(1, 100,200,300,400)*16)-16)

#plot!
ggplot() + 
  geom_ribbon(data=zdat,aes(ymin = q5, ymax = q95,x=time), fill = "grey70") +
  geom_line(data=zdat,aes(y = median,x=time)) +
  geom_vline(data= fireplot, aes(xintercept = time),lty=2,alpha=0.5) +
  geom_point(data=obsplot,aes(y=y,x=time),col='red',alpha=0.2) +
  annotate("rect", xmin = 380, xmax = 440, ymin = 0, ymax = 0.5,
           alpha = .2, fill = "grey") +
  theme_bw() +
  ylim(0,0.5) +
  scale_x_continuous(breaks = c(1,100,200,300,400), labels = dateslab) +
  labs(x='date',y='NDVI')

ggsave('example_ts.png',scale=0.6, width=10, height=6)


