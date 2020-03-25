rm(list=ls(all=TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# require(tidyverse)
require(ggfortify)
require(lamW)
require(GGally)
require(ggcorrplot)
require(pracma)
require(nleqslv)
source("./ggplot_mytheme.R")

load( "./dataestimate.Rdata" )

maintext_figs <- c("gut1", "river", "sludge")

beta_fun <- function( eta , ...){
  pars <- list(...)
  return( pars[1] )
}


### NUMBER OF SPECIES VS NREADS
predicted_sp <- function(N, mu, sigma, beta, sp_tot ){
  bb = beta
  eta = rnorm( 10^6, mean = mu, sd = sigma )
  p0 = mean( (1.+N*exp(eta)/bb)^-bb )
  return( sp_tot*(1-p0) )
}

### CAMBIARE CON SIMPSON?
predicted_shannonindex <- function(N, mu, sigma, beta, sp_tot ){
  s <- c()
  for(i in 1:100){
    eta = rnorm( stot, mean = mu, sd = sigma)   
    p0 =  mean( (beta/(beta+N*exp(eta))  )^beta ) 

    prob = (beta/(beta+N*exp(eta)))
    k <- rnbinom(length(prob), size = beta, prob = prob)
    k <- 1.*k[k>0] / N
    if( -sum( k*log(k) ) > 0 ){
      s <- c(s,  -sum( k*log(k) ) )
    }
  }
  return( mean(s) )
}

predicted_reads <- function(N, mu, sigma, beta, sp_tot ){
  eta = rnorm( 10^6, mean = mu, sd = sigma  )
  p0 = mean( N * exp(eta) )
  return( sp_tot*p0 )
}

simulates_sp <- data.frame()
npoints <- 10
for( id in mean_pars$idall ){
  print(id)
  sname <-  (mean_pars %>% filter(idall == id))$sname[1] 
  nmin <- min(as.integer( (proj %>% filter(idall == id))$nreads ))
  nmax <- max(as.integer( (proj %>% filter(idall == id))$nreads ))
  dn = 1.*(log10(nmax)-log10(nmin))/npoints
  sobs <- ( proj %>% filter(idall == id) %>% summarise( sobs = n_distinct(otu_id) ) )$sobs
  beta <- (mean_pars %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars %>% filter(idall == id))$stot[1]
  for(i in -2:(npoints+2) ){
    n = nmin * 10^( i*dn)
    pred_sp <- predicted_sp(n, mu, sigma, beta, stot )
    pred_nr <- predicted_reads(n, mu, sigma, beta, stot )
    pred_shannon <- predicted_shannonindex(n, mu, sigma, beta, stot  )
    d <- as.data.frame( list( idall = id, sname = sname, pred_sp = pred_sp, pred_si = pred_shannon, tot_sp = stot, nreads = pred_nr, sigma = sigma, mu = mu, sobs = sobs, beta = beta, n0 = n )  )
    simulates_sp <- rbind(d, simulates_sp)
  }
} 

psp <- proj %>% group_by(project_id, sample_id, run_id, sname, scat, nreads) %>% summarise( nsp = n_distinct(otu_id) ) %>% ungroup() %>% 
  ggplot() + mytheme +
  aes(
    x = nreads,
    y = nsp,
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point(size = 1,alpha = 1, color = "gray") +
  geom_line( data = simulates_sp, aes( x = n0, y = pred_sp ), color = "black", size = 1.5 ) +
  geom_hline( data = simulates_sp, aes( yintercept = sobs ) , size = 1, linetype = 2) +
  geom_hline( data = mean_pars, aes( yintercept = stot ) , size = 1, linetype = 4) +
  stat_summary_bin(fun.y = "mean",  size = 2.5, stroke = 1.5, geom="point", bins = 10) + scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none") +
  scale_x_log10( "Number of reads"  ) + scale_y_log10( "Number of observed species"  ) 
psp
#ggsave(psp,filename =    "../SI/pred_numberspecies.pdf", width = 10, height = 8  )


pspm <- proj %>% group_by(project_id, sample_id, run_id, sname, scat, nreads) %>% summarise( nsp = n_distinct(otu_id) ) %>% ungroup() %>% 
  filter( sname %in% maintext_figs ) %>% 
  ggplot() + mytheme_main +
  aes(
    x = nreads,
    y = nsp,
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point(size = 1,alpha = 1, color = "gray") +
  geom_hline( data = simulates_sp %>%  filter( sname %in% maintext_figs ), aes( yintercept = sobs ) , size = 1, linetype = 2) +
#  geom_hline( data = mean_pars %>%  filter( sname %in% maintext_figs ), aes( yintercept = stot ) , size = 1, linetype = 4) +
  stat_summary_bin(fun.y = "mean",  size = 3, stroke = 1.5, geom="point", bins = 10) + scalecols + scaleshapes +
  geom_line( data = simulates_sp  %>%  filter( sname %in% maintext_figs ), aes( x = n0, y = pred_sp ), color = "black", size = 1.5 ) +
  facet_wrap( ~ sname ) + theme(legend.position = "none") +
  scale_x_log10( "Number of reads" , label = fancy_scientific ) + scale_y_log10( "Number of\nobserved species" ,  label = fancy_scientific ) +
  geom_line( data = simulates_sp  %>%  filter( sname %in% maintext_figs ), aes( x = n0, y = pred_sp ), color = "black", size = 1.5 ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
pspm
#ggsave(pspm, filename =    "./num_species.pdf", width = 8.25, height = 2.6  )





psh <- proj %>% group_by(project_id, sample_id, run_id, sname, scat, nreads) %>% summarise( si = -sum( 1.*count*log(1.*count/nreads)/nreads ) ) %>% ungroup() %>% 
  ggplot() + mytheme +
  aes(
    x = nreads,
    y = si,
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point(size = 1,alpha = 1, color = "gray") +
  geom_line( data = simulates_sp %>% filter(pred_si < 1000), aes( x = n0, y = pred_si ), color = "black", size = 1.5 ) +
  stat_summary_bin(fun.y = "mean",  size = 2.5, stroke = 1.5, geom="point", bins = 10) + scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none") +
  scale_x_log10( "Number of reads"  ) + scale_y_continuous( "Shannon Index"  ) 
#ggsave(psh,filename =    "../SI/pred_shannon.pdf", width = 10, height = 8  )
psh


### PREDICTION OCCURRENCE DIST
nbin <- 10
occ_binned <- proj %>% dplyr::select(project_id, idall, sname, otu_id, o) %>% 
  distinct() %>% 
  mutate( b = as.integer(o*nbin)  ) %>% 
  group_by( project_id, idall, sname, b ) %>% 
  summarise( do = 1./nbin, n = n(), o = mean(o) ) %>% 
  ungroup() %>% group_by(project_id, idall, sname) %>% 
  mutate( p = n/sum(n)/do ) %>% ungroup()

predict_occupancys <- function(N, eta ){
  p0 = ( (beta/(beta+N*exp(eta))  )^beta )
  return( 1-p0 )
}




expected_occ <- data.frame()
npoints <- 10
for( id in mean_pars$idall ){
  print(id)
  sname <-  (mean_pars %>% filter(idall == id))$sname[1] 
  nreads <- as.integer( (proj %>% filter(idall == id) %>% dplyr::select(nreads, run_id) %>% distinct() )$nreads ) 
  beta <- (mean_pars %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars %>% filter(idall == id))$stot[1] 
  etas <-  rnorm( stot, mean = mu, sd = sigma )
  occupancys <- outer(nreads, etas, FUN=predict_occupancys)
  weight_eta <- 1.-exp( colSums(log(1.- occupancys )) )
  occ_eta <- colMeans(occupancys )
  rnd_occupancys <-  1.*( matrix( runif( dim(occupancys)[1]*dim(occupancys)[2] ) , dim(occupancys)[1], dim(occupancys)[2]) < occupancys )
  rndocc_eta <- colMeans( rnd_occupancys )
  d <- as.data.frame( list( idall = id, sname = sname, o = occ_eta, ornd = rndocc_eta, f = exp(etas), pobs = weight_eta, sigma = sigma, mu = mu, beta = beta )  )
  expected_occ <- rbind(d, expected_occ)
} 

nbins = 20
occ_binned_exp <- expected_occ %>%  mutate( b = as.integer(o*nbin)  ) %>% 
  group_by( idall, sname, b ) %>% 
  summarise( do = 1./nbin, n = sum(pobs), o = mean(o*pobs)/mean(pobs) ,  pobsm = mean(pobs)) %>% 
  ungroup() %>% group_by(idall, sname) %>% 
  mutate( p = n/sum(n)/do ) %>% ungroup()
occ_rndbinned_exp <- expected_occ %>% mutate(o = ornd) %>%  filter(o>0) %>%   mutate( b = as.integer(o*nbin)  ) %>% 
  group_by( idall, sname, b ) %>% 
  summarise( do = 1./nbin, n = sum(pobs), o = mean(o*pobs)/mean(pobs) ,  pobsm = mean(pobs)) %>% 
  ungroup() %>% group_by(idall, sname) %>% 
  mutate( p = n/sum(n)/do ) %>% ungroup()


pocc <- occ_binned %>% ggplot() + mytheme +
  scalecols + scaleshapes +
  aes(
    x = o,
    y = p, 
    color = (sname) ,
    shape = as.factor(sname)
  ) + 
  geom_line( data = occ_rndbinned_exp, aes(x = o, y = p) , color = "black", size = 1.5 ) +
  geom_line( data = occ_binned_exp, aes(x = o, y = p) , color = "darkgray", size = 1.5, linetype = "dashed" ) +
  geom_point(size = 3, stroke = 1.5) +
  facet_wrap( ~ sname ) + theme( legend.position = "none") +
  scale_y_log10( "Probability density"  ) +
  scale_x_continuous( "Occupancy"  ) 
pocc
#ggsave(pocc,filename =    "../SI/pred_occupancydist.pdf", width = 10, height = 8  )


poccm <- occ_binned %>% filter( sname %in% maintext_figs ) %>% ggplot() + mytheme_main +
  scalecols + scaleshapes +
  aes(
    x = o,
    y = p, 
    color = (sname) ,
    shape = as.factor(sname)
  ) + 
  geom_point(size = 3, stroke = 1.5) +
  facet_wrap( ~ sname ) + theme( legend.position = "none") +
  geom_line( data = occ_binned_exp %>%  filter( sname %in% maintext_figs ), aes(x = o, y = p) , color = "black", size = 1.5 ) +
  scale_y_log10( "Probability density" , label = fancy_scientific ) +
  scale_x_continuous( "Occupancy" , label = fancy_linear  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
poccm
#ggsave(poccm,filename =    "./occupancydist.pdf", width = 8, height = 2.6  )


#PREDICTION OCCURRENCE FREQ
occ_freq <- proj %>% 
  dplyr::select(project_id, idall, sname, otu_id, o, f) %>% 
  distinct()

nbins = 12
occ_freq_binned <- occ_freq %>%  group_by(idall, sname ) %>% 
  mutate( l = log10(f), dl = (max(l)-min(l))/nbins ) %>% mutate(b = as.integer( (l-min(l) )/dl )  ) %>% 
  ungroup() %>% group_by( idall, sname, b ) %>% 
  summarise( o = mean(o), f = mean(f)  ) %>% 
  ungroup()
expected_occ_binned <- expected_occ %>%  group_by(idall, sname ) %>% 
  mutate( l = log10(f), dl = (max(l)-min(l))/20 ) %>% mutate(b = as.integer( (l-min(l) )/dl )  ) %>% 
  ungroup() %>% group_by( idall, sname, b ) %>% 
  summarise( o = mean(o), f = mean(f) , pobs = mean(pobs) ) %>% 
  ungroup()



pof <- occ_freq %>% ggplot() + mytheme + 
  scalecols + scaleshapes +
  aes(
    x = log10(f),
    y = o,
    color = (sname) ,
    shape = as.factor(sname)
  ) + geom_point( color = "gray", alpha = 0.1, size = 0.5 ) +
  # stat_summary_bin( fun.y = "mean", geom = "point", size = 3,
  #                       stroke = 1.5, bins = 15 ) +
  geom_point( data = occ_freq_binned  , aes(
    x = log10(f), y = o, color = as.factor(sname), shape = as.factor(sname) ), size = 3, stroke = 1.5 ) +
  stat_summary_bin( data = expected_occ %>% filter( pobs > 0.22 ), 
                    aes(x = log10(f), y = o), fun.y = "mean", geom = "line",
                    size = 1, color = "black" ) +
  facet_wrap( ~ sname) + theme( legend.position = "none") +
  scale_y_log10( "Occupancy" ) +
  scale_x_continuous( "Average relative abundance" )

pof
#ggsave(pof,filename =    "../SI/pred_occfreq.pdf", width = 10, height = 8  )



pofm <- occ_freq %>% filter( sname %in% maintext_figs ) %>% ggplot() + mytheme_main + 
  scalecols + scaleshapes +
  aes(
    x = f,
    y = o,
    color = (sname) ,
    shape = as.factor(sname)
  ) + geom_point( color = "gray", alpha = 0.1, size = 0.5 ) +
  geom_point( data = occ_freq_binned  %>% filter( sname %in% maintext_figs ), aes(
    x = f, y = o, color = as.factor(sname), shape = as.factor(sname) ), size = 3, stroke = 1.5 ) +
  geom_line( data = expected_occ_binned  %>% filter( sname %in% maintext_figs ) %>% filter(pobs > 0.2), aes(
    x = f, y = o ), size = 1.5, color = "black" ) +
  # stat_summary_bin( fun.y = "mean", geom = "point", size = 3,
  #                   stroke = 1.5, bins = 15 ) +
#  stat_summary_bin( data = expected_occ  %>% filter( pobs > 0.22 ) %>% filter( sname %in% maintext_figs ), 
#                    aes(x = log10(f), y = o), fun.y = "mean", geom = "line",
#                    size = 1.5, color = "black" ) +
  facet_wrap( ~ sname) + theme( legend.position = "none") +
  scale_y_log10( "Occupancy" , label = fancy_scientific )  +
  scale_x_log10( "Average relative abundance", label = fancy_scientific ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
pofm
#ggsave(pofm,filename =    "./occfreq.pdf", width = 8., height = 2.6 )

#PREDICTION RSA
predicted_cumsad <- function(N, mu, sigma, beta, sp_tot, klist ){
  eta = rnorm( 100000, mean = mu, sd = sigma)   
  p0 =  mean( (beta/(beta+N*exp(eta))  )^beta ) 
  res = data.frame()
  for(k in klist){
    prob = (beta/(beta+N*exp(eta))) 
    sk = mean(1.- pnbinom(k-1, size = beta, prob = prob, log.p = FALSE) )/(1.-p0)
    res <- rbind( 
      as.data.frame(list(
        k = k, sk = sk
      )), res
    )
  }
  return( res )
}

sads_expected <- data.frame()
for( id in mean_pars$idall ){
  print(id)
  sname <-  (mean_pars %>% filter(idall == id))$sname[1] 
  nmin <- min(as.integer( (proj %>% filter(idall == id))$nreads ))
  nmax <- max(as.integer( (proj %>% filter(idall == id))$nreads ))
  dn = 1.*(log10(nmax)-log10(nmin))/npoints
  beta <- (mean_pars %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars %>% filter(idall == id))$stot[1]
  
  
  klist <- as.integer(10^((0:(log10(nmin/2)*5))/5.))
  sp_pred <- predicted_sp(nmin, mu, sigma, beta, stot )
  res <- predicted_cumsad(nmin, mu, sigma, beta, stot, klist ) %>% mutate(
    idall = id, sname = sname, nreads = nmin
  ) %>%  mutate( sp = sp_pred )
  sads_expected <- rbind(res, sads_expected)
  
  klist <- as.integer(10^((0:(log10(nmax/2)*5))/5.))
  sp_pred <- predicted_sp(nmax, mu, sigma, beta, stot )
  res <- predicted_cumsad(nmax, mu, sigma, beta, stot, klist ) %>% mutate(
    idall = id, sname = sname, nreads = nmax
  )%>%  mutate( sp = sp_pred )
  sads_expected <- rbind(res, sads_expected)
  
} 


sads <- proj %>% group_by( nreads, idall, run_id, sname) %>%
  mutate( sp = n_distinct(otu_id)  ) %>% ungroup() %>% 
  group_by( count, nreads, idall, run_id, sname, sp) %>% 
  summarise( sk = n() ) %>% ungroup() %>%  arrange(idall, sname, run_id, -count) %>% 
  group_by(idall, sname, run_id) %>% mutate( pcum = cumsum(sk)/sum(sk)  ) %>% ungroup()

psads <- sads %>%   ggplot() + mytheme + 
  scalecols + scaleshapes +
  aes(
    x = count,
    y = pcum,
    color = (sname) ,
    shape = as.factor(sname)
  ) + geom_line( aes( group = as.factor(run_id) ), alpha = 0.2 ) +
  geom_line( data = sads_expected, 
             aes( x = k, y = sk , group = as.factor(nreads) ), color = "black", size = 1
  ) +
  facet_wrap( ~ sname ) +
  scale_x_log10( "Number of reads k") + scale_y_log10( "Fraction of species with more than k reads" ) +
  theme( legend.position = "none") 
psads
#ggsave(psads,filename =    "../SI/pred_sads.pdf", width = 10, height = 8  )

psadsm <- sads %>% filter( sname %in% maintext_figs ) %>%   ggplot() + mytheme_main + 
  scalecols + scaleshapes +
  aes(
    x = count,
    y = pcum,
    color = (sname) ,
    shape = as.factor(sname)
  ) + geom_line( aes( group = as.factor(run_id) ), alpha = 0.2 ) +
  geom_line( data = sads_expected %>% filter( sname %in% maintext_figs ) %>%  filter(sk > 0.0001), 
             aes( x = k, y = sk , group = as.factor(nreads) ), color = "black", size = 1.5
  ) +
  facet_wrap( ~ sname ) +
  scale_x_log10( "Number of reads k", label = fancy_scientific) +
  scale_y_log10( name =  "Fraction of species\nwith \u2265k reads", label = fancy_scientific ) +
  theme( legend.position = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
psadsm
#ggsave(psadsm,filename =    "./sads.pdf", width = 8.25, height = 2.6 ,  )

