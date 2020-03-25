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

load( "./dataestimate_time.Rdata" )



### ADD LAWS
nbin <- 15
statp <- gamma_pars_time %>%  mutate( lf = log(f ) , cutoff = -100  ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup()


simple_par  <- function(x) { - x^2   }  
p1 <- statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>%  left_join(mean_pars_time) %>%  filter(lf > c, n > 10) %>% 
  arrange(xx, sname) %>% ggplot() + mytheme_main +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3 , stroke = 1.5) +
  scalecols + scaleshapes +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  scale_x_continuous( "Rescaled log average\n relative abundance", limits = c(-2.1,2.1) ) +
  scale_y_log10( "Probability density" , labels = fancy_scientific   )  +
  theme( legend.position = "none" ) #+    facet_wrap( ~ sname )
p1
#ggsave( filename =  "./logn_aver_time.pdf", p1, width = 4., height = 3.5  )


plnind <- statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>%  left_join(mean_pars_time) %>% filter(lf > c, n > 10) %>% 
  arrange(xx, sname) %>% ggplot() + mytheme +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3 , stroke = 1.) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average abundance", limits = c(-2.5,2.5) ) +
  scale_y_log10( "Probability density"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname )
plnind
#ggsave( filename =  "../SI/logn_ind_time.pdf", plnind, width = 10, height = 8  )


#rm(statp)
#TAYLOR
ptaylor_ind <- gamma_pars_time %>%  ggplot() + mytheme +
  aes(
    x = f,
    y = vf
  ) +
  geom_abline( slope = 2, size =2  ) +
  geom_point( size = 0.7 , alpha = 0.5, color = "darkgray") +
  scalecols + scaleshapes +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 3, stroke = 1.5, bins = 13 ) +
  scale_x_log10( "Average abundance", limits = c(-2.5,2.5) ) +
  scale_y_log10( "Variance of abundance across samples"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname ) 
#ggsave( filename =  "../SI/taylor_ind_time.pdf", ptaylor_ind, width = 10, height = 8  )
ptaylor_ind

ptaylor<- gamma_pars_time %>%  ggplot() + mytheme_main +
  aes(
    x = f,
    y = vf
  ) +
  geom_abline( slope = 2., size =2  ) +
  scalecols + scaleshapes +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 3, stroke = 1.5, bins = 13 ) +
  scale_x_log10( "Average abundance" , label = fancy_scientific ) +
  scale_y_log10( "Variance of abundance\nacross samples" , label = fancy_scientific  )  +
  theme( legend.position = "none" ) 
#ggsave( filename =  "./taylor_time.pdf", ptaylor, width = 4.25, height = 3.25  )
ptaylor

## GAMMA
dh <- proj_time %>% filter( o > 0.999 ) %>%
  mutate( l = log(count/nreads), f = count/nreads ) %>% group_by( otu_id, idall, sname ) %>% 
  mutate( ml = mean(l), sl = sd(l), k = mean(f)^2/var(f) ) %>% ungroup() %>% 
  mutate( lf = (l-ml)/sl ) %>% 
  group_by( idall , sname ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df, sname ) %>% 
  summarise(  lf = mean(lf),  n = n() , k = mean(k)) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup()

gammalog  <- function(x, k) { (1.13*x - 0.9 * exp(x)) + 0.5 }  
gammalog  <- function(x, k = 1.7) { ( k*trigamma(k)*x - exp( sqrt(trigamma(k))*x+ digamma(k)) ) - log(gamma(k)) + k*digamma(k) + log10(exp(1)) }  
lognlog  <- function(x) { - x ^ 2 - log10(2 * pi)/2 + log10(exp(1)) } 
p2 <- dh  %>% ggplot() + mytheme_main +
  stat_function( fun = lognlog, color = "black", size = 1, linetype = "dashed" ) +
  stat_function( fun = gammalog, color = "black", size = 1 ) +
  #geom_abline( slope = -1  ) +
  #geom_abline( slope = 1  ) +
  aes(
    x =  lf/ sqrt(2) ,
    y = 10^( log10(p)  ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  )  + geom_point( size = 3, stroke = 1.5 ) +   stat_function( fun = gammalog, color = "black", size = 1 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log\n relative abundance") +
  scale_y_log10( "Probability density", limits = c(0.001,1.), labels = fancy_scientific )  +
  theme( legend.position = "none" ) #+ facet_wrap(~sname)
p2
#ggsave( filename =  "./gamma_fluct_time.pdf", p2, width = 4, height = 3.5  )


p2pans <- dh  %>% ggplot() + mytheme +
  stat_function( fun = lognlog, color = "black", size = 1, linetype = "dashed" ) +
  stat_function( fun = gammalog, color = "black", size = 1 ) +
  aes(
    x =  lf/ sqrt(2) ,
    y = 10^( log10(p)  ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3, stroke = 1.3 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average abundance") +
  scale_y_log10( "Fraction of Species", limits = c(0.001,1.1), label = fancy_scientific )  +
  theme( legend.position = "none" ) + facet_wrap(~sname)
p2pans
#ggsave( filename =  "../SI/gamma_fluct_ind_time.pdf", p2pans, width = 10, height = 8  )
rm(dh)



nbin <- 20

tot_reads_time <- proj_time %>% select(idall, run_id, nreads) %>% distinct() %>% 
  group_by(idall) %>% summarise( treads = sum(nreads)  ) %>% ungroup()

statp_time <- proj_time %>% select(tf,otu_id, o, idall, sname, tvpf) %>% distinct() %>%
  mutate( f=o*tf, vf = o*tvpf ) %>% 
  mutate( f = o*tf ) %>% select(-tf, -o,-tvpf) %>%  left_join(tot_reads_time) %>% 
  mutate( lf = log10(f / sqrt(1+vf/f)) , cutoff = -log10(treads) ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup() %>% 
  
  
  estimate_mean_time <- proj_time %>% select(tf,otu_id, o, idall, sname, tvpf) %>% distinct() %>% 
  mutate( f=o*tf, vf = o*tvpf ) %>% 
  mutate( f = o*tf ) %>% select(-tf, -o) %>%  left_join(tot_reads_time) %>% 
  mutate( lf = log10(f / sqrt(1+vf/f)) , cutoff = -log10(treads*2)  ) %>%  group_by(idall, sname) %>% filter( lf > cutoff ) %>% 
  summarise( c = mean(cutoff),  m1 = mean(lf) , m2 = mean(lf^2)  ) %>%  ungroup() %>% rowwise() %>% 
  mutate( mu = estimatemean_f(c,m1,m2) ) %>% mutate( sigma = sqrt(-c*m1 + m2 + c*mu - m1*mu ) ) %>% 
  ungroup() %>%  as.data.frame()

######## OCCURRENCE PREDICTION

pocc_gamma <- proj_time %>% left_join(gamma_pars_time) %>% left_join(mean_pars_time ) %>% group_by(idall, sname, otu_id) %>% 
  summarize( pocc = 1- mean((1.+theta*nreads)^(-beta ) ), o = mean(o)  ) %>% ungroup() %>% 
  ggplot() + mytheme +
  aes(
    x = pocc,
    y = o
  ) + geom_point(  color = "darkgray", alpha = 0.7 , size = 2, aes( shape = as.factor(sname) ) ) +
  geom_abline( slope = 1, size = 2, color = "black" ) +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 4, stroke = 2, bins = 13 ) +
  scale_y_continuous( "Predicted occurrence" ) +
  scale_x_continuous( "Observed occurrence" ) +   scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none"  )
pocc_gamma
#ggsave( filename =  "../SI/predocc_gamma_time.pdf", pocc_gamma, width = 10, height = 8  )


#### OTHER PRED
predicted_sp <- function(N, mu, sigma, beta, sp_tot ){
  bb = beta
  eta = rnorm( 10^6, mean = mu, sd = sigma )
  p0 = mean( (1.+N*exp(eta)/bb)^-bb )
  return( sp_tot*(1-p0) )
}

predicted_reads <- function(N, mu, sigma, beta, sp_tot ){
  eta = rnorm( 10^6, mean = mu, sd = sigma  )
  p0 = mean( N * exp(eta) )
  return( sp_tot*p0 )
}



simulates_sp <- data.frame()
npoints <- 10
for( id in mean_pars_time$idall ){
  print(id)
  sname <-  (mean_pars_time %>% filter(idall == id))$sname[1] 
  nmin <- min(as.integer( (proj_time %>% filter(idall == id))$nreads ))
  nmax <- max(as.integer( (proj_time %>% filter(idall == id))$nreads ))
  sobs <- ( proj_time %>% filter(idall == id) %>% summarise( sobs = n_distinct(otu_id) ) )$sobs
  dn = max(1.*(log10(nmax)-log10(nmin))/npoints, 0.1)
  beta <- (mean_pars_time %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars_time %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars_time %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars_time %>% filter(idall == id))$stot[1]
  for(i in -2:(npoints+2) ){
    n = nmin * 10^( i*dn)
    pred_sp <- predicted_sp(n, mu, sigma, beta, stot )
    pred_nr <- predicted_reads(n, mu, sigma, beta, stot )
    d <- as.data.frame( list( idall = id, sname = sname, pred_sp = pred_sp, tot_sp = stot, nreads = pred_nr, sigma = sigma, mu = mu, beta = beta, n0 = n, sobs = sobs )  )
    simulates_sp <- rbind(d, simulates_sp)
  }
} 

psp <- proj_time %>% group_by(project_id, sample_id, run_id, sname,  nreads) %>% summarise( nsp = n_distinct(otu_id) ) %>% ungroup() %>% 
  ggplot() + mytheme +
  aes(
    x = nreads,
    y = nsp,
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point(size = 1,alpha = 1, color = "gray") +
  geom_line( data = simulates_sp, aes( x = n0, y = pred_sp ), color = "black", size = 1.5 ) +
  geom_hline( data = simulates_sp, aes( yintercept = sobs ) , size = 1, linetype = 2) +
  geom_hline( data = mean_pars_time, aes( yintercept = stot ) , size = 1, linetype = 4) +
  stat_summary_bin(fun.y = "mean",  size = 2.5, stroke = 1.5, geom="point", bins = 10) + scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none") +
  scale_x_log10( "Number of reads" , label = fancy_scientific ) +
  scale_y_log10( "Number of observed OTUs" , label = fancy_scientific  ) 
psp
#ggsave(psp,filename =    "../SI/pred_numberspecies_time.pdf", width = 10, height = 8  )


predict_occurrences <- function(N, eta ){
  p0 = ( (beta/(beta+N*exp(eta))  )^beta )
  return( 1-p0 )
}

nbin <- 7
occ_binned <- proj_time %>% select(project_id, idall, sname, otu_id, o) %>% 
  distinct() %>% 
  mutate( b = as.integer(o*nbin)  ) %>% 
  group_by( project_id, idall, sname, b ) %>% 
  summarise( do = 1./nbin, n = n(), o = mean(o) ) %>% 
  ungroup() %>% group_by(project_id, idall, sname) %>% 
  mutate( p = n/sum(n)/do ) %>% ungroup()


expected_occ <- data.frame()
npoints <- 10
for( id in mean_pars_time$idall ){
  print(id)
  sname <-  (mean_pars_time %>% filter(idall == id))$sname[1] 
  nreads <- as.integer( (proj_time %>% filter(idall == id) %>% select(nreads, run_id) %>% distinct() )$nreads ) 
  beta <- (mean_pars_time %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars_time %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars_time %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars_time %>% filter(idall == id))$stot[1] 
  etas <-  rnorm( stot, mean = mu, sd = sigma )
  occurrences <- outer(nreads, etas, FUN=predict_occurrences)
  weight_eta <- 1.-exp( colSums(log(1.- occurrences )) )
  occ_eta <- colMeans(occurrences )
  rnd_occurrences <-  1.*( matrix( runif( dim(occurrences)[1]*dim(occurrences)[2] ) , dim(occurrences)[1], dim(occurrences)[2]) < occurrences )
  rndocc_eta <- colMeans( rnd_occurrences )
  d <- as.data.frame( list( idall = id, sname = sname, o = occ_eta, ornd = rndocc_eta, f = exp(etas), pobs = weight_eta, sigma = sigma, mu = mu, beta = beta )  )
  expected_occ <- rbind(d, expected_occ)
} 

nbins = 15
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
  scale_x_continuous( "Occurrence"  ) 
pocc
#ggsave(pocc,filename =    "../SI/pred_occurrencedist_time.pdf", width = 10, height = 8  )



#PREDICTION OCCURRENCE FREQ
occ_freq <- proj_time %>% 
  select(project_id, idall, sname, otu_id, o, f) %>% 
  distinct()

pof <- occ_freq %>% ggplot() + mytheme + 
  scalecols + scaleshapes +
  aes(
    x = log10(f),
    y = o,
    color = (sname) ,
    shape = as.factor(sname)
  ) + geom_point( color = "gray", alpha = 0.1, size = 0.5 ) +
  stat_summary_bin( fun.y = "mean", geom = "point", size = 3,
                    stroke = 1.5, bins = 15 ) +
  stat_summary_bin( data = expected_occ %>% filter( pobs > 0.22 ), 
                    aes(x = log10(f), y = o), fun.y = "mean", geom = "line",
                    size = 1, color = "black" ) +
  facet_wrap( ~ sname) + theme( legend.position = "none") +
  scale_y_log10( "Occurrence" ) +
  scale_x_continuous( "Average relative abundance" )

pof
#ggsave(pof,filename =    "../SI/pred_occfreq_time.pdf", width = 10, height = 8  )



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
for( id in mean_pars_time$idall ){
  print(id)
  sname <-  (mean_pars_time %>% filter(idall == id))$sname[1] 
  nmin <- min(as.integer( (proj_time %>% filter(idall == id))$nreads ))
  nmax <- max(as.integer( (proj_time %>% filter(idall == id))$nreads ))
  dn = 1.*(log10(nmax)-log10(nmin))/npoints
  beta <- (mean_pars_time %>% filter(idall == id))$mbeta[1]
  mu <- (mean_pars_time %>% filter(idall == id))$mu[1]
  sigma <- (mean_pars_time %>% filter(idall == id))$sigma[1]
  stot <- (mean_pars_time %>% filter(idall == id))$stot[1]
  
  
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


sads <- proj_time %>% group_by( nreads, idall, run_id, sname) %>%
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
  scale_x_log10( "Number of reads k", label = fancy_scientific) + scale_y_log10( "Fraction of species with more than k reads", label = fancy_scientific ) +
  theme( legend.position = "none") 
psads
#ggsave(psads,filename =    "../SI/pred_sads_time.pdf", width = 10, height = 8  )




#################### TIME CORRELATION

occurrence_threshold <- 0.5
proj_ranks <- proj_time %>%  filter(o > occurrence_threshold) %>%  mutate(ft = count) %>% select(otu_id, ft, idall, experiment_day, sname, nreads) %>% spread( otu_id, ft, fill = 0 ) %>% 
  gather(otu_id, "ft", -idall,-experiment_day, -sname, -nreads) %>% group_by(idall, sname, otu_id) %>% 
  mutate( ctot = sum(ft) ) %>% ungroup() %>%  filter(ctot > 0) %>% 
  mutate(ft = (ft+0.00001)/nreads ) %>% group_by(idall, sname, otu_id) %>% 
  mutate( rf = rank(ft)/n() ) %>%  ungroup() %>%  arrange(otu_id)

s <- proj_ranks %>%   arrange( idall, sname, otu_id, experiment_day ) %>% 
  group_by(idall, sname, otu_id) %>%
  mutate(exp2 = lead(experiment_day), rf2 = lead(rf) ) %>% 
  ungroup(  ) %>% filter(!is.na(rf2)) %>% 
  mutate(  dt= exp2-experiment_day ) %>% 
  filter( dt == 1) %>% select(-dt, -exp2)


binw  <- 0.1
smeans <- s %>% 
  group_by(idall, sname,otu_id) %>%
  mutate( x = rf - 0.5 , y = rf2 -0.5  ) %>% #### <- 
  #mutate( x = log10(f/theta/(k) ) , y = log10(f2/theta)  ) %>% #### <- 
  mutate( b = as.integer( x / binw ) ) %>%  ungroup() %>% 
  group_by(idall, sname, b, otu_id) %>% 
  mutate( eta = y - mean(y), mfy = mean(y), mfx = mean(x),   nd = n() )  %>% ungroup() %>% 
  filter(nd > 0) %>% 
  group_by(idall, sname,  b ) %>% 
  summarise( m_x =  mean( mfx ), m_y = mean( mfy ) , 
             err_y = sd( eta )/sqrt(n()),
             n_y = sd( eta ), err_ny = sd( (eta - mean(eta) )^2  )/sqrt(n()),
             n = n()  ) %>%  ungroup() 



smeans_rnd <- proj_ranks %>% select(idall, experiment_day) %>% distinct() %>% 
  mutate( exp_sh = sample(experiment_day)  ) %>% right_join(proj_ranks) %>% 
  mutate( experiment_day = exp_sh ) %>%  select(-exp_sh) %>% arrange( idall, sname, otu_id, experiment_day ) %>% 
  group_by(idall, sname, otu_id) %>%
  mutate(exp2 = lead(experiment_day), rf2 = lead(rf) ) %>% 
  ungroup(  ) %>% filter(!is.na(rf2)) %>% 
  mutate(  dt= exp2-experiment_day ) %>% 
  filter( dt == 1) %>% select(-dt, -exp2) %>% 
  group_by(idall, sname,otu_id) %>%
  mutate( x = rf - 0.5 , y = rf2 -0.5  ) %>% 
  mutate( b = as.integer( x / binw ) ) %>%  ungroup() %>% 
  group_by(idall, sname, b, otu_id) %>% 
  mutate( eta = y - mean(y), mfy = mean(y), mfx = mean(x),   nd = n() )  %>% ungroup() %>% 
  filter(nd > 0) %>% 
  group_by(idall, sname,  b ) %>% 
  summarise( m_x =  mean( mfx ), m_y = mean( mfy ) , 
             err_y = sd( eta )/sqrt(n()),
             n_y = sd( eta ), err_ny = sd( (eta - mean(eta) )^2  )/sqrt(n()),
             n = n()  ) %>%  ungroup() 

mainsubset <- c( "feces M3", "L_palm M3", "Tongue M3"  )
p1 <- smeans %>%  filter( sname %in% mainsubset ) %>% ggplot() + mytheme_main +
  geom_abline( size = 1.5, slope = 0, color = "#696969" ) +
  aes( y = m_y, ymin = m_y - err_y, ymax = m_y + err_y,  x = m_x,
       color = (sname) ,
       shape = as.factor(sname)
       ) + geom_point( size = 3, stroke = 1.5) +
  scalecols + scaleshapes +
  facet_wrap( ~ sname) +
#  geom_point(  aes( y = m_y , x = m_x ), data = smeans_rnd , color = "red" ) +
 # geom_abline( slope = 0.3, size = 1 ) +
  theme( legend.position = "none"  ) +
  scale_x_continuous( name = "Abundance (quantile) at time t" ) +
  scale_y_continuous( name = "Average abundance\n(quantile)\nat time t+1 [days]", limits = c(-0.25,0.25)  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p1
#ggsave(p1, filename =    "./average_fluct.pdf", width = 8.25*1.15, height = 2.6*1.15  )

p2 <- smeans %>%  filter( sname %in% mainsubset ) %>% ggplot() + mytheme_main +
  geom_hline( size = 1.5, yintercept = 1/sqrt(12), color = "#696969" ) +
  aes( y = n_y, ymin = m_y - err_y, ymax = m_y + err_y,  x = m_x,
       color = (sname) ,
       shape = as.factor(sname)
  ) + geom_point( size = 3, stroke = 1.5) +
  scalecols + scaleshapes +
  facet_wrap( ~ sname) +
  #  geom_point(  aes( y = n_y , x = m_x ), data = smeans_rnd , color = "red" ) +
  # geom_abline( slope = 0.3, size = 1 ) +
  theme( legend.position = "none"  ) +
  scale_x_continuous( name = "Abundance (quantile) at time t" ) +
  scale_y_continuous( name = "Variance of\n(quantile) abundance\nat time t+1 [days]", limits = c(0.2, 0.31), breaks = c(0.2,0.25, 0.3) ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p2
#ggsave(p2, filename =    "./variance_fluct.pdf", width = 8.25*1.15, height = 2.6*1.15  )

