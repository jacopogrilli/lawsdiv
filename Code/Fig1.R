rm(list=ls(all=TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(tidyverse)
require(ggfortify)
require(lamW)
require(GGally)
require(ggcorrplot)
require(pracma)
require(nleqslv)
library(grid)
library(gridExtra) 
source("./ggplot_mytheme.R")

load( "./dataestimate.Rdata" )

##### LOGNORMAL
nbin <- 20
statp <- gamma_pars %>%  mutate( lf = log(f ) , cutoff = -100  ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup()


simple_par  <- function(x) { - x^2   }  
p1 <- statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>%  left_join(mean_pars) %>% filter(lf > c+0.15, n > 10) %>% 
  arrange(xx, sname) %>% ggplot() + mytheme_main +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3 , stroke = 1.5) +   stat_function( fun = simple_par, color = "black", size = 1 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average\n relative abundance", limits = c(-2.5,2.5) ) +
  scale_y_log10( "Probability density" , labels = fancy_scientific   )  +
  theme( legend.position = "none" ) #+    facet_wrap( ~ sname )
p1
#ggsave( filename =  "./logn_aver.pdf", p1, width = 4., height = 3.5  )


plnind <- statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>%  left_join(mean_pars) %>% filter(lf > c, n > 10) %>% 
  arrange(xx, sname) %>% ggplot() + mytheme +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 4 , stroke = 1.5) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average relative abundance", limits = c(-2.5,2.5) ) +
  scale_y_log10( "Probability density"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname )
plnind
#ggsave( filename =  "../SI/logn_ind.pdf", plnind, width = 10, height = 8  )


rm(statp)


####### TAYLOR
nbins <- 15
taylor_binnes <- gamma_pars %>% group_by( idall, sname ) %>% mutate( lf = log(f), dlf = (max(lf)-min(lf))/nbins  ) %>% 
  mutate(b = as.integer( (lf-min(lf))/dlf )  ) %>% ungroup() %>% group_by(idall, sname,b) %>% 
  summarise( vf = mean(vf), f = mean(f)  ) %>% ungroup() %>% mutate(xx = runif(n() ) ) %>% arrange(xx) 

ptaylor_ind <- gamma_pars %>%  ggplot() + mytheme +
  aes(
    x = f,
    y = vf
  ) +
  geom_abline( slope = 2, size =2  ) +
  #geom_abline( slope = 1.5, size =2, intercept = -0.5  ) +
  geom_point( size = 0.5 , alpha = 0.3, color = "gray") +
  scalecols + scaleshapes +
  geom_point( data = taylor_binnes,  aes(x = f, y = vf, color = as.factor(sname), shape = as.factor(sname) ), size = 3, stroke = 1.5, bins = 13 ) +
  scale_x_log10( "Average abundance", limits = c(-2.5,2.5) ) +
  scale_y_log10( "Variance of abundance across samples"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname ) 
#ggsave( filename =  "../SI/taylor_ind.pdf", ptaylor_ind, width = 10, height = 8  )
ptaylor_ind

ptaylor_all <- taylor_binnes  %>% 
  ggplot() + mytheme_main +
  aes(
    x = f,
    y = vf,
    color = as.factor(sname),
    shape = as.factor(sname)
  ) +
  scalecols + scaleshapes +
  geom_point(size = 3, stroke = 1.5 ) +   geom_abline( slope = 2, size = 2, intercept = 1.  ) +
  scale_x_log10( "Average\nrelative abundance", limits = c(-2.5,2.5), labels = fancy_scientificb ) +
  scale_y_log10( "Variance of\nrelative abundance" , labels = fancy_scientific  )  +
  theme( legend.position = "none" ) 
#ggsave( filename =  "./taylor_ind.pdf", ptaylor_all, width = 4.3, height = 3.5  )
ptaylor_all


######### GAMMA
nbin <- 20
dh <- proj %>% filter( o > 0.999 ) %>%
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
p2 <- dh %>%  mutate( xx = runif(dim(dh)[1])  ) %>%   arrange(xx, sname) %>% ggplot() + mytheme_main +
  #stat_function( fun = lognlog, color = "black", size = 1, linetype = "dashed" ) +
  aes(
    x =  lf/ sqrt(2) ,
    y = 10^( log10(p)  ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3, stroke = 1.5 ) +   stat_function( fun = gammalog, color = "black", size = 1 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log\n relative abundance") +
  scale_y_log10( "Probability density", limits = c(0.001,0.8), labels = fancy_scientific )  +
  theme( legend.position = "none" ) #+ facet_wrap(~sname)
p2
#ggsave( filename =  "./gamma_fluct.pdf", p2, width = 4, height = 3.5  )


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
  scale_y_log10( "Fraction of Species", limits = c(0.001,1.1) )  +
  theme( legend.position = "none" ) + facet_wrap(~sname)
#ggsave( filename =  "../SI/gamma_fluct_ind.pdf", p2pans, width = 10, height = 8  )
p2pans
rm(dh)

nbin <- 10
dhind <- proj %>% filter( o > 0.999 ) %>%
  mutate( l = log(count/nreads) ) %>% group_by( otu_id, idall, sname ) %>% 
  mutate( ml = mean(l), sl = sd(l) ) %>% ungroup() %>% 
  mutate( lf = (l-ml)/sl ) %>% 
  group_by( otu_id, idall , sname ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df, sname, otu_id ) %>% 
  summarise(  lf = mean(lf),  n = n() ) %>% 
  ungroup() %>% group_by(idall, otu_id) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup()

p2_ind <- dhind %>%  arrange(lf) %>% ggplot() + mytheme +
  stat_function( fun = lognlog, color = "black", size = 1, linetype = "dashed" ) +
  aes(
    x =  lf/ sqrt(2) ,
    y = 10^( log10(p)  ),
    group = as.factor(otu_id),
  ) + geom_line(   alpha = 0.7  , color = "gray"  ) +
  stat_function( fun = gammalog, color = "black", size = 1 , args = list(k = 1.7)) +
  scalecols + scaleshapes +
  geom_point( data = dh %>% mutate(otu_id = NA),  aes(
    x =  lf/ sqrt(2) ,
    y = 10^( log10(p)  ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ),  size = 2, stroke = 1. ) +
  scale_x_continuous( "Rescaled log average abundance") +
  scale_y_log10( "Fraction of Species", limits = c(0.001,1) )  +
  theme( legend.position = "none" ) + facet_wrap(~sname)
p2_ind
#ggsave( filename =  "../SI/gamma_fluct_ind.pdf", p2_ind, width = 10, height = 8  )
rm(dhind)
