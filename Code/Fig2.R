rm(list=ls(all=TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(tidyverse)
require(ggfortify)
require(lamW)
require(GGally)
require(ggcorrplot)
require(pracma)
source("./ggplot_mytheme.R")

load( "./dataestimate.Rdata" )



################### PREDICT OCCURRENCE
pocc_gamma_ind <- proj %>% left_join(gamma_pars) %>% left_join(mean_pars ) %>%
  mutate( beta = ifelse(is.na(beta), mbeta, beta)  ) %>% mutate( theta = f/beta  ) %>%  
  group_by(idall, sname, otu_id) %>% 
  summarize( pocc = 1- mean((1.+theta*nreads)^(-beta ) ), o = mean(o)  ) %>% ungroup() %>% 
  ggplot() + mytheme +
  aes(
    x = pocc,
    y = o
  ) + geom_point(  color = "gray", alpha = 0.2 , aes( shape = as.factor(sname) ) ) +
  geom_abline( slope = 1, size = 2, color = "black" ) +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 4, stroke = 2, bins = 13 ) +
  scale_y_continuous( "Predicted occupancy" ) +
  scale_x_continuous( "Observed occupancy" ) +   scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none"  )
pocc_gamma_ind
#ggsave( filename =  "./predocc_gamma_ind.pdf", pocc_gamma_ind, width = 10, height = 8  )

pocc_gamma <- proj %>% left_join(gamma_pars) %>% left_join(mean_pars ) %>%
  mutate( beta = ifelse(is.na(beta), mbeta, beta)  ) %>% mutate( theta = f/beta  ) %>%  
  group_by(idall, sname, otu_id) %>% 
  summarize( pocc = 1- mean((1.+theta*nreads)^(-beta ) ), o = mean(o)  ) %>% ungroup() %>% mutate(xx = runif(n() ) ) %>% 
  arrange(xx) %>% 
  ggplot() + mytheme_main +
  aes(
    x = pocc,
    y = o
  ) + #geom_point(  color = "gray", alpha = 0.2 , aes( shape = as.factor(sname) ) ) +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 4, stroke = 2, bins = 13 ) +
  geom_abline( slope = 1, size = 2, color = "black" ) +
  scale_y_continuous( "Predicted occupancy\n(from mean and variance\nof relative abundance)" ) +
  scale_x_continuous( "Observed occupancy\n(fraction of samples\nwhere an OTU is present)" ) +   scalecols + scaleshapes +
  theme(legend.position = "none"  )
pocc_gamma
#ggsave( filename =  "./predocc_gamma.pdf", pocc_gamma, width = 4.5, height = 3.5  )



#### LOGNORMAL DOES NOT PREDICT OCCURRENCE


lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,
                    min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}

logn_occ <- function( x, sigma, mu  ){  Re(   exp(- ( lambertW(-x * sigma^2 * exp(mu) )^2
                                                      + 2 * lambertW( -x * sigma^2 * exp(mu) ) )/ (2*sigma^2) )  / sqrt(1 + lambertW( -x * sigma^2 * exp(mu)  ))  )  }

op <- proj  %>%  filter(vf > 0 ) %>% 
  mutate( mu = log( f / sqrt( 1+vf/f^2 ) ) , s2 = log(1.+vf/f^2)  )  %>% 
  mutate(  oexp = 1. - logn_occ(-1, sqrt(s2), mu + log(nreads)  ) ) %>% 
  group_by( sname, idall, otu_id, o ) %>%  summarise( oexp = mean(oexp)  ) %>%  ungroup()
  
  
pocc_logn <- ggplot( op ) + mytheme +
  aes(x = o, y = oexp , color = as.factor(sname)) +
  geom_point(  color = "gray", alpha = 0.2 , aes( shape = as.factor(sname) ) ) +
  geom_abline( slope = 1, size = 2, color = "black" ) +
  stat_summary_bin(fun.y = "mean", geom = "point", aes( color = as.factor(sname), shape = as.factor(sname) ), size = 4, stroke = 2, bins = 13 ) +
  scale_y_continuous( "Predicted occupancy" ) +
  scale_x_continuous( "Observed occupancy" ) +   scalecols + scaleshapes +
  facet_wrap( ~ sname ) + theme(legend.position = "none"  )
pocc_logn
#ggsave( filename =  "../SI/predocc_logn_ind.pdf", pocc_logn, width = 10, height = 8  )
