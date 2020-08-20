rm(list=ls(all=TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(tidyverse)
require(ggfortify)
require(lamW)
require(GGally)
require(ggcorrplot)
require(pracma)
source("./ggplot_mytheme.R")

remove_runs = c("ERR1104477", "ERR1101508", "SRR2240575") ## duplicate runs for the same sample


### LOAD STATIC DATA
load( "../Data/crosssecdata.RData" )
nreads_cutoff <- 10^1
count_cutoff <- 0


summarydata <- datatax %>% filter(nreads > nreads_cutoff ) %>%  group_by(project_id, classification) %>%
  summarise( n_of_runs = n_distinct(run_id), mean_nreads = mean(nreads)) %>%  ungroup() %>% 
  as.data.frame() %>%  filter( n_of_runs > 00  ) %>% mutate(idall = paste(project_id, classification))
datatax <- datatax %>% filter( paste(project_id, classification) %in% paste(summarydata$project_id,summarydata$classification  )  )
summarydata
names <- c("sludge", "lake", "river", "gut2", "sea vents", "glacier", "seawater", "soil", "gut1", "oral1")
scat <- c( "ref", "water", "water", "host", "water", "water", "water", "soil", "host", "host"  )
shortnames <- as.data.frame(list( idall = summarydata$idall, sname = names, scat = scat  ))
proj <- datatax %>% filter(  nreads > nreads_cutoff, count > count_cutoff ) %>% filter( ! run_id %in% remove_runs ) %>% 
  group_by( project_id, classification, otu_id ) %>% mutate( tf = mean(count/nreads), o = n(), tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% ungroup() %>%
  group_by( project_id, classification ) %>%  mutate(o = o / n_distinct(run_id) ) %>% 
  mutate( f = o*tf, vf = o*tvpf ) %>% mutate(vf = vf - f^2 ) %>%  ungroup() %>% 
  mutate(idall = paste(project_id, classification)) %>% select( -tf, -tvpf )
proj <- proj %>% left_join( shortnames)  %>% filter( sname != "sea vents" )


############ PARAMETERS OF THE GAMMA DISTRIBUTION
gamma_pars <- proj %>% select( idall, sname, otu_id, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta )


gamma_pars %>%  ggplot() + mytheme +
  aes(
    x = f,
    y = vf
  ) + geom_point( alpha = 0.5 ) + facet_wrap(  ~ sname ) +
  scale_x_log10() + scale_y_log10()


gamma_pars %>% filter(f > 2*10^-5) %>% group_by(sname) %>% summarise( cor(beta, log(f)) ) %>%  ungroup()
### think about this correlation (sampling) and variability of cv/beta (could be able to say something from the model)


############## LOGNORMAL PARAMETERS
nbin <- 20
cutoffs <- list( sname = c("sludge", "lake", "river", "gut2", "glacier", "seawater", "soil", "gut1", "oral1" ),
                 c = c(-18.8, -17, -15., -17.50000, -16.25, -17.9, -14.26616, -14.29533, -16.1) ) %>% as.data.frame()

fun_erf <- function( mu, c, m1, m2  ){
  sigma <- sqrt(-c*m1 + m2 + c*mu - m1*mu )
  x <- (c-mu)/sigma/sqrt(2.)
  f <- (mu-m1)*erfc(x) + exp(- x^2) * sqrt(2/pi) * sigma 
  return(f)
}

estimatemean_f <- function( c, m1, m2 ){
  mumin <- (c*m1 - m2)/(c - m1)
  muest <- uniroot(fun_erf, c = c, m1 = m1, m2 = m2, interval = c(-50,mumin-0.001), tol = 0.0001 )$root
  return( muest )
}



estimate_mean <- gamma_pars %>% 
  mutate(lf = log(f)) %>%
  left_join(cutoffs) %>%  ungroup() %>% filter( lf > c ) %>% 
  group_by(idall, sname) %>% 
  summarise( c = mean(c),  m1 = mean(lf) , m2 = mean(lf^2), ns_obs = n_distinct(otu_id), nf = sum(f)  )  %>%   ungroup() %>% rowwise() %>% 
  mutate( mu = estimatemean_f(c,m1,m2) ) %>% mutate( sigma = sqrt(-c*m1 + m2 + c*mu - m1*mu ) ) %>% 
  ungroup() %>%  as.data.frame() %>% mutate( stot = 2*ns_obs / erfc( (c-mu)/sigma/sqrt(2)  ) ) # , np_tot_r = 1./exp(mu+sigma^2/2.)/nf )


statp <-  gamma_pars %>% 
  mutate( lf = log(f ) , cutoff = -100  ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup()

simple_par  <- function(x) { - x^2   }  
p1 <- statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>% filter(lf > cutoff) %>% left_join(estimate_mean) %>% 
  arrange(xx, sname) %>% ggplot() + mytheme +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) ,
    shape = as.factor(sname)
  ) + geom_point( size = 3 , stroke = 1) +
  geom_vline( aes( xintercept = (c-mu)/sqrt(2*sigma^2)  )  ) +
  geom_vline( aes( xintercept = 0  )  ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average abundance" ) +
  scale_y_log10( "Fraction of Species"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname )
p1


mean_pars <- estimate_mean %>% select( idall, sname, nf, c, mu, sigma, stot  ) %>% left_join(
  gamma_pars %>% filter(f > 10^-5, vf > 0 ) %>% group_by(sname) %>% summarise( mbeta = mean(beta) ) %>%  ungroup()
)

save( proj, gamma_pars, mean_pars , file = "./dataestimate.Rdata" )



##### REPEAT WITH TIME

nreads_cutoff <- 10^3
count_cutoff <- 0


load( "../Data/longitudinal.RData")

proj_time <- proj_time %>% filter(nreads > nreads_cutoff) %>% 
  group_by( project_id, classification, host_id, otu_id, samplesite ) %>% mutate( tf = mean(count/nreads), o = n(), tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% ungroup() %>%
  group_by( project_id, classification, host_id, samplesite ) %>%  mutate(o = o / n_distinct(run_id) ) %>% ungroup() %>% 
  mutate( f = o*tf, vf = o*tvpf ) %>% mutate(vf = vf - f^2 ) %>%  
  select( -tf, -tvpf ) %>% mutate( idall = paste( project_id, samplesite, host_id ),  sname = paste(samplesite, host_id) ) 
#  left_join(alternative_taxonomy2 %>% select(otu_id, genus) ) %>%  filter( !is.na(genus)) %>%  select(-genus)


gamma_pars_time <- proj_time %>% select( idall, sname, otu_id, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta )

gamma_pars_time %>%  ggplot() + mytheme +
  aes(
    x = f,
    y = vf
  ) + geom_point( alpha = 0.5 ) + facet_wrap(  ~ sname ) +
  scale_x_log10() + scale_y_log10() +
  geom_abline( slope = 2)

nbin <- 10

cutoffs_time <- list( sname = c("feces F4", "feces M3", "L_palm F4", "L_palm M3", "R_palm F4", "R_palm M3", "Tongue F4", "Tongue M3" ),
                 c = c(-16.25, -16.75, -15.5, -15.7, -14.2, -17, -14.75, -17.0) ) %>% as.data.frame()


statp_time <- gamma_pars_time %>%  
  mutate( lf = log(f ) , cutoff = -100 ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  mutate( p = n / sum(n) / df ) %>%  ungroup() 
  
  
estimate_mean_time <- gamma_pars_time %>%  left_join(cutoffs_time) %>% 
  mutate( lf = log(f ) , cutoff = c  ) %>%  group_by(idall, sname) %>% filter( lf > cutoff ) %>% 
  summarise( c = mean(cutoff),  m1 = mean(lf) , m2 = mean(lf^2),  ns_obs = n_distinct(otu_id), nf = sum(f)) %>%  ungroup() %>% rowwise() %>% 
  mutate( mu = estimatemean_f(c,m1,m2) ) %>% mutate( sigma = sqrt(-c*m1 + m2 + c*mu - m1*mu ) ) %>% 
  ungroup() %>%  mutate( stot = 2*ns_obs / erfc( (c-mu)/sigma/sqrt(2)  ) ) %>%   as.data.frame()


simple_par  <- function(x) { - x^2   }  
p1 <- statp_time %>%  mutate( xx = runif(dim(statp_time)[1])  )  %>% filter(lf > cutoff) %>% left_join(estimate_mean_time) %>% 
  arrange(xx, sname)  %>% ggplot() + mytheme +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  aes(
    x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) )
  ) + geom_point( size = 3 , stroke = 1, shape = 1) +
  geom_vline( aes( xintercept = (c-mu)/sqrt(2*sigma^2)  )  ) +
  geom_vline( aes( xintercept = 0  )  ) +
  #scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average abundance" ) +
  scale_y_log10( "Fraction of Species"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname )
p1

mean_pars_time <- estimate_mean_time %>% select( idall, sname, nf, c, mu, sigma, stot  ) %>% left_join(
  gamma_pars_time %>% filter(f > 10^-5, vf > 0 ) %>% group_by(sname) %>% summarise( mbeta = mean(beta) ) %>%  ungroup()
)

save( proj_time, gamma_pars_time, mean_pars_time , file = "./dataestimate_time.Rdata" )




##### Calculate correlation between averages in cross sectional data
require(corrplot)
cor_test_wrapper = function(col_name1, col_name2, data_frame) {
  a <- data_frame[[col_name1]]
  b <- data_frame[[col_name2]]
  print(a)
  cor.test(data_frame[[col_name1]], data_frame[[col_name2]] )$p.value
}

df <- gamma_pars %>% select(sname, otu_id, f) %>% filter(sname != "seawater") %>% mutate(f = log(f)) %>% 
  spread( sname, f, fill = NA ) %>% select(-otu_id) %>% as.matrix()  
corr_df <- cor(df, use = "pairwise.complete.obs" )
corr_p <- cor.mtest(df, use = "pairwise.complete.obs" , conf.level = .99)

pdf("./SI/corr_averages.pdf", height = 8, width = 8)
corrplot.mixed(corr_df, p.mat = corr_p$p, sig.level = 0.001,  order = "hclust", upper = "circle",
               tl.col = "#444444", tl.cex = 1.3, pch = 4, pch.cex = 2, lower.col = "#888888" )
dev.off()

col_combinations = expand.grid(colnames(df), colnames(df))
p_vals = mapply(cor_test_wrapper, 
                col_name1 = col_combinations[[1]], 
                col_name2 = col_combinations[[2]], 
                MoreArgs = list(data_frame = df))
matrix(p_vals, 3, 3, dimnames = list(names(df), names(df)))



