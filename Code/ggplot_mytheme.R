require(tidyverse)
require(scales)
require(RColorBrewer)

mytheme <- theme_bw() + theme(
  legend.title  = element_text( size=17),
  #  legend.position = "bottom",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text( size=17),
  panel.background = element_rect(fill=NA),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=19),
  panel.border = element_rect( colour = "black", size=2),
  axis.ticks = element_line(size = 1.),
  legend.background = element_rect(fill=NA)
)

mytheme_main <- theme_bw() + theme(
  legend.title  = element_text(family="Helvetica", size=17, color = "#222222"),
  legend.key = element_blank(),
  legend.text  = element_text(family="Helvetica", size=17, color = "#222222"),
  panel.background = element_rect(fill="transparent"),
  #plot.background = element_rect(fill="transparent", colour = NA),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=17, color = "#222222"),
  panel.border = element_blank(),
  axis.title = element_text( family="Helvetica", size=17, color = "#222222"),
  axis.text = element_text( family="Helvetica", size=15, color = "#222222"),
  axis.line = element_line(size = 1., color = "#222222"),
  axis.ticks = element_line(size = 1.,color = "#222222"),
  legend.background = element_rect(fill="transparent", colour = NA)
)


fancy_linear <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = FALSE)
  # return this as an expression
  parse(text=l)
}



fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # e+00 becomes 1
  l <- gsub("e\\+00", "", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}


fancy_scientificb <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}


scientific_10_exp_labels <- trans_format("log10", math_format(10^.x) )
scientific_10_exp_breaks <- trans_format("log10", function(x) 10^x )

scalecols <- scale_colour_manual(values = c(
  "glacier" = "#00d667", "lake" = "#0000FF",  "river" =  "#9933FF",   "sea vents" = "#0000CC",
  "gut2" =  "#FF0000", "gut1" =  "#FF6600", "seawater" = "#6699FF" , "oral1" = "#FFCC00", 
  "sludge" =  "#993333", "soil" =  "#339933",
  "feces F4" = "#FF3300", "feces M3" = "#FF3300", "L_palm F4" = "#ff967c", "L_palm M3" = "#ff967c", "R_palm F4" = "#ff967c", "R_palm M3" = "#ff967c",
  "Tongue F4" = "#d6ab00", "Tongue M3" = "#d6ab00" 
  ) )


scaleshapes <- scale_shape_manual( values = c( "glacier" = 2, "lake" = 3,  "river" = 5,   "sea vents" = 9,
                                               "gut2" =  0, "gut1" =  7, "seawater" = 6, "oral1" = 8, 
                                               "sludge" =  1, "soil" =  4,  "feces F4" = 7, "feces M3" = 0, "L_palm F4" = 1,
                                               "L_palm M3" = 2, "R_palm F4" = 3, "R_palm M3" = 4,
                                               "Tongue F4" = 6, "Tongue M3" = 8  ))
