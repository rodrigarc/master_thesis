setwd("~/Desktop/R_clonotypes")         #set wd
library(dplyr)
library(DEGreport)
library(gridExtra)
library(ggpubr)
library(stringr)
library(doBy)
library(devtools)
library(MASS)
library(RColorBrewer)

dev_mode(TRUE)               # atm you need to use developers mode for ggplot2::annotation_logticks to work properly
install_github("hadley/scales", force = TRUE)       # main branch of development
install_github("hadley/ggplot2", force = TRUE)            # (Optional) install_github("ggplot2", "kohske", "feature/new-guides-with-gtable")
library(scales)
library(ggplot2)                            # load development version of ggplot2 

# Create a new function "get_density" using MASS package (https://slowkow.com/notes/ggplot2-color-by-density/)
# Use this function to calculate the density of dots in each position on the dot plot (create new colum, row 29)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

Clones <- read.csv("Clonotype detection.csv")
Clones$density <- get_density(Clones$sc.size, Clones$NGS.size, n=1000)
Clones$density <- as.numeric(Clones$density)
  
ggplot(Clones, aes(sc.size, NGS.size, fill = density)) +
  geom_point(shape = 21, size = 6, stroke = 0.5) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_x_continuous(trans = log10_trans(), 
                     limits = c(0.8, 30),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = log10_trans(), 
                     limits = c(0.05, 3000),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb", outside = TRUE) + coord_cartesian(clip = "off") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 7))) +
  labs(x = "Clonotype size by single cell", 
       y = "Clonotype size by RepSeq", 
       title = "Clonotype detection")

ggsave("v4 clonotype detection.pdf", height = 5, width = 5)
