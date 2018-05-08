library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

setwd("/home/gmoc/Desktop/pablo/adaptive_project/results/analysis")

data <- as.data.table(read.table(file = "/home/gmoc/Desktop/pablo/adaptive_project/plans/Opt4D/P02_gustr3754511/cbct_6/adapt_raw/vf.dat",
                                    colClasses = "numeric",
                                    header = TRUE))
data <- data %>% mutate(d = sqrt(vx*vx + vy*vy + vz*vz), beamid = as.factor(beamid))

shift_hist <- ggplot(data, aes(d, fill = beamid)) +
                     geom_histogram(binwidth = 0.001) +
                     xlab("Size (cm)")
xy <- ggplot(data, aes(x, y, color = d)) +
        geom_segment(aes(xend = x + vx, yend = y + vy, color = d),
                     arrow = arrow(length = unit(0.1,"cm"))) +
        xlab("X pos (cm)") + ylab("Y pos (cm)")
xz <- ggplot(data, aes(x, z, color = d)) +
        geom_segment(aes(xend = x + vx, yend = z + vz, color = d),
                 arrow = arrow(length = unit(0.1,"cm"))) +
        xlab("X pos (cm)") + ylab("Z pos (cm)")
yz <- ggplot(data, aes(y, z, color = d)) +
        geom_segment(aes(xend = y + vy, yend = z + vz, color = d),
                 arrow = arrow(length = unit(0.1,"cm"))) +
        xlab("Z pos (cm)") + ylab("Y pos (cm)")

grid.arrange(xy, xz, yz, shift_hist, ncol = 2, nrow = 2)
