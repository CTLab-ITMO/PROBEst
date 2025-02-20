library(dplyr)
library(ggplot2)

pwd <- commandArgs(trailingOnly = TRUE)

df <- read.csv(pwd)

outwd <- paste0(dirname(pwd), "/", "plots")
dir.create(outwd)

changed <- df %>% 
  apply(2, function(x) unique(x) %>% length) > 1
changed[c(2:4, ncol(df))] <- F


for(var in names(df)[changed]) {
  gg <- ggplot(df, aes_string(x = "iteration", y = "mean_hits", color = var, group = var, fill = var)) +
    geom_point() +
    geom_smooth(alpha = .3, formula = 'y~x', method = "loess") +
    theme_minimal() +
    theme(legend.position = "bottom")
    
  ggsave(paste0(outwd,"/", var, ".svg"), gg, width = 7, height = 4)
}