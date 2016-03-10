require(paranomeKsR)

message("USAGE: Rscript path/2/paranomeKsR/exec/fitMixedGaussianModelAndPlot.R plot.pdf")

input.args <- commandArgs(trailingOnly = TRUE)

#' Fit a mixed Gaussian model to the peaks observed in the bar plots:
ks.vals <- chi.wghtd.ks.tbl[which(!is.na(chi.wghtd.ks.tbl$V3)), 3]
mix.mdl.em <- normalmixEM(ks.vals, k = 3)
x.range <- seq(0, 5, by = 0.01)
comp1 <- data.frame(x = x.range, y = (dnorm(x.range, mean = mix.mdl.em$mu[[1]], sd = mix.mdl.em$sigma[[1]]) * 
    mix.mdl.em$lambda[[1]]), stringsAsFactors = FALSE)
comp2 <- data.frame(x = x.range, y = (dnorm(x.range, mean = mix.mdl.em$mu[[2]], sd = mix.mdl.em$sigma[[2]]) * 
    mix.mdl.em$lambda[[2]]), stringsAsFactors = FALSE)
comp3 <- data.frame(x = x.range, y = (dnorm(x.range, mean = mix.mdl.em$mu[[3]], sd = mix.mdl.em$sigma[[3]]) * 
    mix.mdl.em$lambda[[3]]), stringsAsFactors = FALSE)

#' Find the median of the Ks values of the PMEI genes found in C. hirsuta's
#' paralogous fams:
pmei.ks <- median(chi.wghtd.ks.tbl[with(chi.wghtd.ks.tbl, which(V1 %in% chi.pmei.fams)), 
    3])

#' Plot the results:
plot.tbl <- chi.wghtd.ks.tbl[which(!is.na(chi.wghtd.ks.tbl$V3)), ]
x <- ggplot() + geom_histogram(data = plot.tbl, aes(x = plot.tbl$V3, y = ..density.., 
    alpha = 0.2), position = "identity") + geom_line(data = comp1, aes(x = comp1$x, 
    y = comp1$y), colour = "darkblue") + geom_line(data = comp2, aes(x = comp2$x, 
    y = comp2$y), colour = "darkred") + geom_line(data = comp3, aes(x = comp3$x, 
    y = comp3$y), colour = "darkgreen") + xlab("Ks") + geom_segment(aes(x = pmei.ks, 
    xend = pmei.ks, y = 0.62, yend = 0.61), arrow = arrow(type = "closed"), colour = "black", 
    show_guide = FALSE) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(colour = "black"), 
    axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour = "black"), 
    axis.ticks.y = element_line(colour = "black"), axis.line = element_line(colour = "black"), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
    panel.background = element_blank())
ggsave(x, file = input.args[[1]]) 
