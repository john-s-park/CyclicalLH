##############

## Code to: 1) create M matrix from system of ODE's, 2) create periodic P matrix from M, and 3) produce fitness landscapes from P

# NOTE: 
# > This code gives the general machinery that is central to all subsequent analyses. It is tuned with some arbitrary parameters and trade-offs between parameters to demonstrate one possible example. These parameters are noted throughout the code.
# > Output of this code produces Figure 2 of the paper.
# > 'all_optimality_curves.R' is the more extensive code that produces a wide range of model variants. (which produces SI Fig. S3)

##############

# packages
library(reshape2)
library(ggplot2)
library(viridis)

# range of life history parameters
mu_seq <- seq(0.01, 0.05, length.out = 100) 
d_seq <- seq(0.10, 0.40, length.out = 100)
f_seq <- seq(4.00, 17.00, length.out = 100)
gamma_seq <- seq(0.01, 0.20, length.out = 100)

# range of cycle periodicity (T) across which to eventually create fitness landscapes
t_seq<- seq(0.5,14,length.out = 100)

# juvenile- and adult-specific mortality associated with cyclical disturbance events (modelled after demographic disturbance data from sampled populations)
Sj = 0.59
Sa = 0.94

# M as a function of the four life history parameters from the system of ODE's (continuous juv and adult abundances)
# named 'fulltradeoff' because this particular example will eventually include three-way trade-offs (producing Fig 2, and SI, Fig. S2 "M")
M.fulltradeoff <- function(mu, d, f, gamma){
  M <- matrix(c(-(mu+d),
                mu, 
                f, 
                -gamma), 
              2,2)
}

# P as a function of: 1) eigenvalues and eigenvector elements of M, 2) T, and 3) Sj & Sa
P.calc <-function(M,t,Sj,Sa){
  lambdas <- eigen(M)$values
  vecs <- eigen(M)$vectors
  P <- matrix(c(Sj*((vecs[2,1]*exp(lambdas[2]*t)*vecs[1,2])-(vecs[2,2]*exp(lambdas[1]*t)*vecs[1,1]))/
                  (vecs[1,2]*vecs[2,1]-vecs[2,2]*vecs[1,1]), #P[1,1]
                Sa*((vecs[2,1]*exp(lambdas[2]*t)*vecs[2,2])-(vecs[2,2]*exp(lambdas[1]*t)*vecs[2,1]))/
                  (vecs[1,2]*vecs[2,1]-vecs[2,2]*vecs[1,1]), #P[2,1]
                Sj*((vecs[1,2]*exp(lambdas[1]*t)*vecs[1,1])-(vecs[1,1]*exp(lambdas[2]*t)*vecs[1,2]))/
                  (vecs[1,2]*vecs[2,1]-vecs[2,2]*vecs[1,1]), #P[1,2]
                Sa*((vecs[1,2]*exp(lambdas[1]*t)*vecs[2,1])-(vecs[1,1]*exp(lambdas[2]*t)*vecs[2,2]))/
                  (vecs[1,2]*vecs[2,1]-vecs[2,2]*vecs[1,1])  #P[2,2]
  ), # vecs [x,y] is 'xth' element of 'yth' eigenvector because vectors are stored in columns (opposite order in equation notations)
  nrow=2, ncol=2) 
} 

# produce a fitness landscape matrix by calculating a vertical gradient of maximum eigenvalue of P associated with a gradient of mu per T, and 2) repeating across T
# (as I scan through the gradient of mu, all other life history parameters (d, f, gamma) change according to trade-off assumptions, thus creating a different P for each value of mu) 
fulltradeoff_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rev(f_seq)[i], rev(gamma_seq)[i]) # as rate of maturity goes up -> juv mortality up; fecundity down; adult mortality down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa)
    fulltradeoff_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

# restructure landscape matrix in order to plot in heatmap form, normalized by column (per "T")
colnames(fulltradeoff_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fulltradeoff_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
fulltradeoff_mu_fitness_scaled <-  apply(fulltradeoff_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fulltradeoff_mu_melted <- melt(fulltradeoff_mu_fitness_scaled)

# track optimal life history value (maximum fitness) across T
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fulltradeoff_mu_max <- apply(fulltradeoff_mu_fitness[,1:ncol(fulltradeoff_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fulltradeoff_mu_max <- as.data.frame(fulltradeoff_mu_max)

# plot scaled fitness landscape for mu across T
mu_lndsc.p <- ggplot(fulltradeoff_mu_melted, aes(x=Var2,y=Var1))+
  geom_tile(aes(fill=value))+
  ggtitle(expression(paste("Fitness Landscape of ", mu)))+
  scale_fill_viridis("Relative Fitness \n", option="inferno")+
  scale_x_continuous(breaks=seq(0,14,2), limits=c(0,14))+
  xlab("Disturbance Period (T)")+ #x-axis will stay the same for all plots
  ylab(expression(paste(mu)))+ #ylab depends on function used from above (turn on this line, or line below)
  geom_line(data=fulltradeoff_mu_max, aes(x=c(t_seq[1:length(t_seq)]), y=c(mu_seq[fulltradeoff_mu_max]), color="Optimal life history"), size=1) + scale_color_manual("",values= ("max"="black"))+
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        plot.title=element_text(size=16, face="bold", hjust=0.5), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=15),
        axis.line = element_blank())

mu_lndsc.p



## Below shows identical procedure for creating fitness landscape of f (by scanning across 'f_seq' instead of 'mu_seq'), using the same 'M.fulltradeoff' and 'P.calc' functions

fulltradeoff_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rev(d_seq)[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fulltradeoff_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fulltradeoff_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fulltradeoff_f_fitness) <- c(f_seq[1:length(f_seq)])
fulltradeoff_f_fitness_scaled <-  apply(fulltradeoff_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fulltradeoff_f_melted <- melt(fulltradeoff_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fulltradeoff_f_max <- apply(fulltradeoff_f_fitness[,1:ncol(fulltradeoff_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fulltradeoff_f_max <- as.data.frame(fulltradeoff_f_max)

f_lndsc.p <- ggplot(fulltradeoff_f_melted, aes(x=Var2,y=Var1))+
  geom_tile(aes(fill=value))+
  ggtitle(expression(paste("Fitness Landscape of ", italic(f))))+
  scale_fill_viridis("Relative Fitness \n", option="inferno")+
  scale_x_continuous(breaks=seq(0,14,2), limits=c(0,14)) +
  xlab("Disturbance Period (T)")+ #x-axis will stay the same for all plots
  ylab("f")+ #ylab depends on function used from above (turn on this line, or line below)
  geom_line(data=fulltradeoff_f_max, aes(x=c(t_seq[1:length(t_seq)]), y=c(f_seq[fulltradeoff_f_max]), color="Optimal life history"), size=1) + scale_color_manual("",values= ("max"="black"))+
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        plot.title=element_text(face = "plain", size=16, hjust=0.5), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=15),
        axis.line = element_blank())

f_lndsc.p




## Below produces Figure 2 of paper

#function to share legend without squishing one plot
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}


fig2 <- grid_arrange_shared_legend(mu_lndsc.p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), 
                                   f_lndsc.p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), 
                                   nrow = 1, position = "right")
