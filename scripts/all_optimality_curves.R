library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(gtable)
library(latex2exp)

mu_seq <- seq(0.01, 0.05, length.out=100)  
d_seq <- seq(0.01,0.4, length.out=100)
f_seq <- seq(4,17, length.out=100) 
gamma_seq <- seq(0.01,0.2, length.out = 100)
t_seq<- seq(0.5,14,length.out = 100)
Sj = 0.3
Sa = 0.9


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
  ), # vecs [x,y] is 'xth' element of 'yth' eigenvector (because vectors are stored in columns) (opposite order in equation notations)
  nrow=2, ncol=2) 
} 
M.fulltradeoff <- function(mu, d, f, gamma){
  M <- matrix(c(-(mu+d),
                mu, 
                f, 
                -gamma), 
              2,2)
}


####1: mu ~ d####
md_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rep(10.5,length(f_seq))[i], rep(0.1, length(gamma_seq))[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    md_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(md_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(md_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
md_mu_fitness_scaled <-  apply(md_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
md_mu_melted <- melt(md_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
md_mu_max <- apply(md_mu_fitness[,1:ncol(md_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
md_mu_max <- as.data.frame(md_mu_max)
md.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[md_mu_max$md_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


md_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rep(10.5, length(f_seq))[i], rep(0.1, length(gamma_seq))[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    md_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(md_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(md_f_fitness) <- c(f_seq[1:length(f_seq)])
md_f_fitness_scaled <-  apply(md_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
md_f_melted <- melt(md_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
md_f_max <- apply(md_f_fitness[,1:ncol(md_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
md_f_max <- as.data.frame(md_f_max)
md.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[md_f_max$md_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####2: mu ~ f####
mf_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], rep(0.2, length(d_seq))[i], rev(f_seq)[i], rep(0.1, length(gamma_seq))[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mf_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mf_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mf_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mf_mu_fitness_scaled <-  apply(mf_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mf_mu_melted <- melt(mf_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mf_mu_max <- apply(mf_mu_fitness[,1:ncol(mf_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mf_mu_max <- as.data.frame(mf_mu_max)
mf.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mf_mu_max$mf_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mf_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rep(0.2, length(d_seq))[i], f_seq[i], rep(0.1, length(gamma_seq))[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mf_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mf_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mf_f_fitness) <- c(f_seq[1:length(f_seq)])
mf_f_fitness_scaled <-  apply(mf_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mf_f_melted <- melt(mf_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mf_f_max <- apply(mf_f_fitness[,1:ncol(mf_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mf_f_max <- as.data.frame(mf_f_max)
mf.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mf_f_max$mf_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


####3: mu ~ gamma####
mg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], rep(0.2, length(d_seq))[i], rep(10.5, length(f_seq))[i], gamma_seq[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mg_mu_fitness_scaled <-  apply(mg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mg_mu_melted <- melt(mg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mg_mu_max <- apply(mg_mu_fitness[,1:ncol(mg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mg_mu_max <- as.data.frame(mg_mu_max)
mg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mg_mu_max$mg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(mu_seq[i], rep(0.2, length(d_seq))[i], rep(10.5, length(f_seq))[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mg_f_fitness) <- c(f_seq[1:length(f_seq)])
mg_f_fitness_scaled <-  apply(mg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mg_f_melted <- melt(mg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mg_f_max <- apply(mg_f_fitness[,1:ncol(mg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mg_f_max <- as.data.frame(mg_f_max)
mg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mg_f_max$mg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####4: d ~ f####
df_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], f_seq[i], rep(0.1, length(gamma_seq))[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    df_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(df_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(df_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
df_mu_fitness_scaled <-  apply(df_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
df_mu_melted <- melt(df_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
df_mu_max <- apply(df_mu_fitness[,1:ncol(df_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
df_mu_max <- as.data.frame(df_mu_max)
df.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[df_mu_max$df_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


df_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], f_seq[i], rep(0.1, length(gamma_seq))[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    df_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(df_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(df_f_fitness) <- c(f_seq[1:length(f_seq)])
df_f_fitness_scaled <-  apply(df_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
df_f_melted <- melt(df_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
df_f_max <- apply(df_f_fitness[,1:ncol(df_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
df_f_max <- as.data.frame(df_f_max)
df.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[df_f_max$df_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####5: d ~ gamma####
dg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], rep(10.5, length(f_seq))[i], rev(gamma_seq)[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    dg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(dg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(dg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
dg_mu_fitness_scaled <-  apply(dg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
dg_mu_melted <- melt(dg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
dg_mu_max <- apply(dg_mu_fitness[,1:ncol(dg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
dg_mu_max <- as.data.frame(dg_mu_max)
dg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[dg_mu_max$dg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


dg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], rep(10.5, length(f_seq))[i], rev(gamma_seq)[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    dg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(dg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(dg_f_fitness) <- c(f_seq[1:length(f_seq)])
dg_f_fitness_scaled <-  apply(dg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
dg_f_melted <- melt(dg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
dg_f_max <- apply(dg_f_fitness[,1:ncol(dg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
dg_f_max <- as.data.frame(dg_f_max)
dg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[dg_f_max$dg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####6: f ~ gamma####
fg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], rep(10.5, length(f_seq))[i], rev(gamma_seq)[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
fg_mu_fitness_scaled <-  apply(fg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fg_mu_melted <- melt(fg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fg_mu_max <- apply(fg_mu_fitness[,1:ncol(fg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fg_mu_max <- as.data.frame(fg_mu_max)
fg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[fg_mu_max$fg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


fg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], rep(10.5, length(f_seq))[i], rev(gamma_seq)[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fg_f_fitness) <- c(f_seq[1:length(f_seq)])
fg_f_fitness_scaled <-  apply(fg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fg_f_melted <- melt(fg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fg_f_max <- apply(fg_f_fitness[,1:ncol(fg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fg_f_max <- as.data.frame(fg_f_max)
fg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[fg_f_max$fg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####7: mu ~ d & mu ~ f####
mdmf_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rev(f_seq)[i], rep(0.1, length(gamma_seq))[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdmf_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdmf_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdmf_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mdmf_mu_fitness_scaled <-  apply(mdmf_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdmf_mu_melted <- melt(mdmf_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdmf_mu_max <- apply(mdmf_mu_fitness[,1:ncol(mdmf_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdmf_mu_max <- as.data.frame(mdmf_mu_max)
mdmf.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mdmf_mu_max$mdmf_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mdmf_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rev(d_seq)[i], f_seq[i], rep(0.1, length(gamma_seq))[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdmf_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdmf_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdmf_f_fitness) <- c(f_seq[1:length(f_seq)])
mdmf_f_fitness_scaled <-  apply(mdmf_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdmf_f_melted <- melt(mdmf_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdmf_f_max <- apply(mdmf_f_fitness[,1:ncol(mdmf_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdmf_f_max <- as.data.frame(mdmf_f_max)
mdmf.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mdmf_f_max$mdmf_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")



####8: mu ~ d & mu ~ gamma####
mdmg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rep(10.5, length(f_seq))[i], gamma_seq[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdmg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdmg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdmg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mdmg_mu_fitness_scaled <-  apply(mdmg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdmg_mu_melted <- melt(mdmg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdmg_mu_max <- apply(mdmg_mu_fitness[,1:ncol(mdmg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdmg_mu_max <- as.data.frame(mdmg_mu_max)
mdmg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mdmg_mu_max$mdmg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mdmg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rep(10.5, length(f_seq))[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdmg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdmg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdmg_f_fitness) <- c(f_seq[1:length(f_seq)])
mdmg_f_fitness_scaled <-  apply(mdmg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdmg_f_melted <- melt(mdmg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdmg_f_max <- apply(mdmg_f_fitness[,1:ncol(mdmg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdmg_f_max <- as.data.frame(mdmg_f_max)
mdmg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mdmg_f_max$mdmg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####9: mu ~ f & mu ~ gamma####
mfmg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], rep(0.2, length(d_seq))[i], rev(f_seq)[i], gamma_seq[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mfmg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mfmg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mfmg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mfmg_mu_fitness_scaled <-  apply(mfmg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mfmg_mu_melted <- melt(mfmg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mfmg_mu_max <- apply(mfmg_mu_fitness[,1:ncol(mfmg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mfmg_mu_max <- as.data.frame(mfmg_mu_max)
mfmg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mfmg_mu_max$mfmg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mfmg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rep(0.2, length(d_seq))[i], f_seq[i], rev(gamma_seq)[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mfmg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mfmg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mfmg_f_fitness) <- c(f_seq[1:length(f_seq)])
mfmg_f_fitness_scaled <-  apply(mfmg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mfmg_f_melted <- melt(mfmg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mfmg_f_max <- apply(mfmg_f_fitness[,1:ncol(mfmg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mfmg_f_max <- as.data.frame(mfmg_f_max)
mfmg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mfmg_f_max$mfmg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####10: f ~ mu & f ~ d####
fmfd_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], rev(d_seq)[i], rev(f_seq)[i], rep(0.1, length(gamma_seq))[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fmfd_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fmfd_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fmfd_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
fmfd_mu_fitness_scaled <-  apply(fmfd_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fmfd_mu_melted <- melt(fmfd_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fmfd_mu_max <- apply(fmfd_mu_fitness[,1:ncol(fmfd_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fmfd_mu_max <- as.data.frame(fmfd_mu_max)
fmfd.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[fmfd_mu_max$fmfd_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


fmfd_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], d_seq[i], f_seq[i], rep(0.1, length(gamma_seq))[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fmfd_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fmfd_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fmfd_f_fitness) <- c(f_seq[1:length(f_seq)])
fmfd_f_fitness_scaled <-  apply(fmfd_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fmfd_f_melted <- melt(fmfd_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fmfd_f_max <- apply(fmfd_f_fitness[,1:ncol(fmfd_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fmfd_f_max <- as.data.frame(fmfd_f_max)
fmfd.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[fmfd_f_max$fmfd_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####11: f ~ mu & f ~ gamma####
fmfg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], rep(0.2, length(d_seq))[i], rev(f_seq)[i], rev(gamma_seq)[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fmfg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fmfg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fmfg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
fmfg_mu_fitness_scaled <-  apply(fmfg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fmfg_mu_melted <- melt(fmfg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fmfg_mu_max <- apply(fmfg_mu_fitness[,1:ncol(fmfg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fmfg_mu_max <- as.data.frame(fmfg_mu_max)
fmfg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[fmfg_mu_max$fmfg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


fmfg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rep(0.2, length(d_seq))[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fmfg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fmfg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fmfg_f_fitness) <- c(f_seq[1:length(f_seq)])
fmfg_f_fitness_scaled <-  apply(fmfg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fmfg_f_melted <- melt(fmfg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fmfg_f_max <- apply(fmfg_f_fitness[,1:ncol(fmfg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fmfg_f_max <- as.data.frame(fmfg_f_max)
fmfg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[fmfg_f_max$fmfg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####12: f ~ d & f ~ gamma####
fdfg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], f_seq[i], gamma_seq[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fdfg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fdfg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fdfg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
fdfg_mu_fitness_scaled <-  apply(fdfg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fdfg_mu_melted <- melt(fdfg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fdfg_mu_max <- apply(fdfg_mu_fitness[,1:ncol(fdfg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fdfg_mu_max <- as.data.frame(fdfg_mu_max)
fdfg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[fdfg_mu_max$fdfg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


fdfg_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rep(0.03, length(mu_seq))[i], d_seq[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    fdfg_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(fdfg_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(fdfg_f_fitness) <- c(f_seq[1:length(f_seq)])
fdfg_f_fitness_scaled <-  apply(fdfg_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
fdfg_f_melted <- melt(fdfg_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
fdfg_f_max <- apply(fdfg_f_fitness[,1:ncol(fdfg_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
fdfg_f_max <- as.data.frame(fdfg_f_max)
fdfg.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[fdfg_f_max$fdfg_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")

####13: mu ~ d & f ~ gamma & mu ~ f####
mdfgmf_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rev(f_seq)[i], rev(gamma_seq)[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdfgmf_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdfgmf_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdfgmf_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mdfgmf_mu_fitness_scaled <-  apply(mdfgmf_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdfgmf_mu_melted <- melt(mdfgmf_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdfgmf_mu_max <- apply(mdfgmf_mu_fitness[,1:ncol(mdfgmf_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdfgmf_mu_max <- as.data.frame(mdfgmf_mu_max)
mdfgmf.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mdfgmf_mu_max$mdfgmf_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mdfgmf_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rev(d_seq)[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mdfgmf_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mdfgmf_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mdfgmf_f_fitness) <- c(f_seq[1:length(f_seq)])
mdfgmf_f_fitness_scaled <-  apply(mdfgmf_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mdfgmf_f_melted <- melt(mdfgmf_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mdfgmf_f_max <- apply(mdfgmf_f_fitness[,1:ncol(mdfgmf_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mdfgmf_f_max <- as.data.frame(mdfgmf_f_max)
mdfgmf.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mdfgmf_f_max$mdfgmf_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")




####14: mu ~ d & d ~ gamma & mu ~ f ####
mddgmf_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rev(f_seq)[i], rev(gamma_seq)[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mddgmf_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mddgmf_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mddgmf_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mddgmf_mu_fitness_scaled <-  apply(mddgmf_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mddgmf_mu_melted <- melt(mddgmf_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mddgmf_mu_max <- apply(mddgmf_mu_fitness[,1:ncol(mddgmf_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mddgmf_mu_max <- as.data.frame(mddgmf_mu_max)
mddgmf.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mddgmf_mu_max$mddgmf_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mddgmf_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rev(d_seq)[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mddgmf_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mddgmf_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mddgmf_f_fitness) <- c(f_seq[1:length(f_seq)])
mddgmf_f_fitness_scaled <-  apply(mddgmf_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mddgmf_f_melted <- melt(mddgmf_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mddgmf_f_max <- apply(mddgmf_f_fitness[,1:ncol(mddgmf_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mddgmf_f_max <- as.data.frame(mddgmf_f_max)
mddgmf.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mddgmf_f_max$mddgmf_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")



####15: mu ~ d & d ~ f $ mu ~ gamma####
mddfmg_mu_fitness <- matrix(NA, nrow=length(mu_seq), ncol=length(t_seq)) 
for (i in 1:length(mu_seq)){
  M <- M.fulltradeoff(mu_seq[i], d_seq[i], rev(f_seq)[i], gamma_seq[i]) #as rate of mat goes up -> 'juv mort' up; 'fecund' down; 'adult mort' down
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mddfmg_mu_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mddfmg_mu_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mddfmg_mu_fitness) <- c(mu_seq[1:length(mu_seq)])
mddfmg_mu_fitness_scaled <-  apply(mddfmg_mu_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mddfmg_mu_melted <- melt(mddfmg_mu_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mddfmg_mu_max <- apply(mddfmg_mu_fitness[,1:ncol(mddfmg_mu_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mddfmg_mu_max <- as.data.frame(mddfmg_mu_max)
mddfmg.muplot <- ggplot() + geom_line(aes(x=t_seq, y=mu_seq[mddfmg_mu_max$mddfmg_mu_max]), size = 0.5, col = "red") + ylim(0.01, 0.05) +
  ggtitle(expression(paste("Optimal ", mu, " ~ T"))) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")


mddgmf_f_fitness <- matrix(NA, nrow=length(f_seq), ncol=length(t_seq)) 
for (i in 1:length(f_seq)){
  M <- M.fulltradeoff(rev(mu_seq)[i], rev(d_seq)[i], f_seq[i], gamma_seq[i]) #as fecund goes up -> 'rate of mat' down; 'juv mort' down; 'adult mort' up
  for (t in 1:length(t_seq)){
    P <- P.calc(M,t_seq[t],Sj, Sa) #mortality rates
    mddgmf_f_fitness[i,t] = max(Re(eigen(P)$values))
  }
}

colnames(mddgmf_f_fitness) <- c(t_seq[1:length(t_seq)])
rownames(mddgmf_f_fitness) <- c(f_seq[1:length(f_seq)])
mddgmf_f_fitness_scaled <-  apply(mddgmf_f_fitness, MARGIN=2, FUN = function(x) (x-min(x))/diff(range(x)))#normalize by column
mddgmf_f_melted <- melt(mddgmf_f_fitness_scaled)
max_id <- function(fit_col){which(fit_col == max(fit_col))}
mddgmf_f_max <- apply(mddgmf_f_fitness[,1:ncol(mddgmf_f_fitness)],2,max_id) #find, for each T, the maximum value of given variable
mddgmf_f_max <- as.data.frame(mddgmf_f_max)
mddgmf.fplot <- ggplot() + geom_line(aes(x=t_seq, y=f_seq[mddgmf_f_max$mddgmf_f_max]), size = 0.5, col = "red") + ylim(4,17) +
  ggtitle("Optimal f ~ T") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust=0.5), legend.position="none")




#### Stitching Plots #####
Model1 <- grid.arrange(md.muplot + theme(plot.title = element_blank()), md.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('\u03BC', '\u2194', "d"))), left = textGrob("A", gp=gpar(fontface="bold")))
Model2 <- grid.arrange(mf.muplot + theme(plot.title = element_blank()), mf.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('\u03BC', '\u2194', "f"))), left = textGrob("B", gp=gpar(fontface="bold")))
Model3 <- grid.arrange(mg.muplot + theme(plot.title = element_blank()), mg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('\u03BC', '\u2194', '\u03B3'))), left = textGrob("C", gp=gpar(fontface="bold")))
Model4 <- grid.arrange(df.muplot + theme(plot.title = element_blank()), df.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("d", '\u2194', "f"))), left = textGrob("D", gp=gpar(fontface="bold")))
Model5 <- grid.arrange(dg.muplot + theme(plot.title = element_blank()), dg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("d", '\u2194', '\u03B3'))), left = textGrob("E", gp=gpar(fontface="bold")))
Model6 <- grid.arrange(fg.muplot + theme(plot.title = element_blank()), fg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("f", '\u2194', '\u03B3'))), left = textGrob("F", gp=gpar(fontface="bold")))
Model7 <- grid.arrange(mdmf.muplot + theme(plot.title = element_blank()), mdmf.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("d", '\u2194', '\u03BC', '\u2194', "f"))), left = textGrob("G", gp=gpar(fontface="bold")))
Model8 <- grid.arrange(mdmg.muplot + theme(plot.title = element_blank()), mdmg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("d", '\u2194', '\u03BC', '\u2194', '\u03B3'))), left = textGrob("H", gp=gpar(fontface="bold")))
Model9 <- grid.arrange(mfmg.muplot + theme(plot.title = element_blank()), mfmg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("f", '\u2194', '\u03BC', '\u2194', '\u03B3'))), left = textGrob("I  ", gp=gpar(fontface="bold")))
Model10 <- grid.arrange(fmfd.muplot + theme(plot.title = element_blank()), fmfd.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('\u03BC', '\u2194', "f", '\u2194', "d"))), left = textGrob("J ", gp=gpar(fontface="bold")))
Model11 <- grid.arrange(fmfg.muplot + theme(plot.title = element_blank()), fmfg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('\u03BC', '\u2194', "f", '\u2194', '\u03B3'))), left = textGrob("K", gp=gpar(fontface="bold")))
Model12 <- grid.arrange(fdfg.muplot + theme(plot.title = element_blank()), fdfg.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste('d', '\u2194', "f", '\u2194', '\u03B3'))), left = textGrob("L", gp=gpar(fontface="bold")))
Model13 <- grid.arrange(mdfgmf.muplot + theme(plot.title = element_blank()), mdfgmf.fplot + theme(plot.title = element_blank()), ncol = 2, top = textGrob(sprintf(paste("d", '\u2194', '\u03BC', '\u2194', "f", '\u2194', '\u03B3'))), left = textGrob("M", gp=gpar(fontface="bold")))


single_tradeoffs <- grid.arrange(Model1, Model2, Model3, Model4, Model5, Model6,
                                 ncol =1,
                                 top=textGrob("Models with Single Tradeoffs", gp=gpar(fontsize=12,font=4)),
                                 bottom=textGrob(expression(paste("    Optimal ", mu, "~T          Optimal f~T")),gp=gpar(fontsize=13,font=3)))

two_tradeoffs <- grid.arrange(Model7, Model8, Model9, Model10, Model11, Model12,
                              ncol =1,
                              top=textGrob("Models with Secondary Tradeoffs", gp=gpar(fontsize=12,font=4)),
                              bottom=textGrob(expression(paste("    Optimal ", mu, "~T          Optimal f~T")),gp=gpar(fontsize=13,font=3)))

blank <- grid.rect(gp=gpar(col="white"))
three_tradeoffs <- grid.arrange(Model13, textGrob(expression(paste("    Optimal ", mu, "~T          Optimal f~T")),gp=gpar(fontsize=13,font=3), vjust = -2.7), blank, blank, blank, blank, 
                                ncol =1,
                                top = textGrob("Model with Tertiary Tradeoffs", gp=gpar(fontsize=12,font=4)))

all_plot <- plot_grid(single_tradeoffs, NULL, two_tradeoffs, NULL, three_tradeoffs,
                      ncol = 5,
                      rel_widths = c(1,0.1,1,0.1,1))

# export picture at 850W x 900H