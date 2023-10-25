################################################################
################### Phylo-beta-diversity #######################
################################################################

#install.packages("devtools")
# devtools::install_github("vanderleidebastiani/PCPS", force = T)
# devtools::install_github("vanderleidebastiani/SYNCSA", force = T)
library(SYNCSA)
library(PCPS)
library(phytools)
library(vegan)
library(nlme)
library(ggplot2)
library(ape)


# Read community matrix and the pruned phylogenetic tree ------------------

comm.data <- read.csv(here::here("data/processed/data.bfly.csv"), sep = ",")
str(comm.data)

# separating the environmental variable "Temperature"
comm.env <- comm.data[, c("Tmean", "Tsd"), drop = F]
rownames(comm.env) <- comm.data$Code
comm.env$Strata <- substr(rownames(comm.env), 1, 1)
comm.env$SU <- substr(rownames(comm.env), 2, 4)
comm.env$Month <- substr(rownames(comm.env), 5,5)
head(comm.env)

for (i in 3:length(comm.env)) {
  comm.env[, i] <- as.numeric(as.factor(comm.env[, i]))
} # make the variables as factor (needed for organize.pcps)

str(comm.env)

# separating the community matrix
comm.bfly <- comm.data[,-c(1,2,37, 38)]
str(comm.bfly)
rownames(comm.bfly) <- comm.data$Code
head(comm.bfly)

# read the phylogenetic tree
tree.bfly <- read.tree(here::here("data/processed/tree_bfly_flona.txt"))
str(tree.bfly)

# matrix of phylogenetic distance
dist.phylo <- cophenetic(tree.bfly)

# organizing the data 
org <- organize.pcps(comm = comm.bfly, phylodist = dist.phylo, envir = comm.env)
envir <- org$environmental
comm <- org$community
dist.phylo <- org$phylodist

## extracting the P matrix - community matrix based on phylogenetic information
P_comm<- SYNCSA::matrix.p(as.matrix(comm), as.matrix(dist.phylo))$matrix.P
P_comm

# performing the PCPS
pcps_comm <- pcps(comm, dist.phylo, method = "bray", squareroot = T)
pcps_comm

# Autovectors (PCPS) 
pcps <- pcps_comm$vectors
pcps

# correlation between PCPS axis and matrix.P
pcps_correl <- pcps_comm$correlations
pcps_correl

# score vectors for sites and species for axis with more than 5% of explanation
scores_pcps <- scores.pcps(pcps_comm, choices = c(1, 2)) # axis 1 and 2
scores_pcps3<- scores.pcps(pcps_comm, choices = c(3, 4)) # axis 3 and 4

# plot(pcps_comm, choices = c(1, 3))

## LME - PCPS (Debastiani & Duarte 2014): testing the effect of environment 
# in the phylogenetic composition
res1 <- pcps.sig(comm = comm, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                formula = pcps.1 ~ Strata , envir = envir, 
                random = ~ 1|Month, choices = 1, runs = 999)
res1
summary(res1$model)

res2 <- pcps.sig(comm = comm, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                formula = pcps.2 ~ Strata , envir = envir, 
                random = ~ 1|Month, choices = 2, runs = 999)
res2
summary(res2$model)

res3 <- pcps.sig(comm = comm, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                formula = pcps.3 ~ Strata , envir = envir, 
                random = ~ 1|Month, choices = 3, runs = 999)
res3
summary(res3$model)

res4 <- pcps.sig(comm = comm, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                formula = pcps.4 ~ Strata , envir = envir, 
                random = ~ 1|Month, choices = 4, runs = 999)
res4
summary(res4$model)

pcps.sig.res <- list(res1, res2, res3, res4)
saveRDS(pcps.sig.res, here::here("output/pcps.sig_res.rds"))

# organizing a table with the results of null models
pcps.sig.res <- readRDS(here::here("output/pcps.sig_res.rds"))

table1 <- data.frame(PCPS1 = c(obs = pcps.sig.res[[1]]$obs.statistic, site = pcps.sig.res[[1]]$p.site.shuffle,
                               taxa = pcps.sig.res[[1]]$p.taxa.shuffle, rel_eig = pcps_comm$values$Relative_eig[1]*100),
                     PCPS2 = c(pcps.sig.res[[2]]$obs.statistic, pcps.sig.res[[2]]$p.site.shuffle,
                               pcps.sig.res[[2]]$p.taxa.shuffle, pcps_comm$values$Relative_eig[2]*100),
                     PCPS3 = c(pcps.sig.res[[3]]$obs.statistic, pcps.sig.res[[3]]$p.site.shuffle, 
                               pcps.sig.res[[3]]$p.taxa.shuffle, pcps_comm$values$Relative_eig[3]*100),
                     PCPS4 = c(pcps.sig.res[[4]]$obs.statistic, pcps.sig.res[[4]]$p.site.shuffle, 
                               pcps.sig.res[[4]]$p.taxa.shuffle, pcps_comm$values$Relative_eig[4]*100))

write.csv(round(table1, digits = 3), here::here("data/processed/table1.csv"))

# model summaries
mod.res <- rbind(summary(pcps.sig.res[[1]]$model)[["tTable"]],
                 summary(pcps.sig.res[[2]]$model)[["tTable"]],
                 summary(pcps.sig.res[[3]]$model)[["tTable"]],
                 summary(pcps.sig.res[[4]]$model)[["tTable"]])

mod.res
write.csv(mod.res, here::here("data/processed/summary_models_pcps.csv"))

# Drawing plots for phylogenetic composition ------------------------------

data.raw <- read.csv(here::here("data/processed/data.raw.csv"), sep = ",")
sp.tribe <- data.raw[, c("Tribe", "Species")]
pcps.sp <- unique(sp.tribe)

spp.data <- pcps.sp$Tribe[match(rownames(scores_pcps$scores.species), pcps.sp$Species)]

pcps.data <- data.frame(Scores = c(rep("Sites", nrow(scores_pcps$scores.sites)),
                                  rep("Species", nrow(scores_pcps$scores.species))),
                       Type = c(rep("Canopy", 30), rep("Understory", 28), spp.data),
                       PCPS1 = c(scores_pcps$scores.sites[, 1], scores_pcps$scores.species[, 1]),
                       PCPS2 = c(scores_pcps$scores.sites[, 2], scores_pcps$scores.species[, 2]),
                       PCPS3 = c(scores_pcps3$scores.sites[, 1], scores_pcps3$scores.species[, 1]),
                       PCPS4 = c(scores_pcps3$scores.sites[, 2], scores_pcps3$scores.species[, 2]))
head(pcps.data)


PCPS_1.2 <- ggplot() + 
  geom_point(data = subset(pcps.data, Scores == "Sites"),
             aes(x = PCPS1, y= PCPS2, shape = Type), size = 3, color = "gray60") +
  scale_shape_manual(values = c(0, 1)) +
  geom_point(data = subset(pcps.data, Scores == "Species"),
             aes(x = PCPS1, y= PCPS2, colour = Type), 
             size = 5, alpha = .7, shape = 18) + 
  scale_color_viridis_d(option = "D", 
                        limits = c("Nymphalini", "Ageroniini",
                                   "Epiphilini", "Callicorini",
                                   "Preponini", "Anaeini",
                                   "Brassolini", "Morphini", "Satyrini")) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  theme_classic() + labs(tag = "a)")

PCPS_1.2


PCPS_1.3 <- ggplot() + 
  geom_point(data = subset(pcps.data, Scores == "Sites"),
             aes(x = PCPS1, y = PCPS3, shape = Type), size = 3, color = "gray60") +
  scale_shape_manual(values = c(0, 1)) +
  geom_point(data = subset(pcps.data, Scores == "Species"),
             aes(x = PCPS1, y = PCPS3, colour = Type), 
             size = 5, alpha = .7, shape = 18) + 
  scale_color_viridis_d(option = "D", limits = c("Nymphalini", "Ageroniini",
                                                 "Epiphilini", "Callicorini",
                                                 "Preponini", "Anaeini",
                                                 "Brassolini", "Morphini", "Satyrini")) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  theme_classic() + labs(tag = "b)")

PCPS_1.3

PCPS_2.3 <- ggplot() + 
  geom_point(data = subset(pcps.data, Scores == "Sites"),
             aes(x = PCPS2, y = PCPS3, shape = Type),  
             size = 3, color = "gray60") +
  scale_shape_manual(values = c(0, 1)) +
  geom_point(data = subset(pcps.data, Scores == "Species"),
             aes(x = PCPS2, y = PCPS3, colour = Type), 
             size = 5, alpha = .7, shape = 18) + 
  scale_color_viridis_d(option = "D", limits = c("Nymphalini", "Ageroniini",
                                                 "Epiphilini", "Callicorini",
                                                 "Preponini", "Anaeini",
                                                 "Brassolini", "Morphini", "Satyrini")) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  theme_classic() + labs(tag = "c)")

PCPS_2.3

leg <- cowplot::get_legend(PCPS_1.2)

p.PCPS <- cowplot::plot_grid(PCPS_1.2 + theme(legend.position = "none"), 
                   PCPS_1.3 + theme(legend.position = "none"), 
                   PCPS_2.3 + theme(legend.position = "none"),
                   leg)

cowplot::save_plot(here::here("output/figures/Fig1_PCPS1.png"), p.PCPS,
                   base_height = 6, base_width = 8)
