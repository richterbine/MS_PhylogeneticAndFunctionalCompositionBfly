#################################################
########## Functional - betadiversity ###########
#################################################

# packages required
library(phytools)
library(SYNCSA)
library(PVR)
library(FD)
library(PCPS)

### Functions ####
source(here::here("R", "functions", "function_ADONIS2.R"))
source(here::here("R", "functions", "EM_mantel_function.R"))
source(here::here("R", "functions", "CWM_glmer_function.R"))
source(here::here("R", "functions", "CWM_physig_function.R"))

###### read data #####
traits.bfly <- read.csv(here::here("data/processed/Mean_traits_bfly.csv"), sep = ",")
head(traits.bfly)

# read the phylogenetic tree
tree.bfly <- read.tree(here::here("data/processed/tree_bfly_flona.txt"))
tree.bfly$tip.label

# read the environmental and species data
comm.data <- read.csv(here::here("data/processed/data.bfly.csv"), sep = ",")
str(comm.data)

# separating the environmental variable "Temperature"
comm.env <- comm.data[, "Code", drop = F]
rownames(comm.env) <- comm.data$Code
comm.env$Strata <- substr(rownames(comm.env), 1, 1)
comm.env$SU <- substr(rownames(comm.env), 2, 4)
comm.env$Month <- substr(rownames(comm.env), 5,5)
comm.env <- comm.env[,-1]
head(comm.env)

# separating the community matrix
comm.bfly <- comm.data[,-c(1,2,37, 38)]
str(comm.bfly)
rownames(comm.bfly) <- comm.data$Code
head(comm.bfly)

# first - remove the allometric component of FEA and FWL
head(traits.bfly)
cor(traits.bfly$TDM, traits.bfly$FEA)
cor(traits.bfly$TDM, traits.bfly$FWL)
cor(traits.bfly$TDM, traits.bfly$AM)
cor(traits.bfly$TDM, traits.bfly$TM)

fwl.af <- lm(FWL ~ TDM, data = traits.bfly)
fea.af <- lm(FEA ~ TDM, data = traits.bfly)
am.af <- lm(AM ~ TDM, data = traits.bfly)
tm.af <- lm(TM ~ TDM, data = traits.bfly)

plot(residuals(fwl.af), traits.bfly$TDM)
plot(residuals(fea.af), traits.bfly$TDM)
plot(residuals(am.af), traits.bfly$TDM)
plot(residuals(tm.af), traits.bfly$TDM)

# include the allometric free FWL and FEA in the traits matrix
traits.bfly$FWL.af <- residuals(fwl.af)
traits.bfly$FEA.af <- residuals(fea.af)
traits.bfly$AM.af <- residuals(am.af)
traits.bfly$TM.af <- residuals(tm.af)
colnames(traits.bfly)

traits <- traits.bfly[,-c(1,2,8:11)] # removing the traits without allometry-correction
colnames(traits)
rownames(traits) <- traits.bfly$Species
head(traits)
saveRDS(traits, here::here("output/traits.rds"))
saveRDS(comm.bfly, here::here("output/community_bfly.rds"))

########################################################################
# function that incorporate the phylogenetic signal, phylogeny correction
# and CWM for each trait

RES <- physig.cwm.function(comm = comm.bfly, traits = traits, tree = tree.bfly,
                           envir = comm.env, binary = TRUE, runs = 999, 
                           AsFactors = c(1,2,3), predictor = "Strata", 
                           random.term = "Month")

saveRDS(RES, here::here("output/res_cwm.rds"))

table2 <- as.data.frame(t(RES$CWM.output))

table2$Statistic <- na.omit(c(RES$Trait.phylosig$k.statistic[match(rownames(table2), 
                                                        rownames(RES$Trait.phylosig$k.statistic)),1],
                   RES$Trait.phylosig$EM.Mantel[match(rownames(table2), 
                                                rownames(RES$Trait.phylosig$EM.Mantel)),1]))
table2$p.value <- na.omit(c(RES$Trait.phylosig$k.statistic[match(rownames(table2), 
                                                           rownames(RES$Trait.phylosig$k.statistic)),2],
                      RES$Trait.phylosig$EM.Mantel[match(rownames(table2), 
                                                         rownames(RES$Trait.phylosig$EM.Mantel)),2]))
table2$p.BM <- RES$Trait.phylosig$EM.Mantel[match(rownames(table2), 
                                                  rownames(RES$Trait.phylosig$EM.Mantel)),3]
tail(table2)
write.csv(table2, here::here("data/processed/Table2_CWM.csv"))



# drawing plots -----------------------------------------------------------

traits <- readRDS(here::here("output/traits.rds"))
comm.bfly <- readRDS(here::here("output/community_bfly.rds"))

org <- organize.syncsa(comm.bfly, traits = traits)

cwm.bfly <- SYNCSA::matrix.t(org$community, org$traits, scale = F)$matrix.T

df.cwm <- data.frame(values = c(cwm.bfly), 
                     traits = rep(colnames(cwm.bfly), each = nrow(cwm.bfly)),
                     code = rep(rownames(cwm.bfly), ncol(cwm.bfly)))
df.cwm$Strata <- substr(df.cwm$code, 1,1)
df.cwm$Month <- substr(df.cwm$code, 5,5)
head(df.cwm)

library(ggplot2)

df.cwm$traits.org <- factor(x = df.cwm$traits, levels = c("AR","WL","WTr","FWL.af","FW.HW",
                                                      "TDM","AM.af","TM.af","Host","FEA.af",
                                                      "Iridescence","Eyespots","Rings",
                                                      "Disruptive","Criptic","Darkness"),
                        labels = c("AR","WL","WTr","FWL**","FW.HW",
                                   "TDM*","AM","TM","Host*","FEA*",
                                   "Iridescence*","Eyespots","Rings",
                                   "Disruptive*","Criptic","Darkness*"))

p.cwm <- ggplot(data = df.cwm, aes(x = Strata, y = values, colour = Strata, fill = Strata)) +
  geom_boxplot() + facet_wrap(~traits.org, scales = "free_y") +
  scale_color_viridis_d(option = "B") +
  scale_fill_viridis_d(option = "B", alpha = .7) +
  labs(y = "CWM values")
p.cwm

pdf(here::here("output", "figures", "Fig2_boxplot_CWM.pdf"), height = 6, width = 8)
p.cwm
dev.off()


# make the ordination
res.cwm <- list()
res.cwm$T <- cwm.bfly
ord.t <- wcmdscale.org(cwm.bfly, method = "euclid", squareroot = TRUE, 
                       eig = TRUE, correlations = TRUE)
res.cwm$values <- ord.t$values # eigenvalues
res.cwm$vectors <- ord.t$vectors # principal axes

colnames(res.cwm$vectors) <- paste("pcps.", seq_len(ncol(res.cwm$vectors)), 
                               sep = "")
row.names(res.cwm$values) <- colnames(res.cwm$vectors)

res.cwm$correlations <- ord.t$correlations
colnames(res.cwm$correlations) <- paste("pcps.", seq_len(ncol(res.cwm$vectors)), sep = "")
rownames(res.cwm$correlations) <- rownames(res.cwm$correlations, 
                                           do.NULL = FALSE, prefix = "spp.")
class(res.cwm) <- "pcps"

scores.cwm <- scores.pcps(res.cwm)
scores.cwm$scores.species

pcfunc.data <- data.frame(Scores = c(rep("Sites", nrow(scores.cwm$scores.sites)),
                                   rep("Traits", nrow(scores.cwm$scores.species))),
                        Type = c(rep("Canopy", 30), rep("Understory", 28), 
                                 rownames(scores.cwm$scores.species)),
                        PC1 = c(scores.cwm$scores.sites[, 1], scores.cwm$scores.species[, 1]),
                        PC2 = c(scores.cwm$scores.sites[, 2], scores.cwm$scores.species[, 2]))
pcfunc.data$trait.label <- ""
pcfunc.data$trait.label[c(63, 64, 65, 69, 70, 71, 72)] <- pcfunc.data$Type[c(63, 64, 65, 69, 70, 71, 72)]
head(pcfunc.data)

library(ggplot2)
library(ggrepel)

PCFS_1.2 <- ggplot(data = pcfunc.data, aes(x = PC1, y= PC2, label = trait.label)) + 
  geom_point(data = subset(pcfunc.data, Scores == "Sites"),
             aes(x = PC1, y= PC2, shape = Type), size = 3, color = "gray60") +
  scale_shape_manual(values = c(0, 1)) +
  geom_point(data = subset(pcfunc.data, Scores == "Traits"),
             aes(x = PC1, y= PC2, colour = Type), 
             size = 5, alpha = .7, shape = 18) +
  geom_text_repel(max.overlaps = 20, size = 3)+
  scale_color_viridis_d(option = "D", limits = c("AR","WL","WTr","FWL.af","FW.HW",
                                                 "TDM","AM.af","TM.af","Host","FEA.af",
                                                 "Iridescence","Eyespots","Rings",
                                                 "Disruptive","Criptic","Darkness")) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray50", alpha = 0.6) +
  theme_classic() + labs(tag = "a)")

PCFS_1.2

pdf(here::here("output", "figures", "Fig2_PC_CWM.pdf"), height = 6, width = 8)
PCFS_1.2
dev.off()
