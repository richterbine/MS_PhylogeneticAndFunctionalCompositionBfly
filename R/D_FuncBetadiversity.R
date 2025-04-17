#################################################
########## Functional - betadiversity ###########
#################################################

# packages required
library(phytools)
library(SYNCSA)
library(PVR)
library(FD)
library(PCPS)
library(openxlsx)
library(tidyverse)

### Functions ####
source(here::here("R", "functions", "function_ADONIS2.R"))
source(here::here("R", "functions", "EM_mantel_function.R"))
source(here::here("R", "functions", "CWM_glmer_function_scaled.R"))
source(here::here("R", "functions", "CWM_physig_function_scaled.R"))

###### read data #####
traits.bfly <- read.xlsx(here::here("data/processed/Mean_traits_bfly_new.xlsx"))
head(traits.bfly)

# read the phylogenetic tree
tree.bfly <- readRDS(here::here("data/processed/tree_Flona.rds"))
tree.bfly$tip.label

# read the environmental and species data
comm.data <- read.xlsx(here::here("data/processed/data.bfly.xlsx"))

comm.data <- comm.data %>% 
  mutate(Strata = substr(Code, 1, 1)) %>% 
  mutate(Month = substr(Code, 5, 5))

# verifying the oscilation in temperature between strata
comm.data %>% 
 ggplot() +
  geom_boxplot(aes(x = Strata, y = Tmean), notch = T)

comm.data %>% 
  ggplot() +
  geom_boxplot(aes(x = Strata, y = Tsd), notch = T)

# separating the environmental variable "Temperature"
comm.env <- comm.data %>% column_to_rownames("Code") %>% 
  select(Strata, Month)
comm.env %>% head

write.xlsx(comm.env, here::here("data/processed/comm_envir.xlsx"))

# separating the community matrix
comm.bfly <- comm.data %>% 
  select(Archaeoprepona_amphimachus:Zaretis_strigosus)
rownames(comm.bfly) <- comm.data$Code
comm.bfly %>% head

write.xlsx(comm.bfly, here::here("data/processed/community_bfly.xlsx"))

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
traits.bfly$FWL <- residuals(fwl.af)
traits.bfly$FEA <- residuals(fea.af)
traits.bfly$AM <- residuals(am.af)
traits.bfly$TM <- residuals(tm.af)

# define rownames
rownames(traits.bfly) <- traits.bfly$Binominal_name
traits.bfly %>% head

traits <- traits.bfly %>% 
  select(-Binominal_name)

traits %>% glimpse()

# define iridescent and eyespot presence as factor
traits$Iridescence_presence_Dorsal <- as.integer(traits$Iridescence_presence_Dorsal)
traits$Eyespot_presence_Ventral <- as.integer(traits$Eyespot_presence_Ventral)

write.xlsx(traits, here::here("data/processed/traits.xlsx"))

traits %>% rownames_to_column("Species") %>% 
  pivot_longer(cols = Mean_contrast_Dorsal:Host,
               names_to = "Trait", 
               values_to = "Values") %>% 
  ggplot() +
  geom_histogram(aes(x = Values)) +
  facet_wrap(~Trait, scales = "free_x")

########################################################################
# function that incorporate the phylogenetic signal, phylogeny correction
# and CWM for each trait

RES <- physig.cwm.function_scaled(comm = comm.bfly, traits = traits, tree = tree.bfly,
                           envir = comm.env, binary = TRUE, runs = 999, 
                           AsFactors = c(1,2), predictor = "Strata", 
                           random.term = "Month")

saveRDS(RES, here::here("output/res_cwm_new.rds"))

# Extracting the results from CWM function ---- 
RES <- readRDS(here::here("output/res_cwm_new.rds"))

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
(table2)

# write.xlsx(table2, here::here("data/processed/Table2_CWM.xlsx"))


# drawing plots -----------------------------------------------------------

traits <- read.xlsx(here::here("data/processed/traits.xlsx")) %>% 
  bind_cols(traits.bfly %>% select(Binominal_name)) %>% 
  select(-Binominal_name)
traits %>% colnames()

comm.bfly <- read.xlsx(here::here("data/processed/community_bfly.xlsx"))

org <- organize.syncsa(comm.bfly, traits = traits)

cwm.bfly <- SYNCSA::matrix.t(org$community, org$traits, scale = F)$matrix.T

df.cwm <- cwm.bfly %>% scale() %>% as.data.frame() %>% 
  rownames_to_column("Code") #%>% 

df.cwm <- df.cwm %>%   
pivot_longer(cols = Mean_contrast_Dorsal:Host,
               names_to = "traits", 
               values_to = "values") %>% 
  mutate(Strata = substr(Code, 1,1), .before = traits) %>% 
  mutate(Month = substr(Code, 5,5), .before = traits) %>% 
  mutate(traits.new = str_replace_all(traits, "Dorsal", "D")) %>% 
  mutate(traits.new = str_replace_all(traits.new, "Ventral", "V")) %>% 
  mutate(traits.new = str_replace_all(traits.new, "luminance", "brightness")) %>% 
  mutate(traits.new = str_replace_all(traits.new, "_presence", "")) %>% 
  mutate(traits.new = str_replace_all(traits.new, "Shannon_cd", "Color_diversity"))

df.cwm %>% head

table2 %>% rownames_to_column("traits") %>% 
  mutate(Significance = case_when(p.site.shuffle.s < 0.05 & 
                                    p.trait.shuffle.s < 0.05 ~ "**",
                                  p.site.shuffle.s < 0.05 & 
                                    p.trait.shuffle.s > 0.05 ~ "*",
                                  .default = "")) %>% 
  mutate(Type = factor(case_when(grepl("Mean_", traits) ~ "Ecophysiological",
                                 grepl("Shannon", traits) ~ "Ecophysiological",
                                 grepl("_presence", traits) ~ "Defense Strategy",
                                 grepl("FEA", traits) ~ "Habitat Perception",
                                 grepl("Host", traits) ~ "Habitat Perception",
                                 .default = "Flight Performance"))) %>% 
  write.xlsx(here::here("output/Table2_CWM.xlsx"))

# Criando um labeller que substitui "_" por " "
custom_labeller <- as_labeller(function(x) str_replace_all(x, "_", " "))

p.cwm <- table2 %>% rownames_to_column("traits") %>% 
  mutate(Significance = case_when(p.site.shuffle.s < 0.05 & 
                                    p.trait.shuffle.s < 0.05 ~ "**",
                                  p.site.shuffle.s < 0.05 & 
                                    p.trait.shuffle.s > 0.05 ~ "*",
                                  .default = "")) %>% 
  select(traits, Significance) %>% 
  mutate(Type = factor(case_when(grepl("Mean_", traits) ~ "Ecophysiological",
                          grepl("Shannon", traits) ~ "Ecophysiological",
                          grepl("_presence", traits) ~ "Defense Strategy",
                          grepl("FEA", traits) ~ "Habitat Perception",
                          grepl("Host", traits) ~ "Habitat Perception",
                          .default = "Flight Performance"))) %>% 
  right_join(df.cwm) %>% 
  mutate(traits.new = paste(traits.new, Significance, sep = " ")) %>% 
  select(-Significance) %>% 
  ggplot(aes(x = Strata, y = values, colour = Strata, 
                                   fill = Strata)) +
  geom_boxplot(alpha = .4) + 
  facet_wrap(Type ~ traits.new, labeller = custom_labeller, 
             scales = "free_y") +
  scale_color_manual(values = c("#A84944", "#556F21")) +
  scale_fill_manual(values = c("#A84944", "#556F21")) +
  labs(y = "CWM values") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(labels = c("Canopy", "Understory")) +
  theme(legend.position = "none") 
p.cwm

cowplot::save_plot(here::here("output/figures/Fig2_boxplot2_CWM.png"), p.cwm,
                   base_height = 6, base_width = 8, dpi = 300)


# NTDR curves ----
traits <- read.xlsx(here::here("data/processed/traits.xlsx")) %>% 
  bind_cols(traits.bfly %>% select(Binominal_name)) %>% 
  select(-Binominal_name) %>% 
  dplyr::rename(Contrast_Dorsal = Mean_contrast_Dorsal,
         Contrast_Ventral = Mean_contrast_Ventral,
         Brightness_Ventral = Mean_luminance_Ventral,
         Brightness_Dorsal = Mean_luminance_Dorsal,
         Iridescence_Dorsal = Iridescence_presence_Dorsal,
         Eyespot_Ventral = Eyespot_presence_Ventral,
         Color_Diversity_Dorsal = Shannon_cd_Dorsal,
         Color_Diversity_Ventral = Shannon_cd_Ventral)

# read the phylogenetic tree
tree.bfly <- readRDS(here::here("data/processed/tree_Flona.rds"))
str(tree.bfly)

# matrix of phylogenetic distance
dist.phylo <- cophenetic(tree.bfly)

# organize data
org2 <- organize.syncsa(comm = comm.bfly, phylodist = dist.phylo, 
                        envir = comm.env, traits = traits)

res.ntdr <- list()

for (i in 1:ncol(org2$traits)) {
    res <- pcps.curve(org2$community, org2$phylodist,
                    org2$traits[,i,drop = FALSE], tree = tree.bfly, 
                    null.model.bm = TRUE, runs = 999)
    res.ntdr[[i]] <- res
}

par(mfrow = c(6, 3))

for (i in 1:length(res.ntdr)) {
  
  plot(res.ntdr[[i]], draw.model = "bm", type = "S",
       col = "red", model.col = "grey50",
       probs = c(0.025, 0.975))
  title(paste0(colnames(org2$traits[, i, drop = F])))
}
colnames(traits)



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
