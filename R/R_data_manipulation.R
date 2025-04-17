############################################################################
##### MS phylo and functional composition of fruit-feeding butterflies #####
############################################################################

# New edition: 29/10/2024
# New analysis considering reviwers comments (ProcB) and a available color traits dataset to the study area (Thomas et al. (2024))

#' Final edition: 13/03/2025

# Script to read raw data, manipulate the phylogenetic tree and alometric relationship among traits --------
# packages required
library(tidyverse)
library(reshape)
library(openxlsx)

# read the raw data
# species information
data.raw <- read.xlsx(here::here("data/raw/Raw_data_FLONA.xlsx"))
data.raw %>% head

# calculation the sampling effort
# sampling days * number of traps by transect * number of transects
data.raw %>% select(Date) %>% 
  unique() %>% count %>% 
  mutate(SE = n*10*6)

# remove Temenis_laothoe - not catched
data.raw <- data.raw[-which(data.raw$Species == "Temenis_laothoe"),]

write.xlsx(data.raw, here::here("data/processed/data.raw.xlsx"))

# environmental variables
data.env.raw <- read.csv(here::here("data/raw/Temperatures_basetrap.csv"), sep = ";")
data.env.raw %>% head

# generate a column called "Code" at both data frames
data.raw$Code <- paste(substr(data.raw$Strata,1,1), data.raw$SU, 
                       substr(data.raw$Month,1,1), sep = "")

data.raw %>% head()  

tmp <- melt(data.raw, id.vars = c("Code", "Species"), measure.vars = "Abundance")
comm.bfly <- cast(tmp, Code ~ Species, sum)
rownames(comm.bfly) <- comm.bfly$Code
comm.bfly %>% head

data.env.raw$Code <- paste(substr(data.env.raw$Strata,1,1), data.env.raw$SU, 
                       substr(data.env.raw$Month,1,1), sep = "")
data.env.raw %>% head()

tmp <- data.frame(Tmean = round(tapply(data.env.raw$Temperature, 
                                             data.env.raw$Code, mean, na.rm = T), 3),
                  Tsd = round(tapply(data.env.raw$Temperature, 
                                       data.env.raw$Code, sd, na.rm = T), 3))
tmp$Code <- rownames(tmp)
tmp %>% head()

data.bfly <- merge(comm.bfly, tmp, by.x = "Code")
data.bfly %>% head()

# save the community matrix to posterior analisys
write.xlsx(data.bfly, here::here("data/processed/data.bfly.xlsx"))

######################################################################
############ Phylogenetic manipulation - pruning tree ################
######################################################################

# Loading the Nymphalidae tree (Chazot et al. 2019)

source(here::here("R/functions/function_treedata_modif.R"))

# Manipulating the proposed phylogenetic hypothesis
Nym.tree <- ape::read.tree(here::here("data/raw/Nymphalidae_2952spp.new"))
# Nym.tree <- ape::read.tree(here::here("data/raw/TimeTree_phylo_species.nwk"))

# As there are unidentified species of Hermeuptychia in our database, we will use the species H. hermes as a base, which we will call Hermeuptychia_sp in the tree.

Nym.tree$tip.label %>% as_tibble() %>% 
  filter(grepl("Hermeuptychia", Nym.tree$tip.label))


Nym.tree$tip.label[which(Nym.tree$tip.label == "Hermeuptychia_hermes")] <- "Hermeuptychia_sp"

df.bfly <- data.bfly %>% select(where(is.integer)) %>%  colnames() %>% 
  as.data.frame()
colnames(df.bfly) <- "spp"

rownames(df.bfly) <- df.bfly$spp

# Cutting the bfly.tree to my specific sampling pool to Araucaria forest
phy.bfly <- treedata_modif(Nym.tree, df.bfly)$phy

plot(phy.bfly)

saveRDS(phy.bfly, here::here("data/processed/tree_Flona.rds"))


######################################################################
################ Functional manipulation procedures ##################
######################################################################

# here we will include a color dataset produced by Nogueira et al. (2024)
# we selected some traits as brightness, contrast and color heterogeneity to describe the community color composition

ventral <- read.csv(here::here("data/raw/ventral_data.csv"), sep = ";")
ventral$Side <- "Ventral" #adiciona coluna que identifica a região das asas
str(ventral)

dorsal <- read.csv(here::here("data/raw/dorsal_data.csv"), sep = ";")
dorsal$Side <- "Dorsal" #adiciona coluna que identifica a região das asas
str(dorsal)

# merging the dataframes
data.color <- bind_rows(dorsal, ventral)
head(data.color)

# manipulating and organizing variables
data.color$Stratum <- as.factor(data.color$Stratum)
data.color$Strata <- factor(data.color$Stratum, levels = c("Sub_bosque", "Dossel", "Borda", "Pinus_Novo", "Pinus_Velho"),
                            labels = c("Understory", "Canopy", "Edge","Young Pinus", "Old Pinus"))
saveRDS(data.color, here::here("data/processed/org_color_traits.rds"))

# filtering only species present in this study (34 spp)

which(is.na(match(df.bfly$spp, unique(data.color$Binominal_name))))
# unique(data.color$Binominal_name) %>% as.data.frame %>% View()

# some species showing differences between species names, let's change the names of the color data
df.bfly$spp[16] # Euptychoides_castrensis e Moneuptychia_castrensis
df.bfly$spp[21] # Hermeuptychia_sp e Hermeuptychia_sp.
df.bfly$spp[33] # Taygetis_ypthima e Pseudodebis_ypthima

data.color %>% 
  filter(grepl("castrensis", Binominal_name) | grepl("Hermeuptychia", Binominal_name) |
           grepl("ypthima", Binominal_name)) %>% 
  select(Binominal_name) %>% unique()

# renomeando manualmente
data.color <- data.color %>% 
  mutate(Binominal_name = case_when(Binominal_name == "Moneuptychia_castrensis" ~ "Euptychoides_castrensis",
                                    Binominal_name == "Hermeuptychia_sp." ~ "Hermeuptychia_sp",
                                    Binominal_name == "Pseudodebis_ypthima" ~ "Taygetis_ypthima",
                                    .default = as.character(Binominal_name)))
x <- comm.bfly %>% 
  select(where(is.numeric)) %>% 
  colnames()

data.color %>% colnames()

df <- data.color %>% 
  filter(Binominal_name %in% x) %>% 
  filter(Strata == "Canopy" | Strata == "Understory")

df %>% select(Binominal_name) %>% unique %>% count()

df.mean  <- df %>% 
  group_by(Side, Binominal_name) %>% 
  summarise(across(Euclidian_distance:Shannon_tr, ~ mean(.x, na.rm = TRUE))) %>% 
  pivot_wider(names_from = Side, values_from = Euclidian_distance:Shannon_tr)
df.mean %>% head()

# read the morphological traits data
traits.raw <- read.csv(here::here("data/raw/traits_cont_bin.csv"), sep = ";")
traits.raw %>% head

df.mean <- df.mean %>% 
  left_join(traits.raw, by = join_by("Binominal_name" == "Species"))
df.mean %>% head

### pca hosts
## to perform the PCA for food specialization

host <- read.csv(here::here("data/raw/Taits_hostplants.csv"), sep = ";",
                 encoding = "latin1") %>% select(-X)
host%>% arrange(Species)

rownames(host) <- host$Binominal_names
host <- host[,-1]

pca.h <- prcomp(host, scale. = T)
summary(pca.h)
biplot(pca.h)
# larger positive values of PC1 indicate species that feed on many species of many families of plants species (generalists), while lower values indicate species that feed only in specific plants (specialists)

PC1 <- pca.h$x[, 1, drop = F] %>% as.data.frame()
PC1$Binominal_name <- rownames(PC1)

df.mean <- df.mean %>% left_join(PC1) %>% 
  dplyr::rename(Host = PC1)


# this data frame have information about morphological and ecophisiological traits for the studied species (34 spp) by wing side (ventral or dorsal). This distinction was made because color traits data was obtained measuring both sides of the wings.
# the traits values represent a mean for species, based in at least 3 individuals per species, and preferentially, males.

# selecting the columns of interest that will be used into analysis
# Iridescence_presence, Mean_lumiance, Eyespot_presence, Mean_contrast and Shannon_cd

df.mean %>% colnames()

df.mean <- df.mean %>% 
  select(Binominal_name, Mean_contrast_Dorsal, Mean_contrast_Ventral,
         Mean_luminance_Ventral, Mean_luminance_Dorsal,
         Iridescence_presence_Dorsal, Eyespot_presence_Ventral,
         Shannon_cd_Dorsal, Shannon_cd_Ventral,
         AR:FEA, Host)

## save the trait matrix for posterior analysis
write.xlsx(df.mean, here::here("data/processed/Mean_traits_bfly_new.xlsx"))

# save packages references
#writeLines(toBibtex(citation("PVR")), con = here::here("doc/Ref_PVR"))

# EXTRA - ploting traits along phylotree ----
# precisa terminar!!!

### PCA for color traits at each side
library(ggtree)

# ploting traits along the phylogeny

# scaling the continuous traits
df.mean %>% glimpse()

d1 <- df.mean %>%  
  select(-Iridescence_presence_Dorsal, -Eyespot_presence_Ventral) %>% 
  tibble::column_to_rownames('Binominal_name') %>% 
  scale()

d2 <- df.mean %>% select(Iridescence_presence_Dorsal, Eyespot_presence_Ventral, Binominal_name) %>% 
  tibble::column_to_rownames('Binominal_name')

new.labs <- phy.bfly$tip.label %>% as_tibble() %>% 
  mutate(Epithet = gsub(".*_", "", value)) %>% 
  mutate(spp = paste(substr(value, 1, 1), ". ", Epithet, sep = ""))

p1 <- ggtree(phy.bfly) + 
  geom_tiplab(offset = 20, align = TRUE) +
  ggnewscale::new_scale_fill()
p1

p2 <- gheatmap(p1, d1,
               offset = .5, width = .5,
               colnames_angle = 90, colnames_offset_y = .25, hjust = 1) +
  scale_fill_viridis_c(option = "H", name = "Standardized values") +
  ggnewscale::new_scale_fill() + theme(legend.position = "bottom")
p2

gheatmap(p1, d2, offset = .5, width = .2,
         colnames_angle = 90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option = "A", name = "discrete\nvalue")  

## Dorsal Side
# we will separet the sides and perform a PCA
tmp <- df.mean %>% ungroup() %>% 
  filter(Side != "Dorsal") %>% 
  select(Binominal_name, 
         Euclidian_distance:Iridescence_presence, 
         Eyespot_presence, Shannon_cd, Shannon_tr) %>% 
  left_join(data.color %>% select(Binominal_name, Tribe), multiple = "any") %>% 
  relocate(Tribe, .before = Euclidian_distance)

tmp %>% 
  ggplot() +
  geom_col(aes(y = reorder(Binominal_name, -Shannon_cd), x = Shannon_cd))

# scaling the continuous traits
d1 <- tmp %>% select(#Iridescence_presence,Eyespot_presence, 
  Mean_luminance, Mean_contrast, Shannon_cd, Binominal_name) %>%
  tibble::column_to_rownames('Binominal_name') %>% 
  scale()

d2 <- tmp %>% select(Iridescence_presence, Eyespot_presence, Binominal_name) %>% 
  tibble::column_to_rownames('Binominal_name')

p1 <- ggtree(phy.bfly, layout = "circular") + 
  geom_tiplab(offset = 20, align = TRUE) +
  ggnewscale::new_scale_fill()
p1

p2 <- gheatmap(p1, d1,
               offset = 0.5, width = .3,
               colnames_angle = 90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option = "H", name = "continuous\nvalue") +
  ggnewscale::new_scale_fill()
p2

gheatmap(p1, d2, offset = .5, width = .2,
         colnames_angle = 90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option = "A", name = "discrete\nvalue")  


