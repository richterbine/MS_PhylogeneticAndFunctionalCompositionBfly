############################################################################
##### MS phylo and functional composition of fruit-feeding butterflies #####
############################################################################

# Script to read raw data, manipulate the phylogenetic tree and alometric relationship among traits --------
# packages required
library(dplyr)
library(reshape)

# read the raw data
# species information
data.raw <- read.csv(here::here("data/raw/Raw_data_FLONA.csv"), sep = ";")
data.raw %>% head

# remove Temenis_laothoe
data.raw <- data.raw[-which(data.raw$Species == "Temenis_laothoe"),]

write.csv(data.raw, here::here("data/processed/data.raw.csv"))

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
write.csv(data.bfly, here::here("data/processed/data.bfly.csv"))

######################################################################
############ Phylogenetic manipulation - pruning tree ################
######################################################################

# Loading the Nymphalidae tree (Chazot et al. 2019)

source(here::here("R", "functions", "genus_to_spp_function.R"))

# Manipulating the proposed phylogenetic hypothesis
Nym.tree <- ape::read.tree(here::here("data/raw/Nymphalidae_tree.txt"))

# Cutting the Nym.tree to my specific sampling pool
tree.bfly <- genus_spp_tree(tree = Nym.tree, 
                            genus.species.list = data.raw[, c("Genus", "Species")])
tree.bfly

ape::write.tree(tree.bfly, here::here("data/processed/tree_bfly_flona.txt"))


######################################################################
################ Functional manipulation procedures ##################
######################################################################

# read the data
traits.raw <- read.csv(here::here("data/raw/traits_cont_bin.csv"), sep = ";")

### pca hosts
## to perform the PCA for food specialization

host <- read.csv(here::here("data/raw/Taits_hostplants.csv"), sep = ";")
rownames(host) <- host$X
host <- host[,-1]

pca.h <- prcomp(host)
summary(pca.h)
PC1 <- pca.h$x[, 1, drop = F]

host$PC1 <- PC1

traits.raw$Host <- PC1

## save the trait matrix for posterior analysis
write.csv(traits.raw, here::here("data/processed/Mean_traits_bfly.csv"))
