library(vegan)
library(readxl)
library(tidyverse)
library(mvabund)
library(parallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- community heritability calculated on arthropod community ----

# ==== Traditional approach ====
# 
# This approach presents several problems at several stages. They will be discussed in the paper (choice of the distance matrix, the ordination method and so forth)

arth_df <- read_excel("../data/Lamit_CWGenotypeForHeritabilityMethods.xlsx", sheet = 2)

arth_matrix <- arth_df %>% 
  select(-2) %>% 
    column_to_rownames("Tree")

arth_NMDS=metaMDS(arth_matrix,k=2,trymax=5000)

(arth_NMDS_scs <- data.frame(arth_NMDS$points) %>% 
  select(-2) %>% 
  rownames_to_column(var="Tree") %>% 
    left_join(select(arth_df, 1:2)) %>% 
    mutate(Genotype=as.factor(Genotype)))


# another thing to consider is that this method of calculating heritability would also require checking the assumptions of linear models
arth_lm <- lm(MDS1~Genotype, data=arth_NMDS_scs)

plot(arth_lm)

summary(arth_lm) # Heritability is nearly 50%

anova(arth_lm) # and it is significant at 0.05

# ====  mvabund approach ====

arth_mvb <- mvabund(arth_matrix)

# model with Genotype. Genotype is significant. Heritability is 20% compared to 47-48%%.  

arth_MGLM <- manyglm(arth_mvb~ Genotype, data = arth_df)

plot(arth_MGLM)

# for significance of heritability
arth_mvb_out <- anova(arth_MGLM, nBoot = 999)

# Heritability calculation

arth_int <- manyglm(arth_mvb~ 1) # intercept only model

(null_dev_arth <- -2*sum(logLik(arth_int)))

res_dev_arth <- -2*sum(logLik(arth_MGLM))

((null_dev_arth-res_dev_arth)/null_dev_arth)*100


# ---- Endophytes ----

# ==== Traditional approach ====

ends_df <- read_excel("../data/Lamit_CWGenotypeForHeritabilityMethods.xlsx", sheet = 4)

ends_matrix <- ends_df %>% 
  select(-2) %>% 
  column_to_rownames("Tree")

ends_NMDS=metaMDS(ends_matrix,k=2,trymax=5000)

(ends_NMDS_scs <- data.frame(ends_NMDS$points) %>% 
    select(-2) %>% 
    rownames_to_column(var="Tree") %>% 
    left_join(select(ends_df, 1:2)) %>% 
    mutate(Genotype=as.factor(Genotype)))


# another thing to consider is that this method of calculating heritability would also require checking the assumptions of linear models

ends_lm <- lm(MDS1~Genotype, data=ends_NMDS_scs)

plot(ends_lm)

summary(ends_lm) # heritability is about 18%

anova(ends_lm) # heritability is non significant


# ====  mvabund approach ====

ends_mvb <- mvabund(ends_matrix)

# model with Genotype. Genotype is significant. Heritability is 20% compared to 47-48%%.  

ends_MGLM <- manyglm(ends_mvb~ Genotype, data = ends_df)

plot(ends_MGLM)

# for significance of heritability

(ends_mvb_out <- anova(ends_MGLM, nBoot = 999, show.time = "all"))

# Heritability calculation

ends_int <- manyglm(ends_mvb~ 1) # intercept only model

(null_dev_ends <- -2*sum(logLik(ends_int)))

res_dev_ends <- -2*sum(logLik(ends_MGLM))

((null_dev_ends-res_dev_ends)/null_dev_ends)*100


# ---- Microbiome ---- 

# ==== Traditional approach ====

mcrb_df <- read_csv("../data/otu_table.csv")

mcrb_matrix <- mcrb_df %>% 
  column_to_rownames(var="Sample_ID") %>% 
  select(-Accession)

mcrb_NMDS=metaMDS(mcrb_matrix,k=2,trymax=5000)

(mcrb_NMDS_scs <- data.frame(mcrb_NMDS$points) %>% 
    select(-2) %>% 
    rownames_to_column(var="Sample_ID") %>% 
    left_join(select(mcrb_df,Sample_ID,Accession)))

# community heritability is low and non significant

mcrb_lm <- lm(MDS1~Accession, data=mcrb_NMDS_scs)

plot(mcrb_lm )

summary(mcrb_lm )

anova(mcrb_lm)

# ====  mvabund approach ====

mcrb_mvb <- mvabund(mcrb_matrix)

# model with Genotype. Genotype is significant. Heritability is 20% compared to 47%

mcrb_MGLM <- manyglm(mcrb_mvb ~ Accession, data = mcrb_df)

plot(mcrb_MGLM)

# for significance of heritability
mcrb_mvb_out <- mclapply (1:7, function(x) {

anova(mcrb_MGLM, nBoot = 30, show.time = "all")
    
}, mc.cores = 7)
  
K=7 # cores
N=30  #bootstrapping

pi_acc <- map_dbl(mcrb_mvb_out, ~ .x$table[2,4])
piMean_acc = mean(pi_acc)
(p_acc = piMean_acc + (piMean_acc-1)*(K-1)/(K*N+1)) # Accession is significant

# Heritability calculation

mcrb_int <- manyglm(mcrb_mvb~ 1) # intercept only model

(null_dev_mcrb <- -2*sum(logLik(mcrb_int)))

res_dev_mcrb <- -2*sum(logLik(mcrb_MGLM))

((null_dev_mcrb-res_dev_mcrb)/null_dev_mcrb)*100


