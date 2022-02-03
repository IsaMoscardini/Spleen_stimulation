######

rm(list = ls())
options(stringsAsFactors = F)
pkgs <- c('dplyr','data.table','reshape2','tidyr','tidyverse', 'ggpubr')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)


###### Pre Process ####
cytok <- read.csv("Data/Bioplex dati elaborati (3)  - Raw data (1).csv")
class(cytok$X)
cytok <- cytok[cytok$X %in% c("1^1", "1^2", "1^3", "1^4", "1^5", "1^6", 
                              "4^1", "4^2", "4^3", 
                              "6^1", "6^2", "6^3", 
                              "8^1", "8^2", "8^3", "8^4", "8^5", "8^6", 
                              "10^1", "10^2", "10^3", "10^4", "10^5", "10^6"),]
View(cytok)

## Select only TIGR4 
cytok <- cytok[cytok$X.2 %in% c("TIGR4", "NS"),]
cytok$X <- gsub("\\^", "_", cytok$X)
cytok$X.2 <- gsub(pattern = "TIGR4", replacement = "T", x = cytok$X.2) 
cytok$samples <- paste(cytok$X, "_", cytok$X.2)
cytok$samples <- gsub(" ", "", cytok$samples)
View(cytok)

pheno <- read.csv2("Data/Mouse_desc_20191024.csv")
cytok$TP <- pheno$Time.point[match(cytok$samples, pheno$Sample.Name)]

rownames(cytok) <- cytok$samples

# Select cytokine names
names_cyto <- colnames(cytok)[4:26]

cytok_values <- apply(cytok, 2, gsub, pattern = "\\*", replacement = "")
cytok_values <- cytok[,c(4:26, 28)]

# Check if everything is numeric
cytok_values <- apply(cytok_values, 2, as.numeric)
rownames(cytok_values) <- paste0("X", cytok$X)

# Check values and min
Cyto_d1 <- cytok %>% filter(TP %in% c(0, 1, 2, 4, 7)) %>% gather(cytokine, value, -c(X, X.1, X.2, samples, TP)) %>%
  mutate(value = as.numeric(gsub("\\*", "", .$value)), 
         TP = as.factor(TP))
min(Cyto_d1$value)

cyto_clean <- Cyto_d1 %>% spread(cytokine, value)

boxplot(cyto_clean[, names_cyto], las=2)



####### Boxplots #####

Cyto_d1$X.2 <- gsub("T", "Stimulated", Cyto_d1$X.2)
Cyto_d1$X.2 <- gsub("NS", "Non_stimulated", Cyto_d1$X.2)


## With the 6 main cytokines

Cyto_d1 %>%
  filter((cytokine %in% c("Il17a","Gm.csf","MIP_1a","IFN.g","KC","Il2"))) %>%
  ggplot() + aes(x = TP, y = value, fill = X.2) + 
  geom_boxplot() +
  facet_wrap(vars(cytokine), scales = "free") +
  scale_fill_manual(values = c(Non_stimulated = "#B5FFF0", Stimulated = "#EABFFF")) +
  #scale_color_manual(values = c(Stimulated = "#EABFFF",Unstimulated = "#B5FFF0")) +
  theme_minimal() +
  theme(legend.position = "top") + labs(x = "Time Point (Days)", y= "Concentration (pn/mL)")


## All cytokines - Panel

Cyto_d1 %>%
  ggplot() + aes(x = TP, y = value, fill = X.2) + 
  geom_boxplot() +
  facet_wrap(vars(cytokine), scales = "free") +
  scale_fill_manual(values = c(Non_stimulated = "#B5FFF0", Stimulated = "#EABFFF")) +
  #scale_color_manual(values = c(Stimulated = "#EABFFF",Unstimulated = "#B5FFF0")) +
  theme_minimal() +
  theme(legend.position = "top") + labs(x = "Time Point (Days)", y= "Concentration (pn/mL)")













