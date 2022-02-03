######

rm(list = ls())
options(stringsAsFactors = F)
pkgs <- c('dplyr','data.table','reshape2','tidyr','tidyverse', 'ggpubr','mixOmics','DESeq2')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

"%ni%" <- Negate("%in%")


##### PRE PROCESS ####
# RNA
countdata <- read.csv2('Data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('Data/Mouse_desc_20191024.csv')
#View(countdata)
descript$Sample.Name <- gsub("T", "S", descript$Sample.Name)
descript$Vitro <- gsub("T", "S", descript$Vitro)
colnames(countdata) <- gsub("T", "S", colnames(countdata))
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))
# Check samples and row/column names
rownames(descript) <- paste0("X", descript$Sample.Name)
common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))
identical(colnames(countdata), rownames(descript))

# Remove outliers from descript table
remove <- c("X1_6_NS") 
descript_clean <- descript[setdiff(rownames(descript), remove),]
nrow(descript_clean)

setdiff(colnames(countdata), rownames(descript_clean)) # "X1_6_NS"
setdiff(rownames(descript_clean), colnames(countdata)) # 0

# Remove outliers from countdata table
common.samples <- intersect(colnames(countdata), rownames(descript_clean))
countdata_clean <- countdata[, common.samples]
identical(colnames(countdata_clean), rownames(descript_clean))

Group <- descript_clean$Group
dds <- DESeqDataSetFromMatrix(countdata_clean, descript_clean, design = ~Group)
dds$Group <- relevel(dds$Group, "NI_NS" )

minCount = 1
high <- assay(dds) > minCount
idx <- (rowSums(high)/ncol(dds)) == 1 # expressed in all samples
dds <- dds[which(idx),]
nrow(dds) # 10757

countdata_clean <- counts(dds)
countdata_clean <- as.data.frame(countdata_clean)

# Normalize and transform
rnames <- rownames(countdata_clean)
countdata_clean <- varianceStabilizingTransformation(as.matrix(countdata_clean))
rownames(countdata_clean) <- rnames

countdata_clean <- as.data.frame(countdata_clean)
countdata_clean$genes <- rownames(countdata_clean)


#### CYTOKINES
cytok <- read.csv("Data/Bioplex dati elaborati (3)  - Raw data (1).csv")
class(cytok$X)
cytok <- cytok[cytok$X %in% c("1^1", "1^2", "1^3", "1^4", "1^5", "1^6", 
                              "4^1", "4^2", "4^3", "6^1", "6^2", "6^3", 
                              "8^1", "8^2", "8^3", "8^4", "8^5", "8^6", 
                              "10^1", "10^2", "10^3", "10^4", "10^5", "10^6"),]

## TIGR4 
cytok <- cytok[cytok$X.2 %in% c("TIGR4", "NS"),]

cytok$X <- gsub("\\^", "_", cytok$X)
cytok$X.2 <- gsub(pattern = "TIGR4", replacement = "T", x = cytok$X.2) 
cytok$samples <- paste(cytok$X, "_", cytok$X.2)
cytok$samples <- gsub(" ", "", cytok$samples)
cytok$samples <- gsub("T", "S", cytok$samples)
cytok$X.2 <- gsub("T", "S", cytok$X.2)

pheno <- descript_clean
cytok$TP <- pheno$Time.point[match(cytok$samples, pheno$Sample.Name)]
cytok$Group <- pheno$Group[match(cytok$samples, pheno$Sample.Name)]
rownames(cytok) <- cytok$samples

cytok_values <- apply(cytok, 2, gsub, pattern = "\\*", replacement = "")
rownames(cytok_values) <- paste0("X",rownames(cytok_values))

##### Bind tables
cp <- countdata_clean[, intersect(rownames(cytok_values), colnames(countdata_clean))]
rna <- t(cp)

cytok_values <- cytok_values[intersect(rownames(cytok_values), colnames(cp)),]
identical(rownames(rna), rownames(cytok_values)) # TRUE
cytok_values <- as.data.frame(cytok_values)
cytok_values$X <- NULL
cytok_values$X.1 <- NULL
cytok_values$samples <- NULL
#write.csv(cytok_values, "Data/Cytok_filtered.csv")

# All numeric
cytok_values[, 2:24] <- apply(cytok_values[, 2:24], 2, as.numeric)
class(cytok_values$Eotaxin) # "numeric"

pheno <- cytok_values[, c(25,26)]
cytok_values$X.2 <- NULL
cytok_values$Group <- NULL
cytok_values$TP <- NULL
cytok_values <- as.matrix(cytok_values)
rna <- as.matrix(rna)

cytok_values <- log10(cytok_values+1)
colnames(cytok_values) <- gsub('\\.','_',gsub('\\.\\..*','',colnames(cytok_values)))
identical(rownames(rna), rownames(cytok_values)) # TRUE
identical(rownames(rna), rownames(pheno)) # TRUE

# Check cytokines
boxplot(cytok_values, las = 2)
min(cytok_values) # 0
head(pheno)

# Check counts
counts.melt = rna %>% as.data.frame %>% rownames_to_column() %>% 
  gather(gene, value, -rowname) %>% 
  merge(pheno, by.y = 'row.names', by.x = 'rowname') 

counts.melt %>% ggplot(aes(value, group = rowname, fill = Group)) + 
  geom_density(alpha = .6) + theme_linedraw() + theme(legend.position = 'top') +
  facet_grid(TP ~., scales = 'free') + labs(x= 'Expression Value')

# Check distributiion
pca <- pca(X = cytok_values, scale = T, center = T) 
plotIndiv(pca, group = pheno$Group)



#### MixOmics PLS and sPLS #####
mRNA <- rna
colnames(cytok_values) <- toupper(colnames(cytok_values))
Cytokines <- cytok_values

MyResult.pls <- pls(mRNA,Cytokines, ncomp = 4)  
perf.pls <- perf(MyResult.pls, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975) # 2 components
perf.pls$Q2.total
# Q2.total
# 1 comp  0.33030771 // 2 comp  0.06093918 // 3 comp -0.01713423 // 4 comp -0.02458433

## Tuning based on MAE
list.keepX <- c(5:20, 70:90)
tune.spls.MAE <- tune.spls(mRNA, Cytokines, ncomp = 2,
                           test.keepX = list.keepX,
                           validation = "Mfold", folds = 5,
                           nrepeat = 100, progressBar = TRUE,
                           measure = 'MAE')

plot(tune.spls.MAE, legend.position = 'topright')
tune.spls.MAE$choice.keepX

pls <- pls(mRNA,Cytokines)  
plotIndiv(pls, group = pheno$Group, ind.names = pheno$TP, legend = T) ## sample plot                   

## sPLS
MyResult.spls <- spls(mRNA,Cytokines, keepX = c(16, 25))
plotIndiv(MyResult.spls, group = pheno$Group, ind.names = pheno$TP, legend = T) ## sample plot                   

plotIndiv(MyResult.spls, group = pheno$Group,
          rep.space = "XY-variate", legend = TRUE,
          legend.title = 'Group',
          ind.names = pheno$TP,
          title = 'sPLS', col.per.group = c("#33CC33", "#CC0099", "#00CCCC", "#CCCC00"))

plotVar(MyResult.spls, cex=c(3,3), legend = TRUE, col = c("#0066CC", "#FF0000"))

X11()
cim(MyResult.spls, comp = 1)
cim(MyResult.spls, comp = 2)



#### Test with 25 features ####
MyResult.spls <- spls(mRNA, Cytokines, keepX = c(25, 25))
plotIndiv(MyResult.spls, group = pheno$Group, ind.names = pheno$TP, legend = T) ## sample plot                   

plotIndiv(MyResult.spls, group = pheno$Group,
          rep.space = "XY-variate", legend = TRUE,
          legend.title = 'Group',
          ind.names = pheno$TP,
          title = 'sPLS', col.per.group = c("#33CC33", "#CC0099", "#00CCCC", "#CCCC00"))

plotVar(MyResult.spls, cex=c(3,3), legend = TRUE, col = c("#0066CC", "#FF0000"))

View(MyResult.spls$loadings$X)
x11()
cim(MyResult.spls, comp = 1)
cim(MyResult.spls, comp = 2)



