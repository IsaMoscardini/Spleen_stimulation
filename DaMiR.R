######
rm(list = ls())
options(stringsAsFactors = F)
#setwd("PROJECTS/Spleen/")
library(S4Vectors)
library(rJava)
library(survival)
library(DaMiRseq)
library(lattice)
library(RColorBrewer)

DaMiR.Clustplot <- function (data, df, type_row = c("euclidean", "correlation"), 
          type_col = c("euclidean", "correlation")) {
  if (missing(data)) 
    stop("'data' argument must be provided")
  if (missing(df)) 
    stop("'df' argument must be provided")
  if (missing(type_row)) {
    type_row <- type_row[1]
  }
  if (missing(type_col)) {
    type_col <- type_col[1]
  }
  if (!(is(data, "SummarizedExperiment") | is.data.frame(data) | 
        is.matrix(data))) 
    stop("'data' must be a 'matrix', a 'data.frame'\n         or a 'SummarizedExperiment' object")
  if (!(is(df, "DataFrame") | is.data.frame(df))) 
    stop("'df' must be a data.frame")
  if (is(data, "SummarizedExperiment")) {
    count_data <- assay(data)
  }
  else {
    count_data <- t(data)
  }
  df <- as.data.frame(df)
  if (any(is.na(count_data))) 
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df))) 
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(count_data))) 
    stop("Inf values are not allowed in the 'data' matrix")
  if (!("class" %in% colnames(df))) 
    stop("'class' info is lacking!\n         Include the variable 'class'\n         in the 'df' data frame and label it 'class'!")
  if (all((count_data%%1) == 0)) 
    stop("This function works with normalized count data")
  if (dim(count_data)[2] != dim(df)[1]) 
    stop("ncol(assay(data)) must be equal to nrow(df)")
  if (type_row == "euclidean") {
    d_r <- "euclidean"
  }
  else if (type_row == "correlation") {
    d_r <- "correlation"
  }
  else {
    stop("Please set 'euclidean' or 'correlation' as distance type.")
  }
  if (type_col == "euclidean") {
    d_c <- "euclidean"
  }
  else if (type_col == "correlation") {
    d_c <- "correlation"
  }
  else {
    stop("Please set 'euclidean' or 'correlation' as distance type.")
  }
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
  pheatmap(count_data, clustering_distance_rows = d_r, clustering_distance_cols = d_c, 
           scale = "row", col = colors, annotation_col = df)
}

#######
countdata <- read.csv2('Data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('Data/Mouse_desc_20191024.csv')
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))
#View(descript)
descript$Class <- paste0(descript$Group, "_", descript$Time.point)
rownames(descript) <- paste0("X", descript$Sample.Name)
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))
identical(colnames(countdata), rownames(descript))

descript <- descript[descript$Sample.Name != "1_6_NS" ,]
#descript <- descript[descript$Sample.Name != "1_3_NS" ,]

descript <- descript[descript$Vitro == "T", ]
common.samples <- intersect(colnames(countdata), rownames(descript))
countdata <- countdata[,common.samples]

minCount = 1
high <- countdata > minCount
idx <- (rowSums(high)/ncol(countdata)) == 1 # expressed in all samples
countdata <- countdata[which(idx),]
nrow(countdata)

######## pre process ########

colnames(descript) <- c("Barcode", "Sample", "Mapped_reads", "Valid_reads", "targets_detected", "Chip", "Vitro", "Vivo", "TP", "group", "class")
#View(descript)
descript[] <- lapply(descript[], as.factor)
descript$Valid_reads <- gsub("%", "", descript$Valid_reads)
descript$targets_detected <- gsub("%", "", descript$targets_detected)
descript <- as.data.frame(descript)
descript$class <- as.factor(descript$class)
descript$Valid_reads <- as.numeric(descript$Valid_reads)
descript$targets_detected <- as.numeric(descript$targets_detected)

str(descript)

descript$class <- gsub("I_T", "I_S", descript$class) 
 
View(descript)

descript[] <- lapply(descript[], as.factor)
str(descript)

countdata <- countdata[, rownames(descript)] 
identical(rownames(descript), colnames(countdata))

descript$Vitro <- NULL
descript$Vivo <- NULL
View(descript)

descript$class <- gsub("I_NS_1", "I_NS_D12", descript$class)
descript$class <- gsub("I_NS_2", "I_NS_D12", descript$class)
descript$class <- gsub("I_NS_4", "I_NS_D47", descript$class)
descript$class <- gsub("I_NS_7", "I_NS_D47", descript$class)

descript$class <- gsub("I_S_1", "I_S_D12", descript$class)
descript$class <- gsub("I_S_2", "I_S_D12", descript$class)
descript$class <- gsub("I_S_4", "I_S_D47", descript$class)
descript$class <- gsub("I_S_7", "I_S_D47", descript$class)

descript$class <- as.factor(descript$class)

#####
# Criar Summarized Experiment object
SE <- DaMiR.makeSE(x = countdata, y = descript) # 23930 Features

# Filtering and applying VST normalization
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "no")    

# Quality check on the correlation of samples
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.8)  

sv <- DaMiR.SV(data_filt)
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)
#pdf("Figures/Quality_Check_DaMiRseq.pdf")
#DaMiR.Allplot(data_filt, colData(data_filt))
#dev.off()

sv <- sv[,c(1,2,4,5,6,7,8)]
data_adjust <- DaMiR.SVadjust(data_filt, sv, n.sv=7)
assay(data_adjust[c(1:5), c(1:5, 21:25)])

pdf("Figures/Paper/Quality_Check_DaMiRseq_pos_QC.pdf")
DaMiR.Allplot(data_adjust, colData(data_adjust))
dev.off()

set.seed(1234567)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)

data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.4)
genes_damir <- colnames(data_reduced$data)

data_reduced <- DaMiR.FReduct(data_reduced$data) # 32 remained

pdf("Figures/Paper/Reduced_25genes.pdf", width = 6, height = 4)
DaMiR.MDSplot(data_reduced, df)
dev.off()

#colnames(data_reduced) <- paste0("a_", colnames(data_reduced)) 
pdf("Figures/Paper/Importance_DaMiR.pdf", width = 8, height = 8)
df.importance <- DaMiR.FSort(data_reduced, df)
dev.off()
View(df.importance)

selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,
                                 n.pred = 11)
df2 <- data.frame(df)
df2$Barcode <- NULL
df2$Sample <- NULL
df2$Mapped_reads <- NULL
df2$Valid_reads <- NULL
df2$targets_detected <- NULL
df2$Chip <- NULL
df2$TP <- NULL
df2$group <- NULL
selected <- selected_features$data
names <- rownames(selected)
names <- gsub("T", "S", names)
names <- gsub("X", "", names)
DaMiR.Clustplot(selected_features$data, df2)
