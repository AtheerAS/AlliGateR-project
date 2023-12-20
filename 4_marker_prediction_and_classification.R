#Marker prediction and classification using SVM


# Load required libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(tidyr)
library(proxy)
library(vegan)
library(corrplot)

# Load Seurat object (broad B cell annotation  and FACS-like annotation included) ----
ess = readRDS(file ="CITE_and_level1_anno.rds")

# Find new gatinbg markers using Seurat `FindAllMarkers` function ----
## Change identity and find all ADT markers within the FACS-like group ----
Idents(object = ess) <- "CITESeq_annot"
adt.markers_CITE <- FindAllMarkers(ess, assay = "ADT", only.pos = TRUE)

## Filter markers ----
adt.markers_CITE <- adt.markers_CITE %>%
  filter(avg_log2FC >= 2) %>%
  filter(p_val_adj < 0.05)

write.csv(adt.markers_CITE, quote = F, row.names = T, file = "Tables/dt_markers_based_on_FACS_label.csv")

# Use SVM classifier to rank markers ----
#Load and curate markers list
markers <- read.csv(file = "/Tables/dt_markers_based_on_FACS_label.csv")
markers$X <- gsub("-", "_", markers$X)
marker_list <- markers$X


# Create dataset from Seurat object 
dataset <- data.frame(CITESeq_annot = ess$CITESeq_annot,
                      Level_one_annotation=ess$Level_one_annotation)

# Get ADT assay data from Seurat object
assay_df <- as.data.frame(GetAssayData(ess,assay = "ADT", layer = "count"))
rownames(assay_df) <- gsub("-", "_", rownames(assay_df))


# SVM classifier functions

err_metric=function(CM){
  TN =CM[1,1]
  TP =CM[2,2]
  FP =CM[1,2]
  FN =CM[2,1]
  precision =(TP)/(TP+FP)
  recall_score =(FP)/(FP+TN)
  f1_score=2*((precision*recall_score)/(precision+recall_score))
  accuracy_model  =(TP+TN)/(TP+TN+FP+FN)
  False_positive_rate =(FP)/(FP+TN)
  False_negative_rate =(FN)/(FN+TP)
  sensitivity = (TP)/(TP+FN)
  print(paste("Precision value of the model: ",round(precision,2)))
  print(paste("Accuracy of the model: ",round(accuracy_model,2)))
  print(paste("Recall value of the model: ",round(recall_score,2)))
  print(paste("False Positive rate of the model: ",round(False_positive_rate,2)))
  print(paste("False Negative rate of the model: ",round(False_negative_rate,2)))
  print(paste("f1 score of the model: ",round(f1_score,2)))
  a = c(precision, accuracy_model, sensitivity, recall_score, False_positive_rate, False_negative_rate,f1_score )
  names(a) = c("precision", "accuracy_model", "sensitivity","recall_score", "False_positive_rate", "False_negative_rate","f1_score")
  return(a)
}

run_classifier<-function(value,factor){
  library(e1071)
  #Did you do the training with all values in the dataset?
  model <-svm(value,factor, type = "C-classification",kernel = "sigmoid")
  logit_P = predict(model , newdata = value ,type = 'response' )
  CM= table(factor, logit_P)
  mf = err_metric(CM)
  return(mf)
}

# Bind True values from both datasets (FACS-like popualtions and multi-omics populations)

data_list <- list(
  list(flow_like = "Naive", omics_label = "Naive"),
  list(flow_like = c("Switched memory PB", "Atypical PB") , omics_label = "Plasmablast"),
  list(flow_like = "Switched memory", omics_label = c("Memory (switched)")),
  list(flow_like = "Unswitched Memory" , omics_label = "Memory (unswitched)")
)

# Create output dataframe
column_names <- c("precision", "accuracy_model", "sensitivity",
                  "recall_score", "False_positive_rate",
                  "False_negative_rate", "f1_score", "cluster", "marker")

classifier_output_df <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
colnames(classifier_output_df) <- column_names

# Loop through datalist (true/false) values to run the classifire in each group.

for (celltype in data_list) {
  markers_F <- markers %>%
    filter(cluster %in% celltype$flow_like)
  
  for (marker in marker_list) {
    omics_marker_data <- assay_df[rownames(assay_df) == marker, ]
    if (nrow(omics_marker_data)>0){
      omics_marker_data <- t(omics_marker_data)
    }
    
    dataset_F <- dataset %>%
      filter(CITESeq_annot %in% celltype$flow_like) %>%
      mutate(value =  ifelse(Level_one_annotation %in% celltype$omics_label , "True", "False")) %>%
      merge(., omics_marker_data, by.x = "row.names", by.y = "row.names") %>%
      as.matrix()
    
    if (nrow(dataset_F)>0){
      
      
      type = dataset_F[,"value"]
      factor = factor(type)
      value = as.numeric(dataset_F[,ncol(dataset_F)])
      mf = run_classifier(value,factor)
      temp_df <- t(mf) %>%
        as.data.frame()
      temp_df$cluster <- paste0(celltype$flow_like, sep = "", collapse = " & ")
      temp_df$marker <- marker
      classifier_output_df <- rbind(classifier_output_df,temp_df)
    }
  }
}

write.csv(classifier_output_df, "Tables/marker_classification_sigmoid.csv", quote = F, row.names = F)

# density map ----

# Create data list with groups of interest

data_list <- list(
  list(flow_like = "Naive", omics_label = "Naive"),
  list(flow_like = "Switched memory", omics_label = c("Memory (switched)")),
  list(flow_like = "Unswitched Memory" , omics_label = "Memory (unswitched)")
)

#Create a character vector of markers of interest

marker_of_interest <- c("ADT_CD24","ADT_CD21","ADT_CD20")

#loop to create density map for each group/marker

for (celltype in data_list) {
  
  markers_F <- markers %>%
    filter(X %in% marker_of_interest) %>%
  filter(cluster == celltype$flow_like)
  marker_list <- markers_F$X
  
  for (i in marker_list) {
    adt_data <- assay_df[rownames(assay_df) == i, ]
    if(!is.data.frame(adt_data) || !nrow(adt_data)) next
    adt_data <- t(adt_data)
    
    dataset_F <- dataset %>%
      filter(CITESeq_annot == celltype$flow_like) %>%
      mutate(value =  ifelse(Level_one_annotation == celltype$omics_label, "True", "False")) %>%
      merge(., adt_data, by.x = "row.names", by.y = "row.names")
    
    
    ggplot(dataset_F, aes(x=dataset_F[,ncol(dataset_F)])) + 
      geom_density(aes(fill = value), linewidth = 0.5, alpha = 0.5) +
      scale_fill_brewer(palette = "Dark2", labels = c("Others", paste(celltype$flow_like)))+
      guides(fill=guide_legend(title="Multi-omics defined:")) +
      ggtitle(paste("FACS-defined:",celltype$flow_like))+
      xlab(paste(i," expression"))+
      theme(axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 18),
            plot.title = element_text(size = 18, hjust = 0.5))
    
    ggsave(filename = paste("markers/",celltype$flow_like,"_",i, ".png"),
           height = 2.5, width = 5, bg = "white")
    
  }
}








