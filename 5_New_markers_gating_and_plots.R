
# Load required libraries ----
library(Seurat)
library(CITEViz)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(tidyr)
library(tidyverse)

# Load Seurat object (broad B cell annotation  and FACS-like annotation included) ----
ess = readRDS(file ="CITE_and_level1_anno.rds")

ess_for_new_markers <- ess


# Original gating subsets ----

DefaultAssay(ess) <-"ADT"
atypnpb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`<4 & `adt_ADT-CD38`<8  )
atyppb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`<4 & `adt_ADT-CD38`>8 & `adt_ADT-CD24`<4 )
mempb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`<4 & `adt_ADT-IgM`<3 & `adt_ADT-CD38`>8 & `adt_ADT-CD24`<4)
unswmem<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`>4 & `adt_ADT-IgM`>3  & `adt_ADT-CD38`<8  )
memnpb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`<4 & `adt_ADT-IgM`<3 & `adt_ADT-CD38`<8  )
naive<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`>4)


# saveRDS(unswmem, "rds_files/unswmem.rds")
# saveRDS(memnpb, "rds_files/memnpb.rds")
# saveRDS(naive, "rds_files/naive.rds")


# Gate using CiteViz then load the data here ----

## Unswitched memory B cells - gated by new markers (CD21 and CD20)----

unsw_gated_by_21_and_20 <- readRDS(file = "rds_files/unsw_cd20_and_cd21_updated_18dec.rds" )
unsw_filtered <- subset(unswmem, subset = barcode_id %in% unsw_gated_by_21_and_20[["gate_1"]]@subset_cells[[1]])

## Switched memory B cells - gated by new markers (CD21 and CD20)----
sw_21_20 <- readRDS(file = "rds_files/memnpp_cd20_and_cd21_updated_18dec.rds")
mem_filtered <- subset(memnpb, subset = barcode_id %in% sw_21_20[["gate_1"]]@subset_cells[[1]])

## Naive B cells - gated by new markers (CD20 and CD24 gate)----
naive_cd20_24 <- readRDS(file = "rds_files/naive_cd20_cd24.rds")
naive_filtered <- subset(naive, subset = barcode_id %in% naive_cd20_24[["gate_1"]]@subset_cells[[1]])
#saveRDS(naive_filtered, "naive_filtered.rds")

## Naive B cells - gated by new markers (CD21 gate)----
naive_cd21 <- readRDS(file = "rds_files/naive_cd21.rds")
naive_cd21_from_all <- subset(naive_filtered, subset = barcode_id %in% naive_cd21[["gate_1"]]@subset_cells[[1]])

## Add gating labels to Seurat metadata ---- 
### Label the cells based on their populations ----

atyppb@meta.data$New_gate<-"Atypical PB"

atypnpb@meta.data$New_gate<-"Atypical non-PB"

mem_filtered@meta.data$New_gate<-"Switched memory"

mempb@meta.data$New_gate<-"Switched memory PB"

naive_cd21_from_all@meta.data$New_gate<-"Naive"

unsw_filtered@meta.data$New_gate<-"Unswitched Memory"

# function that merges the metadata of the subsets
annot_funct<-function(x)
{
  CellsName<-NULL
  for(i in 1: length(x))
  {
    CellsNameInt <- subset(x[[i]]@meta.data, select = c("New_gate"))
    if(is.null(CellsName)==F)
      CellsName <-rbind(CellsName,CellsNameInt) else
        CellsName <- CellsNameInt
  }
  return(CellsName)
}

# call the function using a vector that contains the subsets
metadatasubs<-annot_funct(c(mem_filtered,  atypnpb, unsw_filtered, atyppb, mempb,naive_cd21_from_all))
ess_for_new_markers <-AddMetaData(ess_for_new_markers, metadatasubs)

# New gating markers analysis ----
## Confusion matrix ----

confusion_matrix<-as.data.frame(table(ess_for_new_markers@meta.data$New_gate, ess_for_new_markers@meta.data$Level_one_annotation))
confusion_matrix$text.color <- ifelse(confusion_matrix$Freq > 300, "white", "black")

ggplot(data =  confusion_matrix, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1, color = confusion_matrix$text.color,
            size = 12) +
  scale_fill_gradient2(low =  "#9E9E9E" ,
                       mid =   "#F5C710",
                       high = "#DF536B",
                       space = "Lab",
                       na.value = "white",
                       guide = "colourbar",
                       trans="log") +
  xlab("FACS-defined groups") +
  ylab("Omics-defined groups")+
  theme_bw()+
  theme(axis.text=element_text(size=28),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_text(size = 30),
        legend.position="none")

ggsave("new_markers/new_gate_confusion_mat.png", height = 10, width = 15, bg = "white")


## Final stats (accuracy, specificity, and sensitivity) ----
confucsion_table <- as.data.frame(confusion_matrix)
confucsion_table$text.color <-NULL
colnames(confucsion_table) <- c("FACS-like_label","Multi-omics_label","Freq")

Table_of_analysis <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(Table_of_analysis) <- c("FACS-defined group","Set","Accuracy","Sensitivity","Specificity")

for (celltype in unique(confucsion_table$`FACS-like_label`)) {
  
  if (celltype == "Naive") {set <- "Naive"}
  else if (celltype == "Unswitched Memory") {set <- "Memory (unswitched)"}
  else if (celltype == "Atypical PB") {set <- "Plasmablast"}
  else if (celltype == "Switched memory PB") {set <- "Plasmablast"}
  else if (celltype == "Switched memory") {set <- c("Memory (switched)")}
  else {next}
  
  
  TP <-  as.numeric(confucsion_table %>%
                      filter(`FACS-like_label`== celltype, `Multi-omics_label` %in% set) %>%
                      select(Freq) %>%
                      sum())
  
  TN <- as.numeric(confucsion_table %>%
                     filter(`FACS-like_label` != celltype, `Multi-omics_label` != set) %>%
                     select(Freq) %>%
                     sum())
  
  FP <- as.numeric(confucsion_table %>%
                     filter(`FACS-like_label` == celltype, `Multi-omics_label` != set) %>%
                     select(Freq) %>%
                     sum())
  
  FN <- as.numeric(confucsion_table %>%
                     filter(`FACS-like_label` != celltype, `Multi-omics_label` %in% set) %>%
                     select(Freq) %>%
                     sum())
  
  Specificity <- round((TN/(TN+FP)*100),2)
  Sensitivity <- round((TP/(TP+FN)*100),2)
  Accuracy <- round(((TP+TN)/(TP+TN+FP+FN)*100),2)
  
  
  Table_of_analysis[nrow(Table_of_analysis)+1,] <- c(celltype, paste0(set, collapse = " and "), Accuracy,Sensitivity,Specificity )
  
}

write.csv(Table_of_analysis, "new_markers/New_gate_acc_sens_spec.csv")


# Between diseases analysis Using ggstatsplot ----

cell_comps_percentage <- data.frame(Level_one_annotation = ess_for_new_markers$Level_one_annotation,
                                    new_gate = ess_for_new_markers$New_gate,
                                    scRNASeq_sample_ID = ess_for_new_markers$scRNASeq_sample_ID,
                                    patient_group = ess_for_new_markers$Source_abrev)

Group3loop <- unique(na.omit(cell_comps_percentage$new_gate))
Group2loop <- unique(cell_comps_percentage$Level_one_annotation)

for (celltype in Group3loop) {
  for (omicstype in Group2loop) {
    try( 
      t2 <- cell_comps_percentage %>%
        filter(new_gate == celltype) %>%
        # select(Level_one_annotation, scRNASeq_sample_ID, ) %>%
        group_by(scRNASeq_sample_ID, Level_one_annotation) %>%
        summarize(Frequency = n()) %>%
        group_by(scRNASeq_sample_ID) %>%
        mutate(Normalized_Frequency = Frequency / sum(Frequency)*100) %>%
        left_join(cell_comps_percentage %>% select(scRNASeq_sample_ID, patient_group) %>% distinct(), 
                  by = "scRNASeq_sample_ID") %>%
        ungroup() %>%
        filter(Level_one_annotation == omicstype)
    )
    if (sum(t2$Normalized_Frequency) == 0) { next }
    
    
    try(
      
      ggbetweenstats(data = t2,
                     x = patient_group,
                     y = Normalized_Frequency,
                     type = "parametric", # ANOVA or Kruskal-Wallis
                     var.equal = TRUE, # ANOVA or Welch ANOVA
                     plot.type = "box",
                     pairwise.comparisons = TRUE,
                     pairwise.display = "significant",
                     centrality.plotting = FALSE,
                     violin.args = list(width = 0),
                     bf.message = FALSE,
                     ylab = paste0("Multi-omics defined: ", omicstype, " %"),
                     xlab = "Patient groups",
                     title = paste0("FACS-defined: ", celltype),
                     ggsignif.args =list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                     ggplot.component = list(
                       theme(axis.text=element_text(size=13),
                             axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
                             axis.title = element_text(size = 15),
                             plot.title = element_text(hjust = 0.5, size = 15 ))))
    )
    
    ggsave(filename = paste0("new_markers/between_patients_group_analysis2/",omicstype, "_",celltype,".png"), 
           height = 5.5, width = 5.5, bg = "white")
    
  }
}




## Changes in purity before and after ----

confusion_matrix_new_gate<-as.data.frame(table(ess_for_new_markers@meta.data$New_gate, ess_for_new_markers@meta.data$Level_one_annotation))
confusion_matrix_old_gate<-as.data.frame(table(ess_for_new_markers@meta.data$CITESeq_annot, ess_for_new_markers@meta.data$Level_one_annotation))

confusion_matrix_new_gate$old_gate_freq <- confusion_matrix_old_gate$Freq
confusion_matrix_new_gate$differences <- ((confusion_matrix_old_gate$Freq)-(confusion_matrix_new_gate$Freq))/((confusion_matrix_old_gate$Freq)+(confusion_matrix_new_gate$Freq))

mat <- confusion_matrix_new_gate %>%
  select(Var1, Var2, differences) %>%
  pivot_wider(
    names_from = Var1,
    values_from = differences
  ) %>%
  column_to_rownames(var="Var2")

write.csv(mat, "new_markers/purity_matrix.csv")







