#Gating FACS-like populations

# Load required libraries ----

library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
library(ggrepel)
library(ggtext)
library(cowplot)
library(tidyverse)
library(ggh4x)
library(grid)
library(ggstatsplot)
library(CITEViz)

# Load Seurat Object ----
ess = readRDS(file = "20210623_1_Final_All_bcells_Combat.rds")

# Add multi-omics broad B cell annotation ----
ess@meta.data <- ess@meta.data %>%
  mutate(Level_one_annotation = case_when(
    name.unique %in% c("B.cyc","B.int.1.early.act/sw","B.int.2.early.act.IFN.resp","B.int.2.early.act/sw","B.int.2.IFN.resp","B.mitohi") ~ "Activated (switched)",
    name.unique %in% c("B.int.1.early.act","B.int.1.IFN.resp","B.int.2.unsw") ~ "Activated (unswitched)",
    name.unique %in% c("B.NAIVE","B.NAIVE.CD1c","B.NAIVE.IFN.resp","B.NAIVE.IgDlo") ~ "Naive",
    name.unique %in% c("B.SW.MEM","B.SW.MEM.IFN.resp")~ "Memory (switched)",
    name.unique %in% c("B.TRANSIT.CD10") ~"Transitional",
    name.unique %in% c("B.UNSW.MEM") ~ "Memory (unswitched)",
    name.unique %in% c("PB","PB.cyc","PB.IFN.resp","PB.mitohi") ~ "Plasmablast" ))



# Gate populations using ADT markers ----
DefaultAssay(ess) <-"ADT"
atypnpb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`<4 & `adt_ADT-CD38`<8  )
atyppb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`<4 & `adt_ADT-CD38`>8 & `adt_ADT-CD24`<4 )
mempb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`<4 & `adt_ADT-IgM`<3 & `adt_ADT-CD38`>8 & `adt_ADT-CD24`<4)
unswmem<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`>4 & `adt_ADT-IgM`>3  & `adt_ADT-CD38`<8  )
memnpb<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`>5 & `adt_ADT-IgD`<4 & `adt_ADT-IgM`<3 & `adt_ADT-CD38`<8  )
naive<-subset(ess, subset = `adt_ADT-CD10`<4 & `adt_ADT-CD45`>0 & `adt_ADT-CD3`<4 & `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`>4)

# Add gating labels to Seurat metadata ---- 
## Label the cells based on their populations ----

atyppb@meta.data$CITESeq_annot<-"Atypical PB"

atypnpb@meta.data$CITESeq_annot<-"Atypical non-PB"

memnpb@meta.data$CITESeq_annot<-"Switched memory"

mempb@meta.data$CITESeq_annot<-"Memory PB"

naive@meta.data$CITESeq_annot<-"Naive"

unswmem@meta.data$CITESeq_annot<-"Unswitched Memory"

## function that merges the metadata of the subsets ----
annot_funct<-function(x)
{
  CellsName<-NULL
  for(i in 1: length(x))
  {
    CellsNameInt <- subset(x[[i]]@meta.data, select = c("CITESeq_annot"))
    if(is.null(CellsName)==F)
      CellsName <-rbind(CellsName,CellsNameInt) else
        CellsName <- CellsNameInt
  }
  return(CellsName)
}

## call the function using a vector that contains the subsets ----
metadatasubs<-annot_funct(c(naive, atyppb, memnpb, mempb, atypnpb, unswmem))
ess <-AddMetaData(ess, metadatasubs)










