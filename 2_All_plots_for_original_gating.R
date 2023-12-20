# All generated plots

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
library(scales)
library(scCustomize)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(RColorBrewer)


# Load Seurat object (broad B cell annotation and FACS-like annotation included) ----
ess = readRDS(file ="CITE_and_level1_anno.rds")

# UMAP ----

## colors ----
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# Define a palette of 21 colors
my_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
                "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
                "#ffff99", "#b15928", "#8dd3c7", "#fb8072", "#80b1d3",
                "#ccebc5", "#ffed6f", "#d9d9d9", "#bc80bd", "#ccebc5",
                "#80cdc1")

## UMAP broad B cell annotation ----
DimPlot(ess, reduction = "umap_allBcells",label = F, label.size = 7, pt.size = 1, label.color = "black",
        group.by = "pseudobulk",
        repel = T) +
  ggtitle(label = "") +
  scale_color_manual("Pseudobulk Annotation",values = my_palette)+
  guides(guide_legend(ncol = 3))

ggsave("umap_pseudobulk.png", height = 5, width = 7, bg = "white")

## UMAP FACS-like annotation ----
DimPlot(ess, reduction = "umap_allBcells",label = F, label.size = 7, pt.size = 1, label.color = "black"
        ,group.by = "name.unique",
        repel = T) +
  ggtitle(label = "") +
  scale_color_manual("Detailed Annotation",values = my_palette)+
  guides(guide_legend(ncol = 3))

ggsave("umap_detailed_annotation.png", height = 6, width = 9, bg = "white")

## UMAP detailed annotation ----

DimPlot(ess, reduction = "umap_allBcells",label = F, label.size = 7, pt.size = 1, label.color = "black", 
        group.by = "Level_one_annotation",
        repel=TRUE) +
  ggtitle(label = "") +
  scale_color_manual("Level One Annotation", values = my_palette)


# Confusion matrix ----
## All ----

confusion_matrix<-as.data.frame(table(ess@meta.data$CITESeq_annot, ess@meta.data$Level_one_annotation))
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

ggsave("confusion_mat_F.png", height = 10, width = 15, bg = "white")

## Confusion matrix for: Atypical + Naive + and Switched memory ----
### multi-omics defined: Memory (switched) ----

t <- confusion_matrix %>%
  filter(Var2 == "Memory (switched)")

ggplot(data = t, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = paste(gsub("\\s+", "\n", Var1), (sprintf("%1.0f", Freq)), sep ='\n')), vjust = 0.4, size=12,
            color = t$text.color ) +
  scale_fill_gradient( low = "#F5C710",
                       high = "#DF536B",
                       space = "Lab",
                       na.value = "white",
                       trans="log",
                       guide = "none") +
  xlab("")+
  ylab("FACS-defined:")+
  ggtitle("Multi-omics group: Switched memory B cells") +
  theme_bw()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  coord_fixed(ratio = 1)+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size = 29),
        plot.title = element_text(size = 29))
ggsave("confusion_mat_memory_sw_b_cell.png", height = 4, width = 15, bg = "white")


### FACS-like defined: Atypical non-PB ----

t <- confusion_matrix %>%
  filter(Var1 == "Atypical non-PB")

ggplot(data = t, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = paste(gsub("\\s+", "\n", Var2), (sprintf("%1.0f", Freq)), sep ='\n')), vjust = 0.4, size=10,
            color = t$text.color ) +
  scale_fill_gradient( low = "#F5C710",
                       high = "#DF536B",
                       space = "Lab",
                       na.value = "white",
                       trans="log",
                       guide = "none") +
  ylab("Multi-omics defined") +
  xlab("")+
  ggtitle("FACS-defined lable: Atypical non-PB") +
  theme_bw()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  coord_fixed(ratio=1)+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 29))

ggsave("confusion_mat_atypnonPB.png", height = 4, width = 15, bg = "white")

### FACS-like defined: Naive ----

t <- confusion_matrix %>%
  filter(Var1 == "Naive")

ggplot(data = t, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = paste(gsub("\\s+", "\n", Var2), (sprintf("%1.0f", Freq)), sep ='\n')), vjust = 0.4, size=10,
            color = t$text.color ) +
  scale_fill_gradient(low = "#F5C710",
                      high = "#DF536B",
                      space = "Lab",
                      na.value = "white",
                      trans="log",
                      guide = "none") +
  ylab("Multi-omics defined") +
  xlab("")+
  ggtitle("FACS-defined lable: Naive B cells") +
  theme_bw()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  coord_fixed(ratio=1)+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 29))

ggsave("confusion_mat_naive.png", height = 4, width = 15, bg = "white")

### FACS-like defined: Switched memory ----

t <- confusion_matrix %>%
  filter(Var1 == "Switched memory")

ggplot(data = t, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = paste(gsub("\\s+", "\n", Var2), (sprintf("%1.0f", Freq)), sep ='\n')), vjust = 0.4, size=10,
            color = t$text.color ) +
  scale_fill_gradient( low = "#F5C710",
                       high = "#DF536B",
                       space = "Lab",
                       na.value = "white",
                       trans="log",
                       guide = "none") +
  ylab("Multi-omics defined") +
  xlab("")+
  ggtitle("FACS-defined lable: Switched memory") +
  theme_bw()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  coord_fixed(ratio=1)+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 29))

ggsave("confusion_mat_mem_non+pb.png", height = 4, width = 15, bg = "white")



#  Isotype usages ----
#Function to find the frequency of each Isotype, grouped by broad B cell annotation 
plot_isotypes_distribution <- function(pie_data, plot_title) {
  
  texts <- list(element_text(colour = "white", face = "bold", size = 14.5,  margin = margin(7, 7, 7, 7),
  ))
  background <- list( element_rect(fill ="#01665E",
                                   linetype = 1,
                                   linewidth =unit(1, "npc"),
                                   colour = "white"))
  strip <- strip_themed(background_x =background,
                        by_layer_x= T,
                        clip = "on",
                        text_x = texts)
  
  dim2 <- pie_data %>%
    mutate(csum = rev(cumsum(rev(Percentage))), 
           pos = Percentage/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Percentage/2, pos)) %>%
    filter(Percentage > 0)
  
  iso_plot <- ggplot(pie_data, aes(x = "", y = Percentage, fill = Var2)) +
    geom_col(color = "white", width = 50, just = 1) +
    scale_fill_manual(values = unique(pie_data$ig_color)) + 
    coord_polar(theta = "y") +
    facet_wrap2( ~(paste0("Multi-omics label: ", Var1)), nrow=2, ncol=4,
                 labeller = "label_value",
                 strip = strip, trim_blank = FALSE
    ) +
    theme_void() +
    guides(fill = guide_legend(title = "Isotype")) +
    geom_label_repel(data = dim2,
                     aes(y = pos, label = paste0(Percentage, "%")),
                     size = 5.5, nudge_x = 1, show.legend = FALSE) +
    theme(legend.text = element_text(size = 13))
  
  iso_plot + ggtitle(plot_title)
}

## Isotype usages of FACS-defined subset: Atypical non PB ----

atyp_nonPB_pie_data <- as.data.frame(table(atyp_nonPB$Level_one_annotation, atyp_nonPB$c_gene_HC)) %>%
  group_by(Var1) %>%
  mutate(Percentage = round((Freq/sum(Freq) * 100), 2)) %>%
  mutate(ig_color = case_when(Var2 == "IGHA1" ~ "#DF536B",
                              Var2 == "IGHA2" ~ "#D95F02",
                              Var2 == "IGHD" ~ "#E6AB02",
                              Var2 == "IGHE" ~ "#E7298A",
                              Var2 == "IGHG1" ~ "#66A61E",
                              Var2 == "IGHG2" ~ "#7570B3",
                              Var2 == "IGHG3" ~ "#A6761D",
                              Var2 == "IGHG4" ~ "#2297E6",
                              Var2 == "IGHM" ~ "#1B9E77",
                              Var2 == "None" ~ "#9E9E9E"))

plot_isotypes_distribution(atyp_nonPB_pie_data, "")
ggsave("atyp_iso_dist_F.png", height = 10, width = 20, bg = "white")

## Isotype usages of FACS-defined subset: Switched memory non PB ----
mem_nonpb_pie_data <- as.data.frame(table(mem_nonpb$Level_one_annotation, mem_nonpb$c_gene_HC)) %>%
  group_by(Var1) %>%
  mutate(Percentage = round((Freq/sum(Freq) * 100), 2)) %>%
  mutate(ig_color = case_when(Var2 == "IGHA1" ~ "#DF536B",
                              Var2 == "IGHA2" ~ "#D95F02",
                              Var2 == "IGHD" ~ "#E6AB02",
                              Var2 == "IGHE" ~ "#E7298A",
                              Var2 == "IGHG1" ~ "#66A61E",
                              Var2 == "IGHG2" ~ "#7570B3",
                              Var2 == "IGHG3" ~ "#A6761D",
                              Var2 == "IGHG4" ~ "#2297E6",
                              Var2 == "IGHM" ~ "#1B9E77",
                              Var2 == "None" ~ "#9E9E9E"))

plot_isotypes_distribution(mem_nonpb_pie_data, "")
ggsave("FACS-defined_memory_nonPB/mem_nonpb_iso_dist_F.png", height = 10, width = 20, bg = "white")


## Isotype usages of FACS-defined subset: Naive ----
Naive_pie_data <- as.data.frame(table(Naive_facs$Level_one_annotation, Naive_facs$c_gene_HC)) %>%
  group_by(Var1) %>%
  mutate(Percentage = round((Freq/sum(Freq) * 100), 2)) %>%
  mutate(ig_color = case_when(Var2 == "IGHA1" ~ "#DF536B",
                              Var2 == "IGHA2" ~ "#D95F02",
                              Var2 == "IGHD" ~ "#E6AB02",
                              Var2 == "IGHE" ~ "#E7298A",
                              Var2 == "IGHG1" ~ "#66A61E",
                              Var2 == "IGHG2" ~ "#7570B3",
                              Var2 == "IGHG3" ~ "#A6761D",
                              Var2 == "IGHG4" ~ "#2297E6",
                              Var2 == "IGHM" ~ "#1B9E77",
                              Var2 == "None" ~ "#9E9E9E"))

plot_isotypes_distribution(Naive_pie_data, "")
ggsave("naive_iso_dist_F.png", height = 10, width = 20, bg = "white")


#Function to find the frequency of each Isotype, grouped by FACS-like annotation
plot_isotypes_distribution <- function(pie_data, plot_title) {
  
  texts <- list(element_text(colour = "white", face = "bold", size = 14.5,  margin = margin(7, 7, 7, 7),
  ))
  background <- list( element_rect(fill ="#01665E",
                                   linetype = 1,
                                   linewidth =unit(1, "npc"),
                                   colour = "white"))
  strip <- strip_themed(background_x =background,
                        by_layer_x= T,
                        clip = "on",
                        text_x = texts)
  
  dim2 <- pie_data %>%
    mutate(csum = rev(cumsum(rev(Percentage))), 
           pos = Percentage/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Percentage/2, pos)) %>%
    filter(Percentage > 0)
  
  iso_plot <- ggplot(pie_data, aes(x = "", y = Percentage, fill = Var2)) +
    geom_col(color = "white", width = 50, just = 1) +
    scale_fill_manual(values = unique(pie_data$ig_color)) + 
    coord_polar(theta = "y") +
    facet_wrap2( ~(paste0("FACS label: ", Var1)), nrow=2, ncol=4,
                 labeller = "label_value",
                 strip = strip, trim_blank = FALSE
    ) +
    theme_void() +
    guides(fill = guide_legend(title = "Isotype")) +
    geom_label_repel(data = dim2,
                     aes(y = pos, label = paste0(Percentage, "%")),
                     size = 5.5, nudge_x = 1, show.legend = FALSE) +
    theme(legend.text = element_text(size = 13))
  
  iso_plot + ggtitle(plot_title)
}



## Isotype usages of multi-omics-defined subset: Memory (switched) ----

memory_pie_data <- as.data.frame(table(memory_switched$CITESeq_annot, memory_switched$c_gene_HC)) %>%
  group_by(Var1) %>%
  mutate(Percentage = round((Freq/sum(Freq) * 100), 2)) %>%
  mutate(ig_color = case_when(Var2 == "IGHA1" ~ "#DF536B",
                              Var2 == "IGHA2" ~ "#D95F02",
                              Var2 == "IGHD" ~ "#E6AB02",
                              Var2 == "IGHE" ~ "#E7298A",
                              Var2 == "IGHG1" ~ "#66A61E",
                              Var2 == "IGHG2" ~ "#7570B3",
                              Var2 == "IGHG3" ~ "#A6761D",
                              Var2 == "IGHG4" ~ "#2297E6",
                              Var2 == "IGHM" ~ "#1B9E77",
                              Var2 == "None" ~ "#9E9E9E"))

plot_isotypes_distribution(memory_pie_data, "")
ggsave("mem_iso_dist_F.png", height = 10, width = 20, bg = "white")

# CD27 expression ----

plot_genes_expression <- function(group, group_by) {
  
  VlnPlot(group, features = "CD27", ncol = 1 ,group.by = 
            group_by) +
    ggtitle("")+
    scale_fill_brewer(palette = "Dark2") +
    ylab("CD27") +
    xlab("") +
    ylim(0,3)+
    theme(legend.position = "none",
          plot.tag = element_text(size = 15)) +
    stat_compare_means(method = "anova", vjust = 0)
}

## FACS-defined subset: Atypical non-PB ----
p <- plot_genes_expression(atyp_nonPB, "Level_one_annotation")
p <-  p + ggtitle("FACS-defined: Atypical B cells")
p + xlab("multi-omics defined groups")
ggsave(filename = paste0("CD27_gene_expression_atypicalnonPB",".png"), 
       height = 5, width = 5, bg = "white")


## FACS-defined subset: Naive ----
p <- plot_genes_expression(Naive_facs, "Level_one_annotation")
p <-  p + ggtitle("FACS-defined: Naive")
p + xlab("multi-omics defined groups")
ggsave("CD27_genes_expression_Naive.png", height = 5, width = 5, bg = "white")


## FACS-like defined subset: Switched memory ----

p <- plot_genes_expression(mem_nonpb, "Level_one_annotation")
p <-  p + ggtitle("FACS-defined: Switched memory")
p + xlab("multi-omics defined groups")
ggsave("CD27_genes_expression_facs_Switched_memory.png", height = 5, width = 5, bg = "white")

## multi-omics defined subset: Memory (switched) ----

memory_switched2 <- subset(memory_switched, subset = CITESeq_annot != "NA")
p <- plot_genes_expression(memory_switched2, "CITESeq_annot")
p <-  p + ggtitle("Multi-omics defined: Switched memory")
p + xlab("FACS-defined groups")
ggsave("CD27_genes_expression_omics_Switched_memory.png", height = 5, width = 5, bg = "white")


# Correlation between protein and RNA expression ----
calculate_corr_and_plot <- function(seurat_object, group_annotation, gene, protein,  legend_title) {
  
  corr <- FeatureScatter(seurat_object, feature1 = protein, feature2 = gene, group.by = group_annotation)
  
  correlation_test <- cor.test(as.vector(corr$data[[1]]), corr$data[[2]])
  
  annotations <- data.frame(
    xpos = c(-Inf, -Inf, Inf, Inf),
    ypos = c(-Inf, Inf, -Inf, Inf),
    annotateText = c("", paste0("cor_coefficient= ", round(correlation_test$estimate, digits = 2)), "", paste0("p-value= ", correlation_test$p.value)),
    hjustvar = c(0, -0.2, 1, 1),
    vjustvar = c(0, 1, 0, 1)
  )
  
  result_plot <- corr + geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
                                  size=5) +
    
    ggtitle("") +
    ylab(paste0(gene, " gene expression"))+
    xlab(paste0(sub("ADT-", "", protein), " protein expression"))+
    guides(color = guide_legend(title = legend_title))
  
  return(result_plot)
}

plot_expression_corr <- function(group, group_by, title, legend_title) {
  
  p1 <- calculate_corr_and_plot(group, group_by, "CD27", "ADT-CD27", legend_title)
  p2 <- calculate_corr_and_plot(group, group_by, "CD22", "ADT-CD22", legend_title)
  p3 <- calculate_corr_and_plot(group,group_by, "CR2", "ADT-CD21", legend_title)
  p4 <- calculate_corr_and_plot(group, group_by, "ITGAX", "ADT-CD11c", legend_title)
  
  (p1 + p2 + p3 + p4) + 
    plot_annotation(title = title, tag_levels = 'A') + 
    plot_layout(guides = "collect", ncol = 2, nrow=2) & 
    theme(plot.tag = element_text(size = 15), 
          plot.title = element_text(size =25, hjust = 0.5, family = "Times",
                                    face="bold"))
}

## FACS-defined subset: Atypical non-PB ----
plot_expression_corr(atyp_nonPB, "Level_one_annotation","FACS-defined: Atypical non-PB cells","Multi-omics defined groups")
ggsave("corr_atyp_F.png", height = 10, width = 15, bg = "white")

## FACS-defined subset: Naive ----
plot_expression_corr(Naive_facs, "Level_one_annotation","FACS-defined: Naive B cells","Multi-omics defined groups")
ggsave("corr_naive_F.png", height = 10, width = 15, bg = "white")

## FACS-defined subset: Switched memory non PB ----
plot_expression_corr(mem_nonpb, "Level_one_annotation","FACS-defined: Switched memory","Multi-omics defined groups")
ggsave("FACS-defined_memory_nonPB/corr_memory_nonpb_F.png", height = 10, width = 15, bg = "white")

## multi-omics defined subset: Memory (switched) ----
plot_expression_corr(memory_switched, "CITESeq_annot","Multi-omics defined: Memory switched B cells","FACS-defined groups")
ggsave("Final3/corr_mem_sw_F.png", height = 10, width = 15, bg = "white")

# Healthy only analysis ----
healthy_only <- subset(ess, Source_abrev == "HV")

## Differences in lambda and kappa ratios ----

healthy_only_ka_lam <- data.frame(Level_one_annotation = healthy_only$Level_one_annotation,
                                  CITESeq_annot = healthy_only$CITESeq_annot,
                                  chain_LC = healthy_only$chain_LC,
                                  scRNASeq_sample_ID = healthy_only$scRNASeq_sample_ID
)

Group1loop <- unique(na.omit(healthy_only_ka_lam$CITESeq_annot))

for (celltype in Group1loop) {
  
  t <- healthy_only_ka_lam %>%
    filter(CITESeq_annot == celltype)%>%
    group_by(Level_one_annotation, scRNASeq_sample_ID) %>%
    drop_na() %>%
    count(chain_LC) %>%
    mutate(Sum_Value = sum(n)) %>%
    filter(Sum_Value > 1) %>%
    filter(chain_LC == "IGK") %>%
    mutate(kappa_percentage = n/Sum_Value *100) %>%
    select(Level_one_annotation,kappa_percentage, scRNASeq_sample_ID) %>%
    group_by(Level_one_annotation) %>%
    filter(n() > 2) %>%
    ungroup()
  
  if (nrow(t) == 0) { next }
  
  try(
    stat.test <- t %>%
      t_test(kappa_percentage ~ Level_one_annotation, p.adjust.method = "bonferroni") %>% 
      add_xy_position(x = "Level_one_annotation")
  )
  try(
    pp <- ggboxplot(
      t, x = "Level_one_annotation", y = "kappa_percentage", 
      color = "Level_one_annotation", palette = c("#E69F00", "#0072B2", "#35978F", "#D55E00","#CC79A7", "#1B9E77"),
      bxp.errorbar = TRUE, add = "dotplot", 
      title = paste0("FACS-defined: ", celltype))
  )
  
  pp <- pp + stat_pvalue_manual(
    stat.test,  label = "p", tip.length = 0,  hide.ns = TRUE,
    step.increase = 0.05, color = "Level_one_annotation",
    step.group.by = "Level_one_annotation",
  ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    xlab("")+
    ylim(0,110)+
    ylab("kappa/lambda %")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17 ),
          axis.text.y = element_text(size = 17),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 27, hjust = 0.5),
          legend.position = "none") +
    stat_n_text(y.pos = 110, size = 6)
  
  #save
  ggsave(filename = paste0("k_l_percentage/CITElambda_kappa_ratio_",celltype,".png"), 
         height = 10, width = 15, bg = "white")
  
}

## IGHV4_34 percentage ----

healthy_only_IGHV4_34 <- data.frame(Level_one_annotation = healthy_only$Level_one_annotation,
                                    CITESeq_annot = healthy_only$CITESeq_annot,
                                    scRNASeq_sample_ID = healthy_only$scRNASeq_sample_ID,
                                    v_call_HC = healthy_only@meta.data$v_call_HC)



Group1loop <- unique(healthy_only_IGHV4_34$CITESeq_annot)

for (celltype in Group1loop) {
  try(
    t <- healthy_only_IGHV4_34 %>%
      #loop through omics-defined populations
      filter(CITESeq_annot == celltype )%>%
      group_by(Level_one_annotation, scRNASeq_sample_ID) %>%
      drop_na() %>%
      #filter any sample/FACS-group has one cell only.
      filter(n() > 1) %>%
      ungroup() %>%
      group_by(Level_one_annotation, scRNASeq_sample_ID) %>%
      #find the ratio of IGHV4-34/all per sample per FACS-group
      summarise(IGHV4_34_ratio = sum(str_detect(v_call_HC, "IGHV4-34"))) %>%
      mutate(percentage = IGHV4_34_ratio/sum(IGHV4_34_ratio) *100) %>%
      replace(is.na(.), 0) %>%
      ungroup() %>%
      group_by(Level_one_annotation) %>%
      filter(n() > 2) %>%
      ungroup()
  )
  
  if (nrow(t) == 0) { next }
  
  try(
    #plot
    stat.test <- t %>%
      t_test(percentage ~ Level_one_annotation, p.adjust.method = "bonferroni") %>% 
      add_xy_position(x = "Level_one_annotation")
    
  )
  pp <- ggboxplot(
    t, x = "Level_one_annotation", y = "percentage", 
    color = "Level_one_annotation", palette = c("#E69F00", "#0072B2", "#35978F", "#D55E00","#CC79A7", "#1B9E77"),
    bxp.errorbar = TRUE, add = "dotplot", 
    title = paste0("FACS-defined: ", celltype))
  
  try(
    pp <- pp + stat_pvalue_manual(
      stat.test,  label = "p", tip.length = 0.1,  hide.ns = TRUE,
      step.increase = 0.05,
      y.position = 50, label.size = 15) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      xlab("")+
      ylim(0,100)+
      ylab("IGHV4-34 %")+
      theme(axis.text.x = element_text(size = 25, vjust = 0.5),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 27),
            plot.title = element_text(size = 30, hjust = 0.5),
            legend.position = "none") +
      stat_n_text(y.pos = 90, size = 15)
    
  )
  ggsave(filename = paste0("IGHV4_34_percent/facs_label_",celltype,".png"), 
         height = 7, width =10, bg = "white")
  
}



## Isotypes usages ----


healthy_only_isotypes <- data.frame(Level_one_annotation = healthy_only$Level_one_annotation,
                                    CITESeq_annot = healthy_only$CITESeq_annot,
                                    c_gene_HC =healthy_only$c_gene_HC,
                                    scRNASeq_sample_ID = healthy_only$scRNASeq_sample_ID
)

for (celltype in unique(na.omit(healthy_only_isotypes$CITESeq_annot))) {
  for (IG in unique(healthy_only_isotypes$c_gene_HC)) {
    
    t <- healthy_only_isotypes %>%
      filter(CITESeq_annot == celltype) %>% #loop through multi-omics defined groups, one at a time.
      select(c_gene_HC, scRNASeq_sample_ID,Level_one_annotation) %>%
      #remove any "None" values in isotypes column.
      filter(c_gene_HC != "None") %>%
      #group by sample & FACS-defined groups, so the downstream counting will be conducted on the grouped populations
      group_by(scRNASeq_sample_ID, Level_one_annotation) %>% 
      #count how many isotype we have per patient per cell type (FACS-defined)
      count(c_gene_HC) %>% 
      #sum the total IG per sample per cell type.
      mutate(Sum_Value = sum(n)) %>% 
      #filter out any samples <= 5 no. of cells
      filter(Sum_Value > 5) %>% 
      #divied the count by sum, to normalize the data.
      mutate(normalized_values = n/Sum_Value *100) %>% 
      select(scRNASeq_sample_ID, c_gene_HC, normalized_values,Level_one_annotation) %>%
      group_by(scRNASeq_sample_ID, Level_one_annotation) %>%
      #loop through all isotypes
      filter(c_gene_HC == IG) %>%
      ungroup() %>%
      group_by(Level_one_annotation) %>%
      filter(n() > 2) %>%
      ungroup()
    
    if (nrow(t) == 0) { next }
    
    try(
      stat.test <- t %>%
        t_test(normalized_values ~ Level_one_annotation, p.adjust.method = "bonferroni") %>% 
        add_xy_position(x = "Level_one_annotation")
    )
    
    try(
      
      pp <- ggboxplot(
        t, x = "Level_one_annotation", y = "normalized_values", 
        color = "Level_one_annotation", palette = c("#E69F00", "#0072B2", "#35978F", "#D55E00","#CC79A7", "#1B9E77"),
        bxp.errorbar = TRUE, add = "dotplot", 
        title = paste0("FACS-defined: ", celltype))
    )
    try(
      pp + stat_pvalue_manual(
        stat.test,  label = "p", tip.length = 0.03,  hide.ns = TRUE,
        step.increase = 0.05,
        y.position = 120, label.size = 10) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
        xlab("")+
        ylim(0,130)+
        ylab(paste0(IG," %"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17 ),
              axis.text.y = element_text(size = 17),
              axis.title.y = element_text(size = 18),
              plot.title = element_text(size = 27, hjust = 0.5),
              legend.position = "none") +
        stat_n_text(y.pos = 110, size = 7)
    )
    
    ggsave(filename = paste0("isotypes_usages/facs_label_",celltype,"_", IG,".png"), 
           height = 10, width = 15, bg = "white")
    
    stat.test <- NULL
    pp <- NULL
    
  }
}


# Between diseases analysis Using ggstatsplot ----

cell_comps_percentage <- data.frame(Level_one_annotation = ess$Level_one_annotation,
                                    CITESeq_annot = ess$CITESeq_annot,
                                    scRNASeq_sample_ID = ess$scRNASeq_sample_ID,
                                    patient_group = ess$Source_abrev)

Group3loop <- unique(na.omit(cell_comps_percentage$CITESeq_annot))
Group2loop <- unique(cell_comps_percentage$Level_one_annotation)

for (celltype in Group3loop) {
  for (omicstype in Group2loop) {
    try( 
      t2 <- cell_comps_percentage %>%
        filter(CITESeq_annot == celltype) %>%
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
    
    ggsave(filename = paste0("between_patients_group_analysis3/",omicstype, "_",celltype,".png"), 
           height = 5.5, width = 5.5, bg = "white")
    
  }
}

# Final stats (accuracy, specificity, and sensitivity) ----

confusion_matrix<-as.data.frame(table(ess@meta.data$CITESeq_annot, ess@meta.data$Level_one_annotation))
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

write.csv(Table_of_analysis, paste0("Tables/Table_of_acc_sen_sp.csv"))

# Find how many patients per group ----

disease_groups <- data.frame(scRNASeq_sample_ID = ess$scRNASeq_sample_ID,
                             patient_group = ess$Source_abrev)

disease_groups <- disease_groups[!duplicated(disease_groups$scRNASeq_sample_ID), ]


disease_groups<- disease_groups %>%
  group_by(patient_group) %>%
  count()

write.csv(disease_groups, "Tables/patients_pre_disease_groups.csv")
