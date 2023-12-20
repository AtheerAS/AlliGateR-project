# Gating by literatures markers

# Load required libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(tidyr)
library(proxy)
library(vegan)
library(corrplot)

# Load Seurat object (broad B cell annotation and FACS-like annotation included) ----
ess = readRDS(file ="CITE_and_level1_anno.rds")

# Gate by references markers ----

## Double negative B cells ----
DNB <- subset(data_pbc, subset = `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`<4)

## Anergic B cells 1 ----
AnB1 <- subset(data_pbc, subset = `adt_ADT-CD19`>4 & `adt_ADT-CD27`<5 & `adt_ADT-IgD`> 4 & `adt_ADT-IgM` < 2 )

## Anergic B cells 2 ----
AnB2 <- subset(data_pbc, subset = `adt_ADT-CD19`>4 & `adt_ADT-CD38`< 5 & `adt_ADT-CD21` < 6)

## Naive B cells ----
#Wg gated IgM by 8 range was 0-25, 25/3 (high,mid,low), low = 8.33 ~ 8
NaiveB <- subset(data_pbc, subset =  `adt_ADT-CD27`<5 & `adt_ADT-IgD`> 4 & `adt_ADT-IgM` < 8)

## Atypical B cells ----
AtypB <- subset(data_pbc, subset = `adt_ADT-CD19`>4 &`adt_ADT-CD27`<5 & `adt_ADT-CD20` > 5 & `adt_ADT-CD21` < 6 & `adt_ADT-CD10` < 4)

# Age-associated B cells ----
AgeB <- subset(data_pbc, subset = `adt_ADT-CD19`>4  & `adt_ADT-CD21` < 6 & `adt_ADT-CD11c` > 4 )

# Add multi-omics broad B cell annotation to each datasets ----

## Double negative B cells ----
DNB$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(DNB@meta.data))
for (row_name in common_rows) {
  DNB$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

## Anergic B cells 1 ----
AnB1$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(AnB1@meta.data))
for (row_name in common_rows) {
  AnB1$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

## Anergic B cells 2 ----
AnB2$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(AnB2@meta.data))
for (row_name in common_rows) {
  AnB2$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

## Naive B cells ----
NaiveB$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(NaiveB@meta.data))
for (row_name in common_rows) {
  NaiveB$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

## Atypical B cells ----
AtypB$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(AtypB@meta.data))
for (row_name in common_rows) {
  AtypB$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

# Age-associated B cells ----
AgeB$Multi_omics_label <- NA
common_rows <- intersect(row.names(data_pbc@meta.data), row.names(AgeB@meta.data))
for (row_name in common_rows) {
  AgeB$Multi_omics_label[row_name] <- data_pbc$Level_one_annotation[row_name]
}

list_of_mod_object <- list(AgeB, AtypB, AnB1, AnB2, NaiveB, DNB)

#saveRDS(list_of_mod_object, file = "reference_labels_list.rds")

# Change Seurat object to dataframe ----
#function to convert
process_to_df <- function(Object, reference_id) {
  object_df <- data.frame(multi_omics_label = Object$Multi_omics_label, reference_id = reference_id)
  object_df <- object_df %>% group_by(reference_id, multi_omics_label) %>% count()
  return(object_df)
}


objects_list <- list(AgeB, AtypB, AnB1, AnB2, NaiveB, DNB)

#link objects to loop through them
data_list <- list(
  list(Object = AgeB, reference_id = "Age associated B cells"),
  list(Object = AtypB, reference_id = "CD21 Atypical B cells"),
  list(Object = AnB1, reference_id = "Anergic B cells 1"),
  list(Object = AnB2, reference_id = "Anergic B cells 2"),
  list(Object = NaiveB, reference_id = "IgMlo naïve B cells"),
  list(Object = DNB, reference_id = "Double negative B cells")
)

result_list <- list()

for (pair in data_list) {
  Object <- pair$Object
  reference_id <- pair$reference_id
  result <- process_to_df(Object, reference_id)
  result_list[[length(result_list) + 1]] <- result
}

combined_result <- do.call(rbind, result_list)
#write.csv(combined_result, "cell_density_between_refere_and_omics.csv", quote = F)

# Bubble plot ----

ggplot(combined_result, aes(x = combined_result$reference_id, y = combined_result$multi_omics_label, size = (combined_result$n))) +
  geom_point(aes(color = combined_result$n)) +
  theme_light()+
  labs(x = "Literature label", color = "Cell density", y = "Multi-Omics label")+
  scale_color_gradient(low = "darkseagreen", high = "darkgreen") +
  guides(size = "none") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave("bubble_plot.png", height = 3, width = 4, bg = "white")

#  Correlation matrix ----

#1 Naive
naive <- combined_result %>% filter(reference_id == "IgMlo naïve B cells")
naive <- naive %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)

#DNB
dnbt <- combined_result %>% filter(reference_id == "Double negative B cells")
dnbt <- dnbt %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)
#Age
age <- combined_result %>% filter(reference_id == "Age associated B cells")
age <- age %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)

#atypical
atyp <- combined_result %>% filter(reference_id == "CD21 Atypical B cells")
atyp <- atyp %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)

#Aneg1
aneg1 <- combined_result %>% filter(reference_id == "Anergic B cells 1")
aneg1 <- aneg1 %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)

#Aneg2
aneg2 <- combined_result %>% filter(reference_id == "Anergic B cells 2")
aneg2 <- aneg2 %>% group_by(reference_id) %>% mutate(percentage = (n / sum(n)) * 100)


# jsd_matrix ----

data_frames <- list(aneg2, aneg1, atyp, age, dnbt, naive)

calculate_jsd <- function(df1, df2) {
  # Combine data frames and create a contingency table
  combined_df <- rbind(df1, df2)
  contingency_table <- table(combined_df$multi_omics_label)
  
  # Calculate probability distributions
  prob_dist1 <- df1$n / sum(df1$n)
  prob_dist2 <- df2$n / sum(df2$n)
  
  # Calculate Jensen-Shannon divergence
  m <- 0.5 * (prob_dist1 + prob_dist2)
  jsd <- sum(0.5 * (prob_dist1 * log2(prob_dist1 / m) + prob_dist2 * log2(prob_dist2 / m)))
  
  return(jsd)
}

# Create a vector to store cell group names
cell_group_names <- sapply(data_frames, function(df) unique(df$reference_id))

jsd_matrix <- matrix(0, nrow = length(data_frames), ncol = length(data_frames))
# Populate the matrix with Jensen-Shannon divergence values
for (i in 1:length(data_frames)) {
  for (j in 1:length(data_frames)) {
    jsd_matrix[i, j] <- calculate_jsd(data_frames[[i]], data_frames[[j]])
  }
}

# Print the matrix with cell group names
print(jsd_matrix)
colnames(jsd_matrix) <- cell_group_names
rownames(jsd_matrix) <- cell_group_names
print(jsd_matrix)
#Convert to df
jsd_df <- as.data.frame(as.table(jsd_matrix))
colnames(jsd_df) <- c("Cell_Group_1", "Cell_Group_2", "JSD")

jsd_df$invers <- 1- jsd_df$JSD
jsd_df$percentage <- jsd_df$invers * 100

#Convert to matrox again
jsd_per_mat <- jsd_df %>%
  select(Cell_Group_1, Cell_Group_2, invers) %>%
  spread(key = Cell_Group_1, value = invers) %>%
  tibble::column_to_rownames(var = "Cell_Group_2") %>%
  as.matrix()

#Plot
tiff("jsd_per_mat.tiff", units="in", width=5, height=5, res=300)

corrplot(jsd_per_mat, type="upper", order = 'alphabet', addCoef.col = 'white',
         method="color",
         tl.cex=1,
         outline=TRUE,
         tl.col="black",
         col.lim = c(0, 1),
         col = COL1('YlGn'),
         addgrid.col = 'white')

dev.off()






