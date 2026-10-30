# Air and environmental sampling expand metagenomic viral surveillance
# in Cambodian pig slaughterhouses
#
# Analysis script for capture-enriched metagenomic surveillance.
# Run the setup below, then the numbered analysis sections you need.
# The original analysis calculations, thresholds and export choices are retained.
#
# Input layout, relative to the R working directory:
#   Pig_*.xlsx                        Count, detection, metadata and taxonomy files
#   Working_Data/Pig_Counts/           Raw CZ ID taxon-count CSV exports
#   Working_Data/Pig_Contig_Taxonomy/  Raw contig-taxonomy CSV exports
#   DIAMOND/                          Annotated DIAMOND workbooks
#   Trees/                            Torovirus.Spike.tre and HKU15.Spike.tre
#
# Preprocessing sections create objects in memory. Their checkpoint write.xlsx()
# calls remain commented, as in the source script. To run preprocessing in order,
# enable those checkpoint exports before the next section reads the workbook, or
# provide the existing intermediate workbooks. Later analyses also require the
# separately annotated DIAMOND workbooks and precomputed phylogenetic trees.
# Commented figure/table exports are optional; active exports go to output_dir.
# Data and tree files are not included with this script.

# Setup: portable file locations ----
# Start R in the project/data folder. No personal or absolute paths are required.
data_dir <- "."
output_dir <- "results"
counts_dir <- file.path(data_dir, "Working_Data", "Pig_Counts")
contig_taxonomy_dir <- file.path(data_dir, "Working_Data", "Pig_Contig_Taxonomy")
tree_dir <- file.path(data_dir, "Trees")

input_file <- function(...) {
  file.path(data_dir, ...)
}

output_file <- function(...) {
  destination <- file.path(output_dir, ...)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  destination
}

# Setup: package dependencies ----
# Install missing packages separately; this script does not install software.
# phyloseq, ComplexHeatmap and ggtree are Bioconductor packages.
# plyr and ggpattern are used through explicit namespaces below.
required_packages <- c("ade4", "vegan", "ape", "phyloseq", "ComplexHeatmap", "circlize", "PMCMRplus",
  "ggrepel", "ggdist", "ggalluvial", "fmsb", "openxlsx", "patchwork", "ggtree", "ggExtra",
  "scales", "grid", "tidyverse", "plyr", "ggnewscale", "ggpattern")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1),
  quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install the required R packages before running: ", paste(missing_packages, collapse = ", "))
}

library(ade4)
library(vegan)
library(ape)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(PMCMRplus)
library(ggrepel)
library(ggdist)
library(ggalluvial)
library(fmsb)
library(openxlsx)
library(patchwork)
library(ggtree)
library(ggExtra)
library(scales)
library(grid)
# Attach tidyverse last so its data-manipulation functions resolve consistently.
library(tidyverse)


# 01. Sample-type colour palette ----
group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")


# 02. Import and summarise raw CZ ID taxon counts ----
myFiles <- list.files(path = counts_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(myFiles) == 0L) {
  stop("No CSV files found in the configured raw-input directory. Check the setup section.")
}

Species <- myFiles %>%
  set_names() %>%
  map_dfr(.f = read_delim, delim = ";", .id = "file_name")

Species_Count <- str_split_fixed(Species$`tax_id,tax_level,genus_tax_id,name,common_name,category,is_phage,agg_score,max_z_score,nt_z_score,nt_rpm,nt_count,nt_contigs,nt_contig_r,nt_percent_identity,nt_alignment_length,nt_e_value,nt_bg_mean,nt_bg_stdev,nt_bg_mean_mass_normalized,nt_bg_stdev_mass_normalized,nr_z_score,nr_rpm,nr_count,nr_contigs,nr_contig_r,nr_percent_identity,nr_alignment_length,nr_e_value,nr_bg_mean,nr_bg_stdev,nr_bg_mean_mass_normalized,nr_bg_stdev_mass_normalized,species_tax_ids,known_pathogen`,
  ",", 35)
Species_Count <- as.data.frame(Species_Count)
# Keep the original filename (including its extension) as the sample identifier.
# basename() works regardless of the location or nesting of the project folder.
Samples <- data.frame(Sample_Name = basename(Species$file_name))
Species_Count <- cbind(Samples, Species_Count)

Species_Type <- Species_Count[c(1, 3, 5, 7, 8, 12, 13, 14, 15, 24, 25, 26, 27)]
Species_List <- subset(Species_Type, V2 == "1")
Virus_List <- subset(Species_List, V6 == "viruses")
Clean_Virus_List <- subset(Virus_List, V7 == "false")
NT_Virus <- Clean_Virus_List[c(1, 3, 6, 7, 8, 9)]
names(NT_Virus)[names(NT_Virus) == "V4"] <- "Species"
names(NT_Virus)[names(NT_Virus) == "V11"] <- "RPM"
names(NT_Virus)[names(NT_Virus) == "V12"] <- "Read_Count"
names(NT_Virus)[names(NT_Virus) == "V13"] <- "Contig_Number"
names(NT_Virus)[names(NT_Virus) == "V14"] <- "Contig_Read_Count"
Final_NT_Virus <- NT_Virus[!(is.na(NT_Virus$Contig_Number) | NT_Virus$Contig_Number == ""),
  ]
NR_Virus <- Clean_Virus_List[c(1, 3, 10, 11, 12, 13)]
names(NR_Virus)[names(NR_Virus) == "V4"] <- "Species"
names(NR_Virus)[names(NR_Virus) == "V23"] <- "RPM"
names(NR_Virus)[names(NR_Virus) == "V24"] <- "Read_Count"
names(NR_Virus)[names(NR_Virus) == "V25"] <- "Contig_Number"
names(NR_Virus)[names(NR_Virus) == "V26"] <- "Contig_Read_Count"
Final_NR_Virus <- NR_Virus[!(is.na(NR_Virus$Contig_Number) | NR_Virus$Contig_Number == ""),
  ]
Combined_Virus <- rbind(Final_NT_Virus, Final_NR_Virus)
Combined_Virus$Names <- paste(Combined_Virus$Sample_Name, "_", Combined_Virus$Species)
Combined_Virus$Contig_Number <- as.numeric(Combined_Virus$Contig_Number)
Ordered_Virus <- Combined_Virus[order(Combined_Virus$Names, -abs(Combined_Virus$Contig_Number)),
  ]
Clean_Combined_Virus <- Ordered_Virus[!duplicated(Ordered_Virus$Names), ]
Virus_Count <- Clean_Combined_Virus[c(1, 2, 6)]
Virus_Count$Contig_Read_Count <- as.numeric(Virus_Count$Contig_Read_Count)
Virus_Count <- aggregate(Contig_Read_Count ~ Species + Sample_Name, Virus_Count, sum)
Virus_Count <- stats::reshape(Virus_Count, idvar = "Sample_Name", timevar = "Species", direction = "wide")
Virus_Count[is.na(Virus_Count)] = 0
# write.xlsx(Virus_Count, file=input_file('Pig_Virus_Count_Raw.xlsx'), rowNames = TRUE)


# 03. Extract GenBank identifiers from contig taxonomy ----
myFiles <- list.files(path = contig_taxonomy_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(myFiles) == 0L) {
  stop("No CSV files found in the configured raw-input directory. Check the setup section.")
}

Species <- myFiles %>%
  set_names() %>%
  map_dfr(.f = read_delim, delim = ";", .id = "file_name")

Species_Count <- str_split_fixed(Species$`contig_name,read_count,contig_length,contig_coverage,NT.species_taxid,NT.genus_taxid,NT.family_taxid,NT.order_taxid,NT.class_taxid,NT.phylum_taxid,NT.kingdom_taxid,NT.superkingdom_taxid,NT.superkingdom_name,NT.Accession,NT.Percentage Identity,NT.Alignment Length,NT.Number of mismatches,NT.E-value,NT.Bitscore,NR.species_taxid,NR.genus_taxid,NR.family_taxid,NR.order_taxid,NR.class_taxid,NR.phylum_taxid,NR.kingdom_taxid,NR.superkingdom_taxid,NR.superkingdom_name,NR.Accession,NR.Percentage Identity,NR.Alignment Length,NR.Number of mismatches,NR.E-value,NR.Bitscore`,
  ",", 34)
Species_Count <- as.data.frame(Species_Count)
Species <- subset(Species_Count, V28 == "Viruses")
GenBank_ID <- Species[c(29)]
GenBank_ID <- GenBank_ID[!duplicated(GenBank_ID$V29), ]
GenBank_ID <- as.data.frame(GenBank_ID)
# write.xlsx(GenBank_ID, file = output_file('DIAMOND_GenBank_ID.xlsx'), rowNames = TRUE)


# 04. Filter low-abundance taxa ----
Species_Raw <- read.xlsx(input_file("Pig_Virus_Count_Raw.xlsx"), rowNames = TRUE)

Species_Norm <- t(apply(Species_Raw,1,function(x)(x/sum(x))))
Species_Norm <- as.data.frame(Species_Norm)
Species_Norm <- Species_Norm*100
Species_Norm_Filt <- Species_Norm[, -which(plyr::numcolwise(sum)(Species_Norm) < 0.1)]
Viruses <- colnames(Species_Norm_Filt)
Species_Count <- Species_Raw[colnames(Species_Raw) %in% Viruses]
# write.xlsx(Species_Count, file=input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)


# 05. Select taxa detected in pig samples ----
Species_Count <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

pig_ids <- rownames(Metadata[Metadata$Sample_Type == "Pig", , drop = FALSE])
pig_only <- Species_Count[rownames(Species_Count) %in% pig_ids, , drop = FALSE]
if (nrow(pig_only) == 0) {
  keep_cols <- character(0)
} else {
  keep_cols <- names(which(colSums(as.matrix(pig_only) > 0, na.rm = TRUE) > 0))
}

Swine_Virus_Count <- Species_Count[, keep_cols, drop = FALSE]
# write.xlsx(Swine_Virus_Count, file = input_file('Pig_Only_RAW_Virus_Count.xlsx'), rowNames = TRUE)


# 06. Apply curated taxon exclusions and merge astrovirus labels ----
Species_Count <- read.xlsx(input_file("Pig_Only_RAW_Virus_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)
Taxonomy <- read.xlsx(input_file("Pig_Taxonomy.xlsx"), rowNames = TRUE)

common_cols <- intersect(colnames(Species_Count), rownames(Taxonomy))
keep_cols <- common_cols[Taxonomy[common_cols, "Keep"] != "No"]
Species_Filt <- Species_Count[, keep_cols, drop = FALSE]

Species_Filt$Astrovirus_swine_PoAstV12_Canada_2006 <- rowSums(Species_Filt[, c("Astrovirus_swine_PoAstV12_3_Canada_2006",
  "Astrovirus_swine_PoAstV14_4_Canada_2006")], na.rm = TRUE)

Species_Filt <- Species_Filt[, !colnames(Species_Filt) %in% c("Astrovirus_swine_PoAstV12_3_Canada_2006",
  "Astrovirus_swine_PoAstV14_4_Canada_2006")]
# write.xlsx(Species_Filt, file = input_file('Pig_Only_Species_Count.xlsx'), rowNames = TRUE)


# 07. Select environmental taxa and retain pig-associated taxa ----
Species_Count <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)
Pig_Species_Count <- read.xlsx(input_file("Pig_Only_Species_Count.xlsx"), rowNames = TRUE)

es_ids <- rownames(Metadata[Metadata$Sample_Group == "Environment", , drop = FALSE])
es_only <- Species_Count[rownames(Species_Count) %in% es_ids, , drop = FALSE]
pig_cols_all <- intersect(colnames(Pig_Species_Count), colnames(Species_Count))
non_pig_cols <- setdiff(colnames(Species_Count), pig_cols_all)
non_pig_cols_in_es <- intersect(non_pig_cols, colnames(es_only))

if (length(non_pig_cols_in_es) == 0) {
  keep_by_thresh <- character(0)
} else {
  totals_es <- colSums(as.matrix(es_only[, non_pig_cols_in_es, drop = FALSE]), na.rm = TRUE)
  grand_es <- sum(totals_es, na.rm = TRUE)
  if (grand_es == 0) {
    keep_by_thresh <- character(0)
  } else {
    keep_by_thresh <- names(totals_es[totals_es/grand_es >= 0.001])
  }
}

keep_cols <- union(pig_cols_all, keep_by_thresh)
ES_Species_Count <- Species_Count[, keep_cols, drop = FALSE]
# write.xlsx(ES_Species_Count, file=input_file('Pig_ES_Only_Count.xlsx'), rowNames = TRUE)


# 08. Calculate reads per million (RPM) ----

Species_Count <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Species_Count$Total_Reads <- rowSums(Species_Count)
Species_RPM <- Species_Count[, -ncol(Species_Count)]/Species_Count$Total_Reads * 1e+06
# write.xlsx(Species_RPM, file=input_file('Pig_Species_RPM.xlsx'), rowNames = TRUE)

ES_Species_Count <- read.xlsx(input_file("Pig_ES_Only_Count.xlsx"), rowNames = TRUE)
ES_Species_Count$Total_Reads <- rowSums(ES_Species_Count)
ES_Species_RPM <- ES_Species_Count[, -ncol(ES_Species_Count)]/ES_Species_Count$Total_Reads *
  1e+06
# write.xlsx(ES_Species_RPM, file=input_file('Pig_ES_Species_RPM.xlsx'), rowNames = TRUE)

Pig_Species_Count <- read.xlsx(input_file("Pig_Only_Species_Count.xlsx"), rowNames = TRUE)
Pig_Species_Count$Total_Reads <- rowSums(Pig_Species_Count)
Pig_Species_RPM <- Pig_Species_Count[, -ncol(Pig_Species_Count)]/Pig_Species_Count$Total_Reads *
  1e+06
# write.xlsx(Pig_Species_RPM, file=input_file('Pig_Only_Species_RPM.xlsx'), rowNames = TRUE)


# 09. Create the pig-associated presence/absence matrix ----
Pig_Species_Count <- read.xlsx(input_file("Pig_Only_Species_Count.xlsx"), rowNames = TRUE)
Pig_Species_Detection <- Pig_Species_Count %>%
  mutate_if(is.numeric, ~1 * (. != 0))
# write.xlsx(Pig_Species_Detection, file=input_file('Pig_Species_Detection.xlsx'), rowNames = TRUE)


# 10. Summarise virus detections with bubble plots ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

Binary_Species <- (Species > 0) * 1

virus_dist <- dist(t(Binary_Species), method = "binary")
virus_hclust <- hclust(virus_dist, method = "average")
clustered_viruses <- virus_hclust$labels[virus_hclust$order]

all_experiments <- sort(unique(Metadata$Experiment))
all_viruses <- colnames(Binary_Species)

count_other_types <- function(rows_other, sample_types_other) {
  if (length(rows_other) == 0)
    return(rep(0L, ncol(Binary_Species)))
  idx_list <- split(rows_other, sample_types_other)
  by_type_present <- lapply(idx_list, function(rr) {
    colSums(Binary_Species[rr, , drop = FALSE]) > 0
  })
  by_type_mat <- do.call(rbind, by_type_present)
  colSums(by_type_mat)
}

bubble_list <- vector("list", length(all_experiments))

for (i in seq_along(all_experiments)) {
  exp_id <- all_experiments[i]
  exp_idx <- Metadata$Experiment == exp_id

  pig_idx <- which(exp_idx & Metadata$Sample_Type == "Pig")
  oth_idx <- which(exp_idx & Metadata$Sample_Type != "Pig")

  pig_total <- length(pig_idx)
  other_total <- length(oth_idx)

  pig_counts <- if (pig_total > 0)
    colSums(Binary_Species[pig_idx, , drop = FALSE]) else rep(0L, length(all_viruses))
  other_counts <- if (other_total > 0)
    colSums(Binary_Species[oth_idx, , drop = FALSE]) else rep(0L, length(all_viruses))

  pig_detected <- pig_counts > 0
  other_detected <- other_counts > 0

  status_code <- integer(length(all_viruses))
  status_code[pig_detected & other_detected] <- 1
  status_code[pig_detected & !other_detected] <- 2
  status_code[!pig_detected & other_detected] <- 3

  other_types_count <- count_other_types(rows_other = oth_idx, sample_types_other = Metadata$Sample_Type[oth_idx])

  bubble_list[[i]] <- tibble(Experiment = exp_id, Experiment_label = paste0("Experiment_",
    exp_id), Virus = all_viruses, status = status_code, n_other_types = as.integer(other_types_count))
}

bubble_df <- bind_rows(bubble_list) %>%
  filter(status > 0) %>%
  mutate(Experiment_label = factor(Experiment_label, levels = paste0("Experiment_", all_experiments)),
    Virus = factor(Virus, levels = clustered_viruses), Outcome = factor(dplyr::recode(status,
      `1` = "Co_detection", `2` = "Pig_only", `3` = "Other_only"), levels = c("Co_detection",
      "Pig_only", "Other_only")))

status_cols_named <- c(Co_detection = "#F8766D", Pig_only = "#F6E09B", Other_only = "#00BFC4")

MAX_BUBBLE_SIZE <- 7
PIG_FIXED_SIZE <- MAX_BUBBLE_SIZE

df_co <- filter(bubble_df, Outcome == "Co_detection")
df_oth <- filter(bubble_df, Outcome == "Other_only")
df_pig <- mutate(filter(bubble_df, Outcome == "Pig_only"), size_fixed = PIG_FIXED_SIZE)

Bubble_Plot <- ggplot() + scale_fill_manual(values = status_cols_named, name = "Outcome") +

geom_point(data = df_co, aes(x = Experiment_label, y = Virus, fill = Outcome, size = n_other_types),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9) + scale_size_area(max_size = MAX_BUBBLE_SIZE,
  breaks = sort(unique(df_co$n_other_types)), name = "Distinct other types (Co-det.)") +

ggnewscale::new_scale("size") + 
geom_point(data = df_oth, aes(x = Experiment_label, y = Virus, fill = Outcome, size = n_other_types),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9) + scale_size_area(max_size = MAX_BUBBLE_SIZE,
  breaks = sort(unique(df_oth$n_other_types)), name = "Distinct other types (Other-only)") +

ggnewscale::new_scale("size") + 
geom_point(data = df_pig, aes(x = Experiment_label, y = Virus, fill = Outcome, size = size_fixed),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9, show.legend = TRUE) + scale_size_identity(name = "Pig-only size",
  breaks = PIG_FIXED_SIZE, labels = "Matches largest", guide = "legend") + 
scale_x_discrete(limits = levels(bubble_df$Experiment_label)) + scale_y_discrete(limits = clustered_viruses) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL, title = "") + theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 9, angle = 45,
    hjust = 1), axis.text.y = element_text(color = "black", face = "bold", size = 9), axis.title.x = element_text(color = "black",
    face = "bold", size = 10), axis.title.y = element_text(color = "black", face = "bold",
    size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "right", legend.title = element_text(face = "bold"))

Bubble_Plot
# ggsave(output_file('Pig_Summary_BubblePlot_A4.pdf'), Bubble_Plot, width = 8, height = 13, units = 'in')

# Alternative horizontal bubble-plot layout
exp_dist <- dist(Binary_Species, method = "binary")
exp_hclust <- hclust(exp_dist, method = "average")

exp_order_samples <- rownames(Binary_Species)[exp_hclust$order]
clustered_experiments <- Metadata[exp_order_samples, , drop = FALSE]$Experiment
clustered_experiments <- unique(clustered_experiments)
clustered_experiments_labels <- paste0("Experiment_", clustered_experiments)

bubble_df <- bind_rows(bubble_list) %>%
  filter(status > 0) %>%
  mutate(Experiment_label = factor(Experiment_label, levels = clustered_experiments_labels),
    Virus = factor(Virus, levels = all_viruses), Outcome = factor(dplyr::recode(status,
      `1` = "Co_detection", `2` = "Pig_only", `3` = "Other_only"), levels = c("Co_detection",
      "Pig_only", "Other_only")))

status_cols_named <- c(Co_detection = "#F8766D", Pig_only = "#F6E09B", Other_only = "#00BFC4")

MAX_BUBBLE_SIZE <- 7
PIG_FIXED_SIZE <- MAX_BUBBLE_SIZE

df_co <- filter(bubble_df, Outcome == "Co_detection")
df_oth <- filter(bubble_df, Outcome == "Other_only")
df_pig <- mutate(filter(bubble_df, Outcome == "Pig_only"), size_fixed = PIG_FIXED_SIZE)

Bubble_Plot <- ggplot() + scale_fill_manual(values = status_cols_named, name = "Outcome") +

geom_point(data = df_co, aes(x = Virus, y = Experiment_label, fill = Outcome, size = n_other_types),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9) + scale_size_area(max_size = MAX_BUBBLE_SIZE,
  breaks = sort(unique(df_co$n_other_types)), name = "Distinct other types (Co-det.)") +

ggnewscale::new_scale("size") + 
geom_point(data = df_oth, aes(x = Virus, y = Experiment_label, fill = Outcome, size = n_other_types),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9) + scale_size_area(max_size = MAX_BUBBLE_SIZE,
  breaks = sort(unique(df_oth$n_other_types)), name = "Distinct other types (Other-only)") +

ggnewscale::new_scale("size") + 
geom_point(data = df_pig, aes(x = Virus, y = Experiment_label, fill = Outcome, size = size_fixed),
  shape = 21, color = "black", stroke = 0.4, alpha = 0.9, show.legend = TRUE) + scale_size_identity(name = "Pig-only size",
  breaks = PIG_FIXED_SIZE, labels = "Matches largest", guide = "legend") + 
scale_x_discrete(limits = all_viruses) + scale_y_discrete(limits = clustered_experiments_labels) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL, title = "") + theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 9, angle = 45,
    hjust = 1, vjust = 1), axis.text.y = element_text(color = "black", face = "bold", size = 9),
    axis.title.x = element_text(color = "black", face = "bold", size = 10), axis.title.y = element_text(color = "black",
      face = "bold", size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "right", legend.title = element_text(face = "bold"))

Bubble_Plot
# ggsave(output_file('Pig_Summary_BubblePlot_A4_swapped.pdf'), Bubble_Plot, width = 9, height = 3, units = 'in')


# 11. Summarise virus detections with heatmaps ----

Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

Binary_Species <- (Species > 0) * 1

virus_dist <- dist(t(Binary_Species), method = "binary")
virus_hclust <- hclust(virus_dist, method = "average")
clustered_viruses <- virus_hclust$labels[virus_hclust$order]

all_experiments <- sort(unique(Metadata$Experiment))
all_viruses <- colnames(Binary_Species)

other_types <- sort(unique(Metadata$Sample_Type[Metadata$Sample_Type != "Pig"]))

all_heat_list <- list()

for (other in other_types) {

  heat_list <- vector("list", length(all_experiments))

  for (i in seq_along(all_experiments)) {
    exp_id <- all_experiments[i]
    exp_idx <- Metadata$Experiment == exp_id

    pig_idx <- which(exp_idx & Metadata$Sample_Type == "Pig")
    oth_idx <- which(exp_idx & Metadata$Sample_Type == other)

    pig_total <- length(pig_idx)
    other_total <- length(oth_idx)

    if (pig_total == 0L || other_total == 0L) {
      heat_list[[i]] <- tibble()
      next
    }

    pig_counts <- colSums(Binary_Species[pig_idx, , drop = FALSE])
    other_counts <- colSums(Binary_Species[oth_idx, , drop = FALSE])

    pig_detected <- pig_counts > 0
    other_detected <- other_counts > 0

    status_code <- integer(length(all_viruses))
    status_code[pig_detected & other_detected] <- 1
    status_code[pig_detected & !other_detected] <- 2
    status_code[!pig_detected & other_detected] <- 3

    heat_list[[i]] <- tibble(Comparison = paste0("Pig_vs_", other), Other_Type = other,
      Experiment = exp_id, Experiment_label = paste0("Experiment_", exp_id), Virus = all_viruses,
      status = status_code)
  }

  all_heat_list[[other]] <- bind_rows(heat_list)
}

heat_df <- bind_rows(all_heat_list) %>%
  filter(status > 0) %>%
  mutate(Experiment_label = factor(Experiment_label, levels = paste0("Experiment_", all_experiments)),
    Virus = factor(Virus, levels = c(clustered_viruses, "Group_annot")), Outcome = factor(dplyr::recode(status,
      `1` = "Co_detection", `2` = "Pig_only", `3` = "Other_only"), levels = c("Co_detection",
      "Pig_only", "Other_only")))

comparison_rank <- heat_df %>%
  group_by(Comparison) %>%
  summarise(n_co = sum(Outcome == "Co_detection"), n_oth = sum(Outcome == "Other_only"),
    score = n_co + n_oth, .groups = "drop") %>%
  arrange(desc(score))

heat_df <- heat_df %>%
  mutate(Comparison = factor(Comparison, levels = comparison_rank$Comparison))

status_cols_named <- c(Co_detection = "#F8766D", Pig_only = "#F6E09B", Other_only = "#00BFC4")

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

ann_df <- heat_df %>%
  distinct(Comparison, Other_Type, Experiment_label) %>%
  mutate(Virus = factor("Group_annot", levels = c(clustered_viruses, "Group_annot")))

Combined_Heatmap <- ggplot() + geom_tile(data = ann_df, aes(x = Experiment_label, y = Virus,
  fill = Other_Type), color = NA) + scale_fill_manual(values = group_palette, name = "Sample type") +

ggnewscale::new_scale_fill() + 
geom_tile(data = heat_df, aes(x = Experiment_label, y = Virus, fill = Outcome), color = "black",
  size = 0.1) + scale_fill_manual(values = status_cols_named, name = "Outcome", drop = FALSE) +

scale_y_discrete(limits = c(clustered_viruses, "Group_annot")) + coord_cartesian(clip = "off") +

facet_grid(~Comparison, scales = "free_x", space = "free_x") + labs(x = NULL, y = NULL, title = "Pig vs other sample types: virus detections by experiment") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.y = element_text(color = "black", face = "bold", size = 8), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), legend.position = "right", legend.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold", size = 10))

Combined_Heatmap
# ggsave(output_file('Pig_Detections_Heatmap_A4.pdf'), Combined_Heatmap, width = 11.5, height = 11.5, units = 'in')


# 12. Compare pig and environmental detections by date ----

species_path <- input_file("Pig_Species_Detection.xlsx")
metadata_path <- input_file("Pig_Metadata.xlsx")

group_palette <- c(Pig = "#EE3377", Environment = "#009988")

species <- read.xlsx(species_path, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID")

metadata <- read.xlsx(metadata_path, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID") %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Date = as.character(Date), Reverse_Date = as.character(Reverse_Date),
    Group = if_else(Sample_Type == "Pig", "Pig", "Environment"))

detections_long <- species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Present") %>%
  filter(Present > 0) %>%
  left_join(metadata %>%
    select(Sample_ID, Sample_Type, Group, Date, Reverse_Date), by = "Sample_ID") %>%
  mutate(Plot_Date = if_else(!is.na(Reverse_Date) & Reverse_Date != "", as.Date(Reverse_Date,
    format = "%Y_%m_%d"), as.Date(Date, format = "%d_%m_%Y"))) %>%
  filter(!is.na(Plot_Date))

# Raw sample-level detections:
# one virus detected in one sample at one time point = one detection.
raw_counts <- detections_long %>%
  count(Plot_Date, Date, Group, name = "Raw_Sample_Detections")

# Unique date-level detections:
# one virus counted once within Pig or Environment on a given date,
# regardless of how many samples in that group were positive on that date.
unique_counts <- detections_long %>%
  distinct(Plot_Date, Date, Group, Virus) %>%
  count(Plot_Date, Date, Group, name = "Unique_Virus_Detections")

summary_long <- full_join(raw_counts, unique_counts, by = c("Plot_Date", "Date", "Group")) %>%
  mutate(Group = factor(Group, levels = c("Pig", "Environment")), Raw_Sample_Detections = coalesce(Raw_Sample_Detections,
    0L), Unique_Virus_Detections = coalesce(Unique_Virus_Detections, 0L)) %>%
  arrange(Plot_Date, Group)

pig_date_unique <- detections_long %>%
  filter(Group == "Pig") %>%
  distinct(Plot_Date, Date, Virus)

environment_date_unique <- detections_long %>%
  filter(Group == "Environment") %>%
  distinct(Plot_Date, Date, Virus)

shared_date_unique <- inner_join(pig_date_unique, environment_date_unique, by = c("Plot_Date",
  "Date", "Virus"))

overlap_summary <- pig_date_unique %>%
  count(Plot_Date, Date, name = "Pig_Unique_Detections") %>%
  full_join(environment_date_unique %>%
    count(Plot_Date, Date, name = "Environment_Unique_Detections"), by = c("Plot_Date",
    "Date")) %>%
  full_join(shared_date_unique %>%
    count(Plot_Date, Date, name = "Shared_Unique_Detections"), by = c("Plot_Date", "Date")) %>%
  mutate(Pig_Unique_Detections = coalesce(Pig_Unique_Detections, 0L), Environment_Unique_Detections = coalesce(Environment_Unique_Detections,
    0L), Shared_Unique_Detections = coalesce(Shared_Unique_Detections, 0L), Pig_Only_Unique_Detections = Pig_Unique_Detections -
    Shared_Unique_Detections, Environment_Only_Unique_Detections = Environment_Unique_Detections -
    Shared_Unique_Detections, Percent_Pig_Recaptured_In_Environment = if_else(Pig_Unique_Detections >
    0, 100 * Shared_Unique_Detections/Pig_Unique_Detections, 0), Percent_Environment_Shared_With_Pig = if_else(Environment_Unique_Detections >
    0, 100 * Shared_Unique_Detections/Environment_Unique_Detections, 0), Date_Label = format(Plot_Date,
    "%d %b %Y")) %>%
  arrange(Plot_Date)

summary_table <- summary_long %>%
  mutate(Date_Label = format(Plot_Date, "%d %b %Y")) %>%
  left_join(overlap_summary %>%
    select(Plot_Date, Percent_Pig_Recaptured_In_Environment, Percent_Environment_Shared_With_Pig),
    by = "Plot_Date") %>%
  mutate(Percent_Found_In_Other_Group = if_else(Group == "Pig", Percent_Pig_Recaptured_In_Environment,
    Percent_Environment_Shared_With_Pig)) %>%
  select(Date_Label, Group, Unique_Virus_Detections, Percent_Found_In_Other_Group) %>%
  arrange(match(Date_Label, unique(format(sort(unique(summary_long$Plot_Date)), "%d %b %Y"))),
    Group)

date_lookup <- summary_long %>%
  mutate(Date_Label = format(Plot_Date, "%d %b %Y")) %>%
  distinct(Plot_Date, Date_Label) %>%
  arrange(Plot_Date) %>%
  mutate(Date_Index = row_number())

summary_wide <- summary_table %>%
  pivot_wider(names_from = Group, values_from = c(Unique_Virus_Detections, Percent_Found_In_Other_Group)) %>%
  left_join(date_lookup, by = "Date_Label") %>%
  arrange(Plot_Date) %>%
  select(-Date_Index, -Plot_Date)

overlap_plot_df <- overlap_summary %>%
  select(Plot_Date, Date_Label, Shared_Unique_Detections, Pig_Only_Unique_Detections) %>%
  pivot_longer(cols = c(Shared_Unique_Detections, Pig_Only_Unique_Detections), names_to = "Overlap_Status",
    values_to = "Count") %>%
  mutate(Overlap_Status = factor(Overlap_Status, levels = c("Shared_Unique_Detections", "Pig_Only_Unique_Detections"),
    labels = c("Recaptured in environment", "Pig only")))

summary_long <- summary_long %>%
  left_join(date_lookup, by = "Plot_Date")

overlap_plot_df <- overlap_plot_df %>%
  left_join(date_lookup, by = "Plot_Date") %>%
  arrange(Date_Index, Overlap_Status)

x_breaks <- date_lookup$Date_Index
x_labels <- date_lookup$Date_Label

counts_plot <- ggplot(summary_long, aes(x = Date_Index, y = Unique_Virus_Detections, color = Group,
  group = Group)) + geom_line(linewidth = 1.1) + geom_point(size = 2.6) + scale_color_manual(values = group_palette) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) + labs(x = NULL, y = "Unique virus detections",
  color = NULL, title = "Unique virus detections in pig and environmental samples over time") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1,
  vjust = 1, color = "black", face = "bold"), axis.text.y = element_text(color = "black",
  face = "bold"), axis.title.y = element_text(color = "black", face = "bold"), panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

percent_plot <- ggplot(overlap_plot_df, aes(x = Date_Index, y = Count, fill = Overlap_Status,
  group = Overlap_Status)) + geom_area(position = "fill", alpha = 0.9, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(`Recaptured in environment` = "#009988", `Pig only` = "#EE3377"),
    breaks = c("Pig only", "Recaptured in environment")) + scale_x_continuous(breaks = x_breaks,
  labels = x_labels) + scale_y_continuous(labels = percent_format()) + labs(x = NULL, y = "Percentage of pig unique detections",
  fill = NULL, title = "Pig unique detections recaptured in the environment over time") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1,
  vjust = 1, color = "black", face = "bold"), axis.text.y = element_text(color = "black",
  face = "bold"), axis.title.y = element_text(color = "black", face = "bold"), panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

combined_plot <- percent_plot/counts_plot + plot_layout(ncol = 1, heights = c(1, 1.05), guides = "collect") &
  theme(legend.position = "bottom")

write.csv(summary_table, output_file("Pig_vs_Environment_UniqueDetections_ByDate.csv"), row.names = FALSE)

write.xlsx(list(Long_Summary = summary_table, Wide_Summary = summary_wide, Overlap_Summary = overlap_summary),
  file = output_file("Pig_vs_Environment_UniqueDetections_ByDate.xlsx"), overwrite = TRUE)

ggsave(filename = output_file("Pig_vs_Environment_UniqueDetections_Counts.pdf"), plot = counts_plot,
  width = 8.5, height = 4.8, units = "in")

ggsave(filename = output_file("Pig_vs_Environment_UniqueDetections_Percent.pdf"), plot = percent_plot,
  width = 8.5, height = 4.8, units = "in")

ggsave(filename = output_file("Pig_vs_Environment_UniqueDetections_Combined.pdf"), plot = combined_plot,
  width = 8.5, height = 8.5, units = "in")
combined_plot


# 13. Summarise detections unique to each environmental sample type ----

species_path <- input_file("Pig_Species_Detection.xlsx")
metadata_path <- input_file("Pig_Metadata.xlsx")

sample_type_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

species <- read.xlsx(species_path, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID")

metadata <- read.xlsx(metadata_path, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID") %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Date = as.character(Date), Reverse_Date = as.character(Reverse_Date),
    Plot_Date = as.Date(Reverse_Date, format = "%Y_%m_%d"))

detections_long <- species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Present") %>%
  filter(Present > 0) %>%
  left_join(metadata %>%
    select(Sample_ID, Sample_Type, Date, Reverse_Date, Plot_Date), by = "Sample_ID") %>%
  filter(!is.na(Plot_Date), !is.na(Sample_Type), Sample_Type != "")

unique_by_type_date <- detections_long %>%
  distinct(Plot_Date, Date, Virus, Sample_Type)

exclusive_environment_detections <- unique_by_type_date %>%
  group_by(Plot_Date, Date, Virus) %>%
  mutate(N_Sample_Types = n_distinct(Sample_Type)) %>%
  ungroup() %>%
  filter(N_Sample_Types == 1, Sample_Type != "Pig")

summary_table <- exclusive_environment_detections %>%
  count(Sample_Type, name = "Exclusive_Unique_Detections") %>%
  mutate(Percentage_of_Exclusive_Environmental_Detections = 100 * Exclusive_Unique_Detections/sum(Exclusive_Unique_Detections)) %>%
  arrange(desc(Exclusive_Unique_Detections)) %>%
  mutate(Sample_Type = factor(Sample_Type, levels = Sample_Type))

date_levels <- exclusive_environment_detections %>%
  distinct(Plot_Date) %>%
  arrange(Plot_Date) %>%
  mutate(Date_Label = format(Plot_Date, "%d %b %Y")) %>%
  pull(Date_Label)

summary_by_date <- exclusive_environment_detections %>%
  mutate(Date_Label = format(Plot_Date, "%d %b %Y")) %>%
  count(Plot_Date, Date_Label, Sample_Type, name = "Exclusive_Unique_Detections") %>%
  mutate(Date_Label = factor(Date_Label, levels = date_levels), Sample_Type = factor(Sample_Type,
    levels = levels(summary_table$Sample_Type))) %>%
  arrange(Plot_Date, Sample_Type)

detail_table <- exclusive_environment_detections %>%
  mutate(Date_Label = format(Plot_Date, "%d %b %Y")) %>%
  select(Date_Label, Sample_Type, Virus) %>%
  arrange(Sample_Type, Date_Label, Virus)

bar_plot <- ggplot(summary_table, aes(x = Sample_Type, y = Exclusive_Unique_Detections, fill = Sample_Type)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.25) + geom_text(aes(label = Exclusive_Unique_Detections),
  vjust = -0.35, size = 3.8, fontface = "bold") + scale_fill_manual(values = sample_type_palette[names(sample_type_palette) !=
  "Pig"]) + expand_limits(y = max(summary_table$Exclusive_Unique_Detections) * 1.12) + labs(x = NULL,
  y = "Exclusive unique detections", fill = NULL, title = "Environmental sample-type-exclusive detections not found in pig") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1,
  vjust = 1, color = "black", face = "bold"), axis.text.y = element_text(color = "black",
  face = "bold"), axis.title.y = element_text(color = "black", face = "bold"), panel.grid.minor = element_blank(),
  legend.position = "none", plot.title = element_text(face = "bold"))

date_plot <- ggplot(summary_by_date, aes(x = Date_Label, y = Exclusive_Unique_Detections, fill = Sample_Type)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.8) + scale_fill_manual(values = sample_type_palette[names(sample_type_palette) !=
  "Pig"]) + labs(x = NULL, y = "Exclusive unique detections", fill = NULL, title = "Environmental sample-type-exclusive detections by sampling date") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1,
  vjust = 1, color = "black", face = "bold"), axis.text.y = element_text(color = "black",
  face = "bold"), axis.title.y = element_text(color = "black", face = "bold"), panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

combined_plot <- bar_plot/date_plot + plot_layout(ncol = 1, heights = c(1, 1.1))

write.csv(summary_table, output_file("Environment_Exclusive_UniqueDetections_BySampleType.csv"),
  row.names = FALSE)

write.xlsx(list(Summary = summary_table, By_Date = summary_by_date, Details = detail_table),
  file = output_file("Environment_Exclusive_UniqueDetections_BySampleType.xlsx"), overwrite = TRUE)

ggsave(filename = output_file("Environment_Exclusive_UniqueDetections_BySampleType.pdf"), plot = bar_plot,
  width = 8, height = 5, units = "in")

ggsave(filename = output_file("Environment_Exclusive_UniqueDetections_BySampleType_Combined.pdf"),
  plot = combined_plot, width = 5, height = 9, units = "in")

combined_plot


# 14. Plot detection outcomes with a Sankey diagram ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

trim_char_cols <- function(df) {
  df %>%
    mutate(across(where(is.character), ~trimws(as.character(.x))))
}
ensure_col <- function(df, target, candidates) {
  if (target %in% names(df))
    return(df)
  hit <- intersect(candidates, names(df))
  if (length(hit) >= 1)
    return(df %>%
      rename(!!target := all_of(hit[1])))
  df[[target]] <- NA
  df
}

Metadata <- Metadata %>%
  as_tibble(rownames = "Sample_ID") %>%
  trim_char_cols() %>%
  ensure_col("Sample_Type", c("Sample_Type", "Sample Group", "Sample_Group", "sample_type",
    "sample.group", "SampleType", "Group", "Type")) %>%
  ensure_col("Experiment", c("Experiment", "Exp", "Experiment_ID", "ExperimentId", "experiment")) %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Experiment = trimws(as.character(Experiment))) %>%
  filter(!is.na(Sample_ID) & Sample_ID != "")

Species <- Species %>%
  as_tibble(rownames = "Sample_ID")
Species_long <- Species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "value") %>%
  mutate(Present = as.integer(replace_na(value, 0) > 0)) %>%
  filter(Present == 1L) %>%
  select(Sample_ID, Virus)

df <- left_join(Species_long, Metadata %>%
  select(Sample_ID, Sample_Type, Experiment), by = "Sample_ID") %>%
  trim_char_cols() %>%
  ensure_col("Sample_Type", c("Sample_Type", "Sample Group", "Sample_Group", "sample_type",
    "sample.group", "SampleType", "Group", "Type")) %>%
  ensure_col("Experiment", c("Experiment", "Exp", "Experiment_ID", "ExperimentId", "experiment")) %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Experiment = trimws(as.character(Experiment))) %>%
  filter(!is.na(Sample_Type), !is.na(Experiment))

if (nrow(df) == 0) {
  present_by_type <- tibble(Experiment = character(), Virus = character(), Sample_Type = character(),
    present = integer())
} else {
  present_by_type <- df %>%
    distinct(Experiment, Virus, Sample_Type) %>%
    mutate(present = 1L)
}

wide <- present_by_type %>%
  tidyr::pivot_wider(names_from = Sample_Type, values_from = present, values_fill = 0)

all_types <- sort(setdiff(unique(Metadata$Sample_Type), "Pig"))
for (t in all_types) if (!t %in% names(wide)) wide[[t]] <- integer(nrow(wide))
if (!"Pig" %in% names(wide)) wide$Pig <- integer(nrow(wide))

counts <- map_dfr(all_types, function(t) {
  co_det <- sum(wide$Pig == 1 & wide[[t]] == 1, na.rm = TRUE)
  pig_only <- sum(wide$Pig == 1 & wide[[t]] == 0, na.rm = TRUE)
  other_only <- sum(wide$Pig == 0 & wide[[t]] == 1, na.rm = TRUE)
  tibble(type = t, co_det = co_det, pig_only = pig_only, other_only = other_only)
})

stages <- c("Other (L)", "Pig (M)", "Other (R)", "Pig (R)")

mk_flow <- function(type, outcome, v_co, v_pig, v_other) {
  if (outcome == "Co_detection") {
    tibble(flow_id = paste(type, outcome, sep = "::"), stage = stages, node = c(paste0(type,
      "_L"), "Pig_M", paste0(type, "_R"), "Pig_R"), value = c(v_co, v_co, 0, 0), Outcome = outcome)
  } else if (outcome == "Pig_only") {
    tibble(flow_id = paste(type, outcome, sep = "::"), stage = stages, node = c(paste0(type,
      "_L"), "Pig_M", paste0(type, "_R"), "Pig_R"), value = c(0, v_pig, v_pig, 0), Outcome = outcome)
  } else {
    tibble(flow_id = paste(type, outcome, sep = "::"), stage = stages, node = c(paste0(type,
      "_L"), "Pig_M", paste0(type, "_R"), "Pig_R"), value = c(0, 0, v_other, v_other),
      Outcome = outcome)
  }
}

alluvial_long <- counts %>%
  pmap_dfr(function(type, co_det, pig_only, other_only) {
    bind_rows(mk_flow(type, "Co_detection", co_det, pig_only, other_only), mk_flow(type,
      "Pig_only", co_det, pig_only, other_only), mk_flow(type, "Other_only", co_det,
      pig_only, other_only))
  }) %>%
  group_by(flow_id) %>%
  filter(sum(value) > 0) %>%
  ungroup()

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

status_cols <- c(Co_detection = "#F8766D", Pig_only = "#F6E09B", Other_only = "#00BFC4")

stratum_labels <- function(x) gsub("(_L|_R|_M)$", "", x)

strata_names <- c("Pig_M", "Pig_R", paste0(all_types, "_L"), paste0(all_types, "_R"))
strata_colors <- setNames(c(rep(group_palette["Pig"], 2), group_palette[all_types], group_palette[all_types]),
  strata_names)

GUTTER_COL <- "white"
RIBBON_EDGE <- 1
STRATUM_EDGE <- 1

Sankey_3stage <- ggplot(alluvial_long, aes(x = stage, stratum = node, alluvium = flow_id, y = value)) +
  geom_alluvium(aes(fill = Outcome), width = 0.15, alpha = 1, color = GUTTER_COL, size = RIBBON_EDGE,
    lineend = "round") + scale_fill_manual(values = status_cols, name = "Outcome") + ggnewscale::new_scale_fill() +

geom_stratum(aes(fill = after_stat(stratum)), width = 0.3, color = GUTTER_COL, linewidth = STRATUM_EDGE) +
  scale_fill_manual(values = strata_colors, guide = "none") + 
geom_text(stat = "stratum", aes(label = stratum_labels(after_stat(stratum))), size = 2.9, fontface = "bold") +

scale_x_discrete(limits = stages, expand = c(0.08, 0.08)) + labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 13) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
  panel.grid = element_blank(), plot.title = element_blank(), legend.title = element_text(face = "bold"),
  legend.position = "right", plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

print(Sankey_3stage)
# ggsave(output_file('Sankey_Plot_Detection_Rates.pdf'), Sankey_3stage, width = 250, height = 250, units = 'mm')


# 15. Plot co-detections with a chord diagram ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Metadata <- Metadata %>%
  rownames_to_column("Sample_ID")
Species <- Species %>%
  rownames_to_column("Sample_ID")

if (!"Experiment" %in% names(Metadata)) {
  Metadata <- Metadata %>%
    mutate(Experiment = "All")
}
if (!"Sample_Type" %in% names(Metadata) && "Sample.Type" %in% names(Metadata)) {
  Metadata <- Metadata %>%
    rename(Sample_Type = Sample.Type)
}

meta_sel <- Metadata %>%
  mutate(Sample_Type = trimws(Sample_Type)) %>%
  select(Sample_ID, Sample_Type, Experiment)

Species_long <- Species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Present") %>%
  filter(Present == 1) %>%
  select(-Present)

df <- left_join(Species_long, meta_sel, by = "Sample_ID")

if (!"Sample_Type" %in% names(df)) {
  if ("Sample.Type" %in% names(df)) {
    df <- df %>%
      rename(Sample_Type = Sample.Type)
  } else {
    stop("No `Sample_Type` column after join. Check Metadata column names and Sample_ID matching.")
  }
}

df <- df %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Experiment = if ("Experiment" %in%
    names(.)) as.character(Experiment) else "All") %>%
  filter(!is.na(Sample_Type) & Sample_Type != "")

present_by_type <- df %>%
  distinct(Experiment, Virus, Sample_Type) %>%
  mutate(present = 1L)

wide <- present_by_type %>%
  pivot_wider(names_from = Sample_Type, values_from = present, values_fill = 0)

all_types <- names(group_palette)
for (t in all_types) {
  if (!t %in% names(wide))
    wide[[t]] <- 0L
}

type_pairs <- t(combn(all_types, 2)) %>%
  as.data.frame()
names(type_pairs) <- c("A", "B")

M <- as.matrix(wide[all_types])
storage.mode(M) <- "integer"

co <- t(M) %*% M

idx <- which(upper.tri(co), arr.ind = TRUE)
pair_counts <- tibble::tibble(A = rownames(co)[idx[, 1]], B = colnames(co)[idx[, 2]], value = as.integer(co[idx])) %>%
  dplyr::filter(value > 0)

if (nrow(pair_counts) == 0) stop("No pairwise co-detections found across Sample_Types.")

dat <- pair_counts %>%
  dplyr::transmute(from = A, to = B, value = value)

if (nrow(pair_counts) == 0) {
  stop("No pairwise co-detections found across Sample_Types.")
}

dat <- pair_counts %>%
  transmute(from = A, to = B, value)

sectors <- all_types

grid.col <- group_palette[sectors]

rng <- range(dat$value, na.rm = TRUE)
alpha_lo <- 0.25
alpha_hi <- 0.95
alpha_vals <- if (diff(rng) == 0) rep((alpha_lo + alpha_hi)/2, nrow(dat)) else (dat$value -
  rng[1])/diff(rng) * (alpha_hi - alpha_lo) + alpha_lo

from_cols <- group_palette[as.character(dat$from)]
link_col <- vapply(seq_len(nrow(dat)), function(i) {
  rgb(t(col2rgb(from_cols[i]))/255, alpha = alpha_vals[i])
}, character(1))

lwd_lo <- 1
lwd_hi <- 6
link_lwd <- if (diff(rng) == 0) rep((lwd_lo + lwd_hi)/2, nrow(dat)) else (dat$value - rng[1])/diff(rng) *
  (lwd_hi - lwd_lo) + lwd_lo

pdf(output_file("Circos_Codetections.pdf"), width = 8, height = 11.5)

circos.clear()
circos.par(gap.after = rep(6, length(sectors)), track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)

chordDiagram(x = dat[, c("from", "to", "value")], order = sectors, grid.col = grid.col, col = link_col,
  transparency = 0, directional = 0, annotationTrack = "none", link.lwd = link_lwd, link.lty = 1,
  link.border = NA, preAllocateTracks = list(track.height = mm_h(6)))

circos.trackPlotRegion(track.index = 1, bg.col = grid.col, bg.border = NA, panel.fun = function(x,
  y) {
  sector <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[2] + mm_y(2.5), sector, facing = "bending.outside", cex = 0.85,
    font = 2)
})
dev.off()


# 16. Plot unique detections with radar charts ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Metadata <- Metadata %>%
  rownames_to_column("Sample_ID")
Species <- Species %>%
  rownames_to_column("Sample_ID")

if (!"Experiment" %in% names(Metadata)) {
  Metadata <- Metadata %>%
    mutate(Experiment = "All")
}
if (!"Sample_Type" %in% names(Metadata) && "Sample.Type" %in% names(Metadata)) {
  Metadata <- Metadata %>%
    rename(Sample_Type = Sample.Type)
}

meta_sel <- Metadata %>%
  mutate(Sample_Type = trimws(Sample_Type)) %>%
  select(Sample_ID, Sample_Type, Experiment)

Species_long <- Species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Present") %>%
  filter(Present == 1) %>%
  select(-Present)

df <- left_join(Species_long, meta_sel, by = "Sample_ID")

if (!"Sample_Type" %in% names(df)) {
  if ("Sample.Type" %in% names(df)) {
    df <- df %>%
      rename(Sample_Type = Sample.Type)
  } else {
    stop("No `Sample_Type` column after join. Check Metadata column names and Sample_ID matching.")
  }
}

df <- df %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Experiment = if ("Experiment" %in%
    names(.)) as.character(Experiment) else "All") %>%
  filter(!is.na(Sample_Type) & Sample_Type != "")

present_by_type <- df %>%
  distinct(Experiment, Virus, Sample_Type) %>%
  mutate(present = 1L)

wide <- present_by_type %>%
  pivot_wider(names_from = Sample_Type, values_from = present, values_fill = 0)

all_types <- names(group_palette)
for (t in all_types) {
  if (!t %in% names(wide))
    wide[[t]] <- 0L
}

type_pairs <- t(combn(all_types, 2)) %>%
  as.data.frame()
names(type_pairs) <- c("A", "B")

M <- as.matrix(wide[all_types])
storage.mode(M) <- "integer"

co <- t(M) %*% M

create_beautiful_radarchart <- function(data, color = "#00AFBB", vlabels = colnames(data),
  vlcex = 0.7, caxislabels = NULL, title = NULL, ...) {
  radarchart(data, axistype = 1, pcol = color, pfcol = scales::alpha(color, 0.2), plwd = 2,
    plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8, axislabcol = "grey", vlcex = vlcex,
    vlabels = vlabels, caxislabels = caxislabels, title = title, ...)
}

all_types <- names(group_palette)

co_sub <- co[all_types, all_types, drop = FALSE]

totals <- colSums(M)[all_types]

radar_mat <- co_sub
diag(radar_mat) <- totals[rownames(radar_mat)]

max_val <- max(radar_mat)
min_val <- 0

par(mfrow = c(3, 3), mar = c(1, 2, 2, 1))

for (focal in all_types) {
  vars <- c(focal, setdiff(all_types, focal))
  row_vals <- radar_mat[focal, vars]
  df_plot <- rbind(rep(max_val, length(vars)), rep(min_val, length(vars)), row_vals)
  df_plot <- as.data.frame(df_plot)
  colnames(df_plot) <- vars
  cax <- pretty(c(min_val, max_val), n = 5)
  create_beautiful_radarchart(data = df_plot, color = group_palette[focal], vlabels = vars,
    vlcex = 0.7, caxislabels = cax, title = focal)
}


# 17. Compare alpha diversity ----
Species <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- t(Species_bin)
Species_bin <- as.data.frame(Species_bin)

common_ids <- intersect(colnames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[, common_ids, drop = FALSE]
Metadata <- Metadata[common_ids, , drop = FALSE]

ASV_Alpha <- otu_table(Species_bin, taxa_are_rows = TRUE)
Meta_ps <- sample_data(Metadata)
physeq <- phyloseq(ASV_Alpha, Meta_ps)

Alpha_Diversity <- estimate_richness(physeq, split = TRUE, measures = c("Observed", "Shannon",
  "Simpson", "Chao1"))

Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))

if (!"Sample_Type" %in% names(Alpha_Diversity)) {
  stop("Sample_Type column not found in Metadata.")
}
Alpha_Diversity$Sample_Type <- as.factor(Alpha_Diversity$Sample_Type)

desired_order <- names(group_palette)

present_groups <- desired_order[desired_order %in% Alpha_Diversity$Sample_Type]

Alpha_Diversity$Sample_Type <- factor(Alpha_Diversity$Sample_Type, levels = present_groups)

plot_palette <- group_palette[present_groups]

Kruskall_Observed <- kwAllPairsDunnTest(Alpha_Diversity$Observed, Alpha_Diversity$Sample_Type,
  p.adjust.method = "fdr")
Kruskall_Chao1 <- kwAllPairsDunnTest(Alpha_Diversity$Chao1, Alpha_Diversity$Sample_Type, p.adjust.method = "fdr")
Kruskall_Shannon <- kwAllPairsDunnTest(Alpha_Diversity$Shannon, Alpha_Diversity$Sample_Type,
  p.adjust.method = "fdr")
Kruskall_Simpson <- kwAllPairsDunnTest(Alpha_Diversity$Simpson, Alpha_Diversity$Sample_Type,
  p.adjust.method = "fdr")

# write.xlsx(Kruskall_Observed$p.value, output_file('Supplementary_Kruskall_Alpha_Observed_Pig.xlsx'), rowNames = TRUE)
# write.xlsx(Kruskall_Chao1$p.value,    output_file('Supplementary_Kruskall_Alpha_Chao1_Pig.xlsx'),    rowNames = TRUE)
# write.xlsx(Kruskall_Shannon$p.value,  output_file('Supplementary_Kruskall_Alpha_Shannon_Pig.xlsx'),  rowNames = TRUE)
# write.xlsx(Kruskall_Simpson$p.value,  output_file('Supplementary_Kruskall_Alpha_Simpson_Pig.xlsx'),  rowNames = TRUE)

alpha_long <- pivot_longer(Alpha_Diversity, cols = c("Observed", "Shannon", "Simpson", "Chao1"),
  names_to = "Metric", values_to = "Value")

alpha_long$Sample_Type <- factor(alpha_long$Sample_Type, levels = present_groups)

plot_list <- alpha_long %>%
  split(.$Metric) %>%
  map(~ggplot(.x, aes(x = Sample_Type, y = Value, fill = Sample_Type)) + geom_boxplot(outlier.shape = NA,
    color = "black") + geom_jitter(aes(color = Sample_Type), width = 0.2, alpha = 0.6,
    size = 1.5) + scale_fill_manual(values = plot_palette) + scale_color_manual(values = plot_palette) +
    theme_minimal() + labs(title = unique(.x$Metric), x = NULL, y = "Alpha Diversity") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
      axis.line.x = element_line(color = "grey50"), axis.line.y = element_line(color = "grey50"),
      axis.text.y = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"), strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, size = 1), legend.position = "bottom"))

final_alpha_plot <- plot_list$Chao1 + plot_list$Shannon + plot_layout(nrow = 2)

final_alpha_plot
# ggsave(filename = output_file('Pig_Alpha_Diversity_Boxplots.pdf'), plot = final_alpha_plot, width = 200, height = 280, units = 'mm')
# ggsave(filename = output_file('Pig_Observed_Species_Boxplots.pdf'), plot = plot_list$Observed, width = 160, height = 160, units = 'mm')


# 18. Compare beta diversity using Bray-Curtis dissimilarity ----
Species <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species <- Species[rowSums(Species) != 0, , drop = FALSE]
Metadata <- Metadata[rownames(Species), , drop = FALSE]

OTU <- otu_table(t(Species), taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
comm_mat <- t(as(otu_table(physeq), "matrix"))

Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
pal <- group_palette[present_groups]

Bray <- vegdist(comm_mat, method = "bray")
Bray <- as.dist(Bray)

Bray_PCoA <- dudi.pco(Bray, scannf = FALSE, nf = 3)
Bray_Li <- Bray_PCoA$li
Bray_Li$Group <- Metadata[rownames(Bray_Li), "Sample_Type", drop = TRUE]

Bray_Adonis <- adonis2(Bray ~ Metadata$Sample_Type, permutations = 1000)
p_overall <- Bray_Adonis$`Pr(>F)`[1]
r2_overall <- Bray_Adonis$R2[1]

axis1_perc <- Bray_PCoA$eig[1]/sum(Bray_PCoA$eig) * 100
axis2_perc <- Bray_PCoA$eig[2]/sum(Bray_PCoA$eig) * 100

.run_pairwise <- function(dist_obj, meta, group_col, pair, nperm = 999) {
  rows <- which(meta[[group_col]] %in% pair)
  if (length(rows) < 3) {
    return(tibble(group1 = pair[1], group2 = pair[2], n1 = sum(meta[[group_col]] == pair[1]),
      n2 = sum(meta[[group_col]] == pair[2]), N = length(rows), Df = NA_real_, F = NA_real_,
      R2 = NA_real_, p = NA_real_, error = "Too few samples"))
  }
  dsub <- as.dist(as.matrix(dist_obj)[rows, rows])
  msub <- droplevels(meta[rows, , drop = FALSE])
  a2 <- adonis2(dsub ~ msub[[group_col]], permutations = nperm)
  tibble(group1 = pair[1], group2 = pair[2], n1 = sum(msub[[group_col]] == pair[1]), n2 = sum(msub[[group_col]] ==
    pair[2]), N = nrow(msub), Df = as.numeric(a2$Df[1]), F = as.numeric(a2$F[1]), R2 = as.numeric(a2$R2[1]),
    p = as.numeric(a2$`Pr(>F)`[1]), error = NA_character_)
}

run_all_pairwise <- function(dist_obj, meta, group_col = "Sample_Type", palette = pal, nperm = 999) {
  present <- intersect(names(palette), unique(as.character(meta[[group_col]])))
  pairs <- combn(present, 2, simplify = FALSE)
  map_dfr(pairs, ~.run_pairwise(dist_obj, meta, group_col, .x, nperm = nperm)) %>%
    mutate(comparison = paste(group1, "vs", group2), p_adj = p.adjust(p, method = "BH"),
      sig = case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
    select(comparison, group1, group2, n1, n2, N, Df, F, R2, p, p_adj, sig, error) %>%
    arrange(p_adj, p)
}

pairwise_table <- run_all_pairwise(dist_obj = Bray, meta = Metadata, group_col = "Sample_Type",
  palette = pal, nperm = 999)
# openxlsx::write.xlsx(pairwise_table, output_file('Pairwise_PERMANOVA_Bray_Pig.xlsx'), rowNames = TRUE)

centroids <- Bray_Li %>%
  group_by(Group) %>%
  summarise(Centroid_X = mean(A1, na.rm = TRUE), Centroid_Y = mean(A2, na.rm = TRUE), .groups = "drop")

Bray_Li2 <- Bray_Li %>%
  left_join(centroids, by = "Group") %>%
  mutate(Group = factor(Group, levels = names(pal)))

subtitle_txt <- sprintf("P = %.4f%s, R² = %.2f", p_overall, if (p_overall < 0.001) "***" else if (p_overall <
  0.01) "**" else if (p_overall < 0.05) "*" else "", r2_overall)

PCoA_Bray <- ggplot(Bray_Li2, aes(x = A1, y = A2, colour = Group)) + geom_point(size = 2) +
  stat_ellipse(aes(group = Group, fill = Group, colour = Group), type = "t", level = 0.1,
    geom = "polygon", alpha = 0.5, linewidth = 2) + coord_fixed() + theme_classic() + scale_colour_manual(values = pal,
  drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) + labs(title = "Bray–Curtis dissimilarity",
  subtitle = subtitle_txt, x = sprintf("PCo1 (%.2f%%)", axis1_perc), y = sprintf("PCo2 (%.2f%%)",
    axis2_perc)) + theme(panel.background = element_blank(), plot.title = element_text(hjust = 0,
  face = "bold", size = 15), plot.subtitle = element_text(hjust = 0, face = "bold", size = 10),
  axis.line = element_line(size = 0.6, colour = "black"), axis.title = element_text(face = "bold",
    size = 12), axis.text = element_text(face = "bold", size = 12, colour = "black"), legend.position = "bottom")

print(PCoA_Bray)
# ggsave(filename=output_file('Pig_PCoA_Bray_Curtis.pdf'), plot = PCoA_Bray, width=200, height=200, units='mm')

PC_scores <- data.frame(PC1 = Bray_PCoA$li$A1, PC2 = Bray_PCoA$li$A2, Group = Metadata[rownames(Bray_PCoA$li),
  "Sample_Type", drop = TRUE])
PC_scores$Group <- factor(PC_scores$Group, levels = names(pal))

pc1_lim <- range(PC_scores$PC1, na.rm = TRUE) + c(-1, 1) * 0.05 * diff(range(PC_scores$PC1,
  na.rm = TRUE))
pc2_lim <- range(PC_scores$PC2, na.rm = TRUE) + c(-1, 1) * 0.05 * diff(range(PC_scores$PC2,
  na.rm = TRUE))

PC1_Box <- ggplot(PC_scores, aes(x = Group, y = PC1)) + geom_boxplot(aes(fill = Group), alpha = 0.9,
  outlier.shape = NA, width = 0.6) + geom_jitter(aes(colour = Group), width = 0.15, size = 1.5,
  alpha = 0.6) + scale_fill_manual(values = pal, name = "Sample Type") + scale_colour_manual(values = pal,
  name = "Sample Type") + scale_y_continuous(limits = pc1_lim) + labs(x = NULL, y = "PCo1") +
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(colour = "grey50", fill = NA),
  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.x = element_blank(), legend.position = "bottom")

PC2_Box <- ggplot(PC_scores, aes(x = Group, y = PC2)) + geom_boxplot(aes(fill = Group), alpha = 0.9,
  outlier.shape = NA, width = 0.6) + geom_jitter(aes(colour = Group), width = 0.15, size = 1.5,
  alpha = 0.6) + scale_fill_manual(values = pal, name = "Sample Type") + scale_colour_manual(values = pal,
  name = "Sample Type") + scale_y_continuous(limits = pc2_lim) + labs(x = NULL, y = "PCo2") +
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(colour = "grey50", fill = NA),
  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.x = element_blank(), legend.position = "bottom")

PC1_Dunn <- kwAllPairsDunnTest(PC_scores$PC1, PC_scores$Group, p.adjust.method = "fdr")
PC2_Dunn <- kwAllPairsDunnTest(PC_scores$PC2, PC_scores$Group, p.adjust.method = "fdr")
# openxlsx::write.xlsx(PC1_Dunn$p.value, output_file('Pig_PC1_Dunn.xlsx'), rowNames = TRUE)
# openxlsx::write.xlsx(PC2_Dunn$p.value, output_file('Pig_PC2_Dunn.xlsx'), rowNames = TRUE)

final_plot <- PCoA_Bray/(PC1_Box | PC2_Box) + plot_layout(heights = c(2, 1))
print(final_plot)
# ggsave(output_file('Pig_PCoA_Bray_Panels.pdf'), final_plot, width = 210, height = 290, units = 'mm')


# 19. Compare beta diversity using Jaccard distance ----
Species <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]

OTU <- otu_table(t(Species_bin), taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
comm_mat <- t(as(otu_table(physeq), "matrix"))

Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
pal <- group_palette[present_groups]

Jaccard <- vegdist(comm_mat, method = "jaccard", binary = TRUE)
Jaccard <- as.dist(Jaccard)

Jaccard_PCoA <- dudi.pco(Jaccard, scannf = FALSE, nf = 3)
Jac_Li <- Jaccard_PCoA$li
Jac_Li$Group <- Metadata[rownames(Jac_Li), "Sample_Type", drop = TRUE]

Jac_Adonis <- adonis2(Jaccard ~ Metadata$Sample_Type, permutations = 1000)
p_overall <- Jac_Adonis$`Pr(>F)`[1]
r2_overall <- Jac_Adonis$R2[1]

axis1_perc <- Jaccard_PCoA$eig[1]/sum(Jaccard_PCoA$eig) * 100
axis2_perc <- Jaccard_PCoA$eig[2]/sum(Jaccard_PCoA$eig) * 100

.run_pairwise <- function(dist_obj, meta, group_col, pair, nperm = 999) {
  rows <- which(meta[[group_col]] %in% pair)
  if (length(rows) < 3) {
    return(tibble(group1 = pair[1], group2 = pair[2], n1 = sum(meta[[group_col]] == pair[1]),
      n2 = sum(meta[[group_col]] == pair[2]), N = length(rows), Df = NA_real_, F = NA_real_,
      R2 = NA_real_, p = NA_real_, error = "Too few samples"))
  }
  dsub <- as.dist(as.matrix(dist_obj)[rows, rows])
  msub <- droplevels(meta[rows, , drop = FALSE])
  a2 <- adonis2(dsub ~ msub[[group_col]], permutations = nperm)
  tibble(group1 = pair[1], group2 = pair[2], n1 = sum(msub[[group_col]] == pair[1]), n2 = sum(msub[[group_col]] ==
    pair[2]), N = nrow(msub), Df = as.numeric(a2$Df[1]), F = as.numeric(a2$F[1]), R2 = as.numeric(a2$R2[1]),
    p = as.numeric(a2$`Pr(>F)`[1]), error = NA_character_)
}

run_all_pairwise <- function(dist_obj, meta, group_col = "Sample_Type", palette = pal, nperm = 999) {
  present <- intersect(names(palette), unique(as.character(meta[[group_col]])))
  pairs <- combn(present, 2, simplify = FALSE)
  map_dfr(pairs, ~.run_pairwise(dist_obj, meta, group_col, .x, nperm = nperm)) %>%
    mutate(comparison = paste(group1, "vs", group2), p_adj = p.adjust(p, method = "BH"),
      sig = case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
    select(comparison, group1, group2, n1, n2, N, Df, F, R2, p, p_adj, sig, error) %>%
    arrange(p_adj, p)
}

pairwise_table <- run_all_pairwise(dist_obj = Jaccard, meta = Metadata, group_col = "Sample_Type",
  palette = pal, nperm = 999)
# write.xlsx(pairwise_table, output_file('Supplementary_Pig_Pairwise_PERMANOVA_Jaccard.xlsx'), rowNames = TRUE)

Jac_Li <- as.data.frame(Jaccard_PCoA$li)
Jac_Li$Group <- factor(Metadata[rownames(Jac_Li), "Sample_Type", drop = TRUE], levels = names(pal))

centroids <- dplyr::group_by(Jac_Li, Group) |>
  dplyr::summarise(Centroid_X = mean(A1, na.rm = TRUE), Centroid_Y = mean(A2, na.rm = TRUE),
    .groups = "drop")

Jac_Li2 <- dplyr::left_join(Jac_Li, centroids, by = "Group") |>
  dplyr::mutate(Group = factor(Group, levels = names(pal)))

subtitle_txt <- sprintf("P = %.4f%s, R² = %.2f", p_overall, if (p_overall < 0.001) "***" else if (p_overall <
  0.01) "**" else if (p_overall < 0.05) "*" else "", r2_overall)

PCoA_Jaccard <- ggplot(Jac_Li2, aes(x = A1, y = A2, colour = Group)) + geom_point(size = 2) +
  stat_ellipse(aes(group = Group, fill = Group, colour = Group), type = "t", level = 0.1,
    geom = "polygon", alpha = 0.5, linewidth = 2) + coord_fixed() + theme_classic() + scale_colour_manual(values = pal,
  drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) + labs(title = "Jaccard dissimilarity",
  subtitle = subtitle_txt, x = sprintf("PCo1 (%.2f%%)", axis1_perc), y = sprintf("PCo2 (%.2f%%)",
    axis2_perc)) + theme(panel.background = element_blank(), plot.title = element_text(hjust = 0,
  face = "bold", size = 15), plot.subtitle = element_text(hjust = 0, face = "bold", size = 10),
  axis.line = element_line(size = 0.6, colour = "black"), axis.title = element_text(face = "bold",
    size = 12), axis.text = element_text(face = "bold", size = 12, colour = "black"), legend.position = "bottom")

print(PCoA_Jaccard)
# ggsave(filename=output_file('Pig_PCoA_Jaccard.pdf'), plot = PCoA_Jaccard, width=200, height=200, units='mm')

PC_scores <- data.frame(PC1 = Jaccard_PCoA$li$A1, PC2 = Jaccard_PCoA$li$A2, Group = Metadata[rownames(Jaccard_PCoA$li),
  "Sample_Type", drop = TRUE])
PC_scores$Group <- factor(PC_scores$Group, levels = names(pal))

desired_order <- names(group_palette)
present <- desired_order[desired_order %in% PC_scores$Group]
pal <- group_palette[present]
PC_scores$Group <- factor(PC_scores$Group, levels = present)

pc1_lim <- range(PC_scores$PC1, na.rm = TRUE) + c(-1, 1) * 0.05 * diff(range(PC_scores$PC1,
  na.rm = TRUE))
pc2_lim <- range(PC_scores$PC2, na.rm = TRUE) + c(-1, 1) * 0.05 * diff(range(PC_scores$PC2,
  na.rm = TRUE))

PC1_Box <- ggplot(PC_scores, aes(x = Group, y = PC1)) + geom_boxplot(aes(fill = Group), alpha = 0.9,
  outlier.shape = NA, width = 0.6) + geom_jitter(aes(colour = Group), width = 0.15, size = 1.5,
  alpha = 0.6) + scale_fill_manual(values = pal, name = "Sample Type") + scale_colour_manual(values = pal,
  name = "Sample Type") + scale_y_continuous(limits = pc1_lim) + labs(x = NULL, y = "PCo1") +
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(colour = "grey50", fill = NA),
  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.x = element_blank(), legend.position = "bottom")

PC2_Box <- ggplot(PC_scores, aes(x = Group, y = PC2)) + geom_boxplot(aes(fill = Group), alpha = 0.9,
  outlier.shape = NA, width = 0.6) + geom_jitter(aes(colour = Group), width = 0.15, size = 1.5,
  alpha = 0.6) + scale_fill_manual(values = pal, name = "Sample Type") + scale_colour_manual(values = pal,
  name = "Sample Type") + scale_y_continuous(limits = pc2_lim) + labs(x = NULL, y = "PCo2") +
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(colour = "grey50", fill = NA),
  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.x = element_blank(), legend.position = "bottom")

PC1_Dunn <- kwAllPairsDunnTest(PC_scores$PC1, PC_scores$Group, p.adjust.method = "fdr")
PC2_Dunn <- kwAllPairsDunnTest(PC_scores$PC2, PC_scores$Group, p.adjust.method = "fdr")
# openxlsx::write.xlsx(PC1_Dunn$p.value, output_file('Pig_PC1_Dunn_Jaccard.xlsx'), rowNames = TRUE)
# openxlsx::write.xlsx(PC2_Dunn$p.value, output_file('Pig_PC2_Dunn_Jaccard.xlsx'), rowNames = TRUE)

final_plot <- PCoA_Jaccard/(PC1_Box | PC2_Box) + plot_layout(heights = c(2, 1))
print(final_plot)
# ggsave(output_file('Pig_PCoA_Jaccard_Panels.pdf'), final_plot, width = 210, height = 290, units = 'mm')


# 20. Compare paired pig and environmental samples using Jaccard distance ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]

OTU <- otu_table(t(Species_bin), taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
comm_mat <- t(as(otu_table(physeq), "matrix"))

Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
pal <- group_palette[present_groups]

Jaccard <- vegdist(comm_mat, method = "jaccard", binary = TRUE)
Jaccard <- as.dist(Jaccard)

Jaccard_PCoA <- dudi.pco(Jaccard, scannf = FALSE, nf = 3)

if (!"Experiment" %in% names(Metadata)) {
  Metadata$Experiment <- "All"
}
Metadata$Experiment <- as.factor(Metadata$Experiment)

desired_order <- names(group_palette)
others_order <- setdiff(desired_order, "Pig")
present_others <- others_order[others_order %in% Metadata$Sample_Type]

D <- as.matrix(Jaccard)
stopifnot(identical(rownames(D), rownames(Metadata)))

experiments <- levels(Metadata$Experiment)
dist_list <- vector("list", length(experiments))
names(dist_list) <- experiments

for (e in experiments) {
  idx_e <- which(Metadata$Experiment == e)
  if (length(idx_e) < 2)
    next

  ids_e <- rownames(Metadata)[idx_e]
  pigs_e <- ids_e[Metadata$Sample_Type[idx_e] == "Pig"]

  if (length(pigs_e) > 1) {
    d_pig <- D[pigs_e, pigs_e, drop = FALSE]
    d_pig <- d_pig[upper.tri(d_pig, diag = FALSE)]
    dist_list[[e]] <- rbind(dist_list[[e]], data.frame(Experiment = e, Sample_Type = "Pig",
      distance = as.vector(d_pig), stringsAsFactors = FALSE))
  }

  if (length(pigs_e) == 0)
    next
  for (t in present_others) {
    others_e <- ids_e[Metadata$Sample_Type[idx_e] == t]
    if (length(others_e) == 0)
      next

    d_mat <- D[pigs_e, others_e, drop = FALSE]
    dist_list[[e]] <- rbind(dist_list[[e]], data.frame(Experiment = e, Sample_Type = t,
      distance = as.vector(d_mat), stringsAsFactors = FALSE))
  }
}

dist_df <- dplyr::bind_rows(dist_list)
if (nrow(dist_df) == 0) stop("No distances computed. Check Experiment/Sample_Type coverage.")

all_groups_in_df <- intersect(names(group_palette), unique(as.character(dist_df$Sample_Type)))
dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = all_groups_in_df)

ord_stats <- dist_df %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::summarise(mean_dist = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Sample_Type = as.character(Sample_Type))

pig_first <- "Pig"

others_ord <- ord_stats %>%
  dplyr::filter(.data$Sample_Type != pig_first) %>%
  dplyr::arrange(mean_dist) %>%
  dplyr::pull(Sample_Type)

ord <- c(pig_first, others_ord)

dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = ord)
pal_ord <- group_palette[ord]

Kruskall_Jaccard <- kwAllPairsDunnTest(dist_df$distance, dist_df$Sample_Type, p.adjust.method = "fdr")
# openxlsx::write.xlsx(Kruskall_Jaccard$p.value, output_file('Supplementary_Kruskall_Jaccard_Distance_Including_Pig.xlsx'), rowNames = TRUE)

p_rain <- ggplot(dist_df, aes(Sample_Type, distance, fill = Sample_Type, colour = Sample_Type)) +
  stat_halfeye(adjust = 0.7, width = 0.7, justification = -0.2, point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.9) + geom_jitter(width = 0.12,
  size = 1.1, alpha = 0.45) + scale_fill_manual(values = pal_ord, drop = FALSE) + scale_colour_manual(values = pal_ord,
  drop = FALSE) + scale_x_discrete(limits = ord, drop = FALSE) + scale_y_continuous(limits = c(0.3,
  1), breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.02, 0.05))) + labs(title = expression("Distance Between Paired Samples (Jaccard " *
  beta * "-diversity)"), subtitle = "", x = NULL, y = "Jaccard Distance") + theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(colour = "grey60", fill = NA), axis.text.x = element_text(angle = 45,
    hjust = 1, face = "bold", colour = "black"), axis.text.y = element_text(colour = "black"),
    legend.position = "bottom")

print(p_rain)
# ggsave(output_file('Pig_Jaccard_Distance_to_Pig_ALL_Including_Pig.pdf'), p_rain, width = 200, height = 140, units = 'mm')


# 21. Compare paired pig and environmental samples using Bray-Curtis dissimilarity ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]

OTU <- otu_table(t(Species_bin), taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
comm_mat <- t(as(otu_table(physeq), "matrix"))

Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
pal <- group_palette[present_groups]

Bray <- vegdist(comm_mat, method = "bray", binary = TRUE)
Bray <- as.dist(Bray)

Bray_PCoA <- dudi.pco(Bray, scannf = FALSE, nf = 3)

if (!"Experiment" %in% names(Metadata)) {
  Metadata$Experiment <- "All"
}
Metadata$Experiment <- as.factor(Metadata$Experiment)

desired_order <- names(group_palette)
others_order <- setdiff(desired_order, "Pig")
present_others <- others_order[others_order %in% Metadata$Sample_Type]

D <- as.matrix(Bray)
stopifnot(identical(rownames(D), rownames(Metadata)))

experiments <- levels(Metadata$Experiment)
dist_list <- vector("list", length(experiments))
names(dist_list) <- experiments

for (e in experiments) {
  idx_e <- which(Metadata$Experiment == e)
  if (length(idx_e) < 2)
    next

  ids_e <- rownames(Metadata)[idx_e]
  pigs_e <- ids_e[Metadata$Sample_Type[idx_e] == "Pig"]

  if (length(pigs_e) > 1) {
    d_pig <- D[pigs_e, pigs_e, drop = FALSE]
    d_pig <- d_pig[upper.tri(d_pig, diag = FALSE)]
    dist_list[[e]] <- rbind(dist_list[[e]], data.frame(Experiment = e, Sample_Type = "Pig",
      distance = as.vector(d_pig), stringsAsFactors = FALSE))
  }

  if (length(pigs_e) == 0)
    next
  for (t in present_others) {
    others_e <- ids_e[Metadata$Sample_Type[idx_e] == t]
    if (length(others_e) == 0)
      next

    d_mat <- D[pigs_e, others_e, drop = FALSE]
    dist_list[[e]] <- rbind(dist_list[[e]], data.frame(Experiment = e, Sample_Type = t,
      distance = as.vector(d_mat), stringsAsFactors = FALSE))
  }
}

dist_df <- dplyr::bind_rows(dist_list)
if (nrow(dist_df) == 0) stop("No distances computed. Check Experiment/Sample_Type coverage.")

all_groups_in_df <- intersect(names(group_palette), unique(as.character(dist_df$Sample_Type)))
dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = all_groups_in_df)

ord_stats <- dist_df %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::summarise(mean_dist = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Sample_Type = as.character(Sample_Type))

pig_first <- "Pig"

others_ord <- ord_stats %>%
  dplyr::filter(.data$Sample_Type != pig_first) %>%
  dplyr::arrange(mean_dist) %>%
  dplyr::pull(Sample_Type)

ord <- c(pig_first, others_ord)

dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = ord)
pal_ord <- group_palette[ord]

Kruskall_Bray <- kwAllPairsDunnTest(dist_df$distance, dist_df$Sample_Type, p.adjust.method = "fdr")
openxlsx::write.xlsx(Kruskall_Bray$p.value, output_file("Supplementary_Kruskall_Bray_Distance_Including_Pig.xlsx"),
  rowNames = TRUE)

p_rain <- ggplot(dist_df, aes(Sample_Type, distance, fill = Sample_Type, colour = Sample_Type)) +
  stat_halfeye(adjust = 0.7, width = 0.7, justification = -0.2, point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.9) + geom_jitter(width = 0.12,
  size = 1.1, alpha = 0.45) + scale_fill_manual(values = pal_ord, drop = FALSE) + scale_colour_manual(values = pal_ord,
  drop = FALSE) + scale_x_discrete(limits = ord, drop = FALSE) + scale_y_continuous(limits = c(0.3,
  1), breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.02, 0.05))) + labs(title = expression("Distance Between Paired Samples (Bray " *
  beta * "-diversity)"), subtitle = "", x = NULL, y = "Bray Distance") + theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(colour = "grey60", fill = NA), axis.text.x = element_text(angle = 45,
    hjust = 1, face = "bold", colour = "black"), axis.text.y = element_text(colour = "black"),
    legend.position = "right")

print(p_rain)
# ggsave(output_file('Pig_Bray_Distance_to_Pig_ALL_Including_Pig.pdf'), p_rain, width = 200, height = 140, units = 'mm')


# 22. Plot relative abundance by sampling date ----
Species <- read.xlsx(input_file("Pig_Species_RPM.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

tableau20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D",
  "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
  "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")

other_color <- "#D9D9D9"

Species_df <- Species %>%
  rownames_to_column(var = "Sample_ID")

Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_Type = factor(Sample_Type, levels = names(group_palette)), Date_Parsed = as.Date(Date,
    format = "%d_%m_%Y")) %>%
  arrange(Date_Parsed, Experiment)

shared_samples <- intersect(Species_df$Sample_ID, Metadata_df$Sample_ID)

Species_df <- Species_df %>%
  filter(Sample_ID %in% shared_samples)

Metadata_df <- Metadata_df %>%
  filter(Sample_ID %in% shared_samples)

date_levels <- Metadata_df %>%
  arrange(Date_Parsed, Experiment) %>%
  distinct(Date) %>%
  pull(Date)

top_viruses <- Species_df %>%
  dplyr::select(-Sample_ID) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Virus", values_to = "Total_Abundance") %>%
  arrange(desc(Total_Abundance)) %>%
  slice_head(n = 19) %>%
  pull(Virus)

Species_long <- Species_df %>%
  pivot_longer(cols = -Sample_ID, names_to = "Virus", values_to = "Abundance") %>%
  mutate(Virus = if_else(Virus %in% top_viruses, Virus, "Other")) %>%
  group_by(Sample_ID, Virus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(Relative_Abundance = Abundance/sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(Metadata_df %>%
    dplyr::select(Sample_ID, Experiment, Sample_Type, Date, Date_Parsed), by = "Sample_ID")

virus_levels <- c(top_viruses, "Other")
virus_palette <- c(setNames(tableau20[seq_along(top_viruses)], top_viruses), Other = other_color)

Species_long <- Species_long %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID, Virus) %>%
  mutate(Virus = factor(Virus, levels = virus_levels), Date = factor(Date, levels = date_levels),
    Sample_ID = factor(Sample_ID, levels = unique(Sample_ID)))

plot <- ggplot(Species_long, aes(x = Sample_ID, y = Relative_Abundance, fill = Virus)) + geom_col(width = 1,
  color = NA, linewidth = 0) + facet_grid(~Date, scales = "free_x", space = "free_x") + scale_fill_manual(values = virus_palette) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  expand = c(0, 0)) + labs(x = "Sample", y = "Relative abundance", fill = "Virus") + theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(),
    legend.position = "bottom", panel.spacing.x = grid::unit(0, "pt"))

print(plot)
# ggsave(output_file('Relative_Abund_Stacked_Barplot_Sampling_Date.pdf'), plot = plot, width = 200, height = 150, units = 'mm')


# 23. Relate virus detections to Jaccard centroid distance ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Species_RPM <- read.xlsx(input_file("Pig_Only_Species_RPM.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

# Manually edit these vectors to control which labels appear in each panel.
label_air <- c("Coronavirus_HKU15", "Mamastrovirus_2", "Rotavirus_A", "Bastrovirus_VietNam_Porcine_17489_85",
  "Porcine_bastrovirus", "Porcine_mastadenovirus_C", "Rotavirus_B", "Porcine_picobirnavirus",
  "Porcine_torovirus", "Porcine_astrovirus_2", "Porcine_adenovirus", "Rotavirus_C", "Porcine_astrovirus_4",
  "Pig_astrovirus_ahast_2", "Porcine_epidemic_diarrhea_virus", "Mamastrovirus_3", "Bastrovirus_like_virus_VietNam_Bat_17819_21",
  "Porcine_bocavirus", "Rotavirus_G", "Enterovirus_G", "Astrovirus_wild_boar_WBAstV_1_2011_HUN")

label_slaughter_box <- c("Coronavirus_HKU15", "Mamastrovirus_2", "Rotavirus_A", "Bastrovirus_VietNam_Porcine_17489_85",
  "Porcine_bastrovirus", "Porcine_mastadenovirus_C", "Rotavirus_B", "Porcine_picobirnavirus",
  "Porcine_torovirus", "Porcine_astrovirus_2", "Porcine_adenovirus", "Rotavirus_C", "Porcine_astrovirus_4",
  "Pig_astrovirus_ahast_2", "Porcine_epidemic_diarrhea_virus", "Mamastrovirus_3", "Bastrovirus_like_virus_VietNam_Bat_17819_21",
  "Porcine_bocavirus", "Rotavirus_G", "Enterovirus_G", "Astrovirus_wild_boar_WBAstV_1_2011_HUN",
  "Mammalian_orthoreovirus", "Mamastrovirus_4", "Astrovirus_swine_PoAstV12_4_Canada_2006",
  "Sapporo_virus", "Astrovirus_swine_PoAstV12_Canada_2006", "Pig_astrovirus_CX1", "Porcine_stool_associated_circular_virus_BEL_15V010",
  "Porcine_astrovirus_5", "Mamastrovirus_pig_HL_China")

label_pig <- c("Coronavirus_HKU15", "Mamastrovirus_2", "Rotavirus_A", "Bastrovirus_VietNam_Porcine_17489_85",
  "Porcine_bastrovirus", "Porcine_mastadenovirus_C", "Rotavirus_B", "Porcine_picobirnavirus",
  "Porcine_torovirus", "Porcine_astrovirus_2", "Porcine_adenovirus", "Rotavirus_C", "Porcine_astrovirus_4",
  "Pig_astrovirus_ahast_2", "Porcine_epidemic_diarrhea_virus", "Mamastrovirus_3", "Bastrovirus_like_virus_VietNam_Bat_17819_21",
  "Porcine_bocavirus", "Rotavirus_G", "Enterovirus_G", "Astrovirus_wild_boar_WBAstV_1_2011_HUN",
  "Mamastrovirus_4", "Pig_astrovirus_CX1", "Sapporo_virus", "Rotavirus_D", "Astrovirus_swine_PoAstV12_4_Canada_2006",
  "Porcine_astrovirus_3", "Mamastrovirus_pig_HL_China", "Porcine_astrovirus_5")

label_wash_water <- c("Coronavirus_HKU15", "Mamastrovirus_2", "Rotavirus_A", "Bastrovirus_VietNam_Porcine_17489_85",
  "Porcine_bastrovirus", "Porcine_mastadenovirus_C", "Rotavirus_B", "Porcine_picobirnavirus",
  "Porcine_torovirus", "Porcine_astrovirus_2", "Porcine_adenovirus", "Rotavirus_C", "Porcine_astrovirus_4",
  "Pig_astrovirus_ahast_2", "Porcine_epidemic_diarrhea_virus", "Mamastrovirus_3", "Bastrovirus_like_virus_VietNam_Bat_17819_21",
  "Porcine_bocavirus", "Rotavirus_G", "Enterovirus_G", "Astrovirus_wild_boar_WBAstV_1_2011_HUN",
  "Mamastrovirus_4", "Pig_astrovirus_CX1", "Astrovirus_swine_PoAstV12_4_Canada_2006")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]

OTU <- otu_table(t(Species_bin), taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
comm_mat <- t(as(otu_table(physeq), "matrix"))

Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
pal <- group_palette[present_groups]

Jaccard <- vegdist(comm_mat, method = "jaccard", binary = TRUE)
Jaccard <- as.dist(Jaccard)

Jaccard_PCoA <- dudi.pco(Jaccard, scannf = FALSE, nf = 3)

coords_df <- Jaccard_PCoA$li[, 1:3, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID")

Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID")

coords_df <- coords_df %>%
  left_join(Metadata_df %>%
    dplyr::select(Sample_ID, Sample_Type, Experiment), by = "Sample_ID")

axis_cols <- c("Axis1", "Axis2", "Axis3")
colnames(coords_df)[2:4] <- axis_cols

build_centroid_volcano <- function(target_group, label_viruses = character(0), y_axis_position = "left") {
  if (!target_group %in% coords_df$Sample_Type) {
    stop(paste0("No samples with Sample_Type == '", target_group, "' were found."))
  }

  target_centroid <- coords_df %>%
    filter(Sample_Type == target_group) %>%
    summarise(across(all_of(axis_cols), median, na.rm = TRUE))

  distance_col <- paste0("Distance_From_", target_group, "_Median")
  mean_presence_col <- paste0("Mean_", target_group, "_Presence")

  distance_df <- coords_df %>%
    mutate(!!distance_col := sqrt((Axis1 - target_centroid$Axis1)^2 + (Axis2 - target_centroid$Axis2)^2 +
      (Axis3 - target_centroid$Axis3)^2)) %>%
    dplyr::select(Sample_ID, Sample_Type, Experiment, all_of(distance_col))

  shared_samples <- intersect(rownames(Species_bin), distance_df$Sample_ID)

  Species_sub <- Species_bin[shared_samples, , drop = FALSE]
  Distance_vec <- distance_df %>%
    filter(Sample_ID %in% shared_samples) %>%
    arrange(match(Sample_ID, shared_samples)) %>%
    pull(all_of(distance_col))

  target_samples <- distance_df %>%
    filter(Sample_Type == target_group) %>%
    pull(Sample_ID)

  target_samples <- intersect(rownames(Species_bin), target_samples)

  Mean_Presence <- Species_bin[target_samples, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample_ID") %>%
    pivot_longer(cols = -Sample_ID, names_to = "Virus", values_to = "Presence") %>%
    group_by(Virus) %>%
    summarise(!!mean_presence_col := mean(Presence, na.rm = TRUE), .groups = "drop")

  Spearman_df <- lapply(colnames(Species_sub), function(virus_name) {
    abundance_vec <- Species_sub[, virus_name]
    keep <- is.finite(Distance_vec) & is.finite(abundance_vec)

    if (sum(keep) < 3 || length(unique(abundance_vec[keep])) < 2) {
      return(data.frame(Virus = virus_name, Spearman_Rho = NA_real_, P_Value = NA_real_))
    }

    test_out <- suppressWarnings(cor.test(x = Distance_vec[keep], y = abundance_vec[keep],
      method = "spearman", exact = FALSE))

    data.frame(Virus = virus_name, Spearman_Rho = unname(test_out$estimate), P_Value = test_out$p.value)
  }) %>%
    bind_rows() %>%
    left_join(Mean_Presence, by = "Virus") %>%
    mutate(FDR_P_Value = p.adjust(P_Value, method = "fdr"), NegLog10_FDR = -log10(FDR_P_Value),
      Significant = if_else(!is.na(FDR_P_Value) & FDR_P_Value < 0.05, "FDR < 0.05", "NS")) %>%
    arrange(FDR_P_Value, desc(abs(Spearman_Rho)))

  significant_taxa <- Spearman_df %>%
    filter(Significant == "FDR < 0.05")

  volcano_labels <- Spearman_df %>%
    filter(Significant == "FDR < 0.05", Virus %in% label_viruses)

  point_jitter_width <- 0.1
  point_jitter_height <- 0.1

  plot_obj <- ggplot(Spearman_df, aes(x = Spearman_Rho, y = NegLog10_FDR)) + geom_point(aes(size = .data[[mean_presence_col]],
    color = Significant), alpha = 0.8, position = position_jitter(width = point_jitter_width,
    height = point_jitter_height)) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") + geom_text_repel(data = volcano_labels,
    aes(label = Virus), size = 3.75, color = "black", max.overlaps = 50) + scale_color_manual(values = c(`FDR < 0.05` = unname(group_palette[target_group]),
    NS = "grey70"), breaks = c("FDR < 0.05", "NS")) + scale_size_continuous(range = c(2.5,
    12.5)) + scale_x_continuous() + scale_y_continuous(position = y_axis_position) + labs(title = target_group,
    x = "Spearman correlation coefficient", y = "-log10(FDR p-value)", size = paste0("Mean presence in ",
      target_group)) + theme_classic(base_size = 15) + theme(legend.position = "bottom",
    legend.box = "horizontal", plot.title = element_text(face = "bold"))

  list(results = Spearman_df, significant_taxa = significant_taxa, labels = volcano_labels,
    plot = plot_obj)
}

Air_Output <- build_centroid_volcano(target_group = "Air", label_viruses = label_air, y_axis_position = "left")

Slaughter_Box_Output <- build_centroid_volcano(target_group = "Slaughter_Box", label_viruses = label_slaughter_box,
  y_axis_position = "right")

Pig_Output <- build_centroid_volcano(target_group = "Pig", label_viruses = label_pig, y_axis_position = "left")

Wash_Water_Output <- build_centroid_volcano(target_group = "Wash_Water", label_viruses = label_wash_water,
  y_axis_position = "right")

Air_Significant_Taxa <- Air_Output$significant_taxa
Slaughter_Box_Significant_Taxa <- Slaughter_Box_Output$significant_taxa
Pig_Significant_Taxa <- Pig_Output$significant_taxa
Wash_Water_Significant_Taxa <- Wash_Water_Output$significant_taxa

plot_air <- Air_Output$plot
plot_slaughter_box <- Slaughter_Box_Output$plot
plot_pig <- Pig_Output$plot
plot_wash_water <- Wash_Water_Output$plot

Air_Slaughter_Box_Centroid_Volcano_Plot <- (plot_air + plot_slaughter_box) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

Pig_Wash_Water_Centroid_Volcano_Plot <- (plot_pig + plot_wash_water) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(Air_Slaughter_Box_Centroid_Volcano_Plot)
print(Pig_Wash_Water_Centroid_Volcano_Plot)

# ggsave(output_file('Distance_Air_SlaughterBox_Spearman_Binary_Volcano.pdf'), Air_Slaughter_Box_Centroid_Volcano_Plot, width = 320, height = 200, units = 'mm')
# ggsave(output_file('Distance_Pig_WashWater_Spearman_Binary_Volcano.pdf'), Pig_Wash_Water_Centroid_Volcano_Plot, width = 320, height = 200, units = 'mm')


# 24. Summarise viral contributions to Jaccard similarity ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

shared_samples <- intersect(rownames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[shared_samples, , drop = FALSE]
Metadata <- Metadata[shared_samples, , drop = FALSE]

Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = names(group_palette))

present_groups <- intersect(names(group_palette), unique(as.character(Metadata$Sample_Type)))
comparison_groups <- setdiff(present_groups, "Pig")

Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID")

build_pairings <- function(metadata_df, comparison_group) {
  strata_df <- metadata_df %>%
    filter(Sample_Type %in% c("Pig", comparison_group)) %>%
    group_by(Experiment, Location, Sample_Type) %>%
    summarise(Sample_IDs = list(Sample_ID), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample_Type, values_from = Sample_IDs) %>%
    filter(!is.na(Pig), !is.na(.data[[comparison_group]]))

  lapply(seq_len(nrow(strata_df)), function(i) {
    list(pig_ids = unlist(strata_df$Pig[[i]]), other_ids = unlist(strata_df[[comparison_group]][[i]]))
  })
}

jaccard_similarity_pairs <- function(mat_a, mat_b) {
  if (nrow(mat_a) == 0 || nrow(mat_b) == 0 || ncol(mat_a) == 0) {
    return(numeric(0))
  }

  mat_a <- as.matrix(mat_a)
  mat_b <- as.matrix(mat_b)

  intersection <- tcrossprod(mat_a, mat_b)
  union <- outer(rowSums(mat_a), rowSums(mat_b), "+") - intersection
  similarity <- intersection/union
  similarity[union == 0] <- NA_real_

  as.vector(similarity)
}

mean_group_similarity <- function(species_matrix, pairings) {
  similarity_values <- unlist(lapply(pairings, function(pair_info) {
    jaccard_similarity_pairs(species_matrix[pair_info$pig_ids, , drop = FALSE], species_matrix[pair_info$other_ids,
      , drop = FALSE])
  }), use.names = FALSE)

  similarity_values <- similarity_values[is.finite(similarity_values)]

  if (length(similarity_values) == 0) {
    return(c(mean_similarity = NA_real_, n_pairs = 0))
  }

  c(mean_similarity = mean(similarity_values, na.rm = TRUE), n_pairs = length(similarity_values))
}

comparison_pairings <- lapply(comparison_groups, function(comparison_group) {
  build_pairings(Metadata_df, comparison_group)
})
names(comparison_pairings) <- comparison_groups

Baseline_Jaccard <- lapply(comparison_groups, function(comparison_group) {
  summary_stats <- mean_group_similarity(Species_bin, comparison_pairings[[comparison_group]])

  tibble(Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = unname(summary_stats["mean_similarity"]),
    N_Pairs = unname(summary_stats["n_pairs"]))
}) %>%
  bind_rows()

virus_order <- colSums(Species_bin, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

Virus_Jaccard_Contribution_Long <- lapply(comparison_groups, function(comparison_group) {
  baseline_similarity <- Baseline_Jaccard %>%
    filter(Comparison == paste0("Pig_vs_", comparison_group)) %>%
    pull(Baseline_Similarity)

  lapply(virus_order, function(virus_name) {
    reduced_matrix <- Species_bin[, colnames(Species_bin) != virus_name, drop = FALSE]

    reduced_stats <- mean_group_similarity(reduced_matrix, comparison_pairings[[comparison_group]])

    similarity_without_virus <- unname(reduced_stats["mean_similarity"])
    contribution_value <- baseline_similarity - similarity_without_virus

    tibble(Virus = virus_name, Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = baseline_similarity,
      Similarity_Without_Virus = similarity_without_virus, Contribution = contribution_value,
      Interpretation = case_when(is.na(contribution_value) ~ NA_character_, contribution_value >
        0 ~ "Contributes to similarity", contribution_value < 0 ~ "Contributes to dissimilarity",
        TRUE ~ "Neutral"))
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Wide <- Virus_Jaccard_Contribution_Long %>%
  mutate(Comparison = factor(Comparison, levels = paste0("Pig_vs_", comparison_groups))) %>%
  arrange(match(Virus, virus_order), Comparison) %>%
  dplyr::select(Virus, Comparison, Contribution) %>%
  pivot_wider(names_from = Comparison, values_from = Contribution)

Virus_Jaccard_Contribution_Mean <- Virus_Jaccard_Contribution_Long %>%
  group_by(Virus) %>%
  summarise(Mean_Contribution = mean(Contribution, na.rm = TRUE), Mean_Abs_Contribution = mean(abs(Contribution),
    na.rm = TRUE), Max_Contribution = max(Contribution, na.rm = TRUE), Min_Contribution = min(Contribution,
    na.rm = TRUE), N_Comparisons = sum(!is.na(Contribution)), .groups = "drop") %>%
  arrange(desc(Mean_Contribution))

Contribution_Matrix <- Virus_Jaccard_Contribution_Wide %>%
  as.data.frame()

rownames(Contribution_Matrix) <- Contribution_Matrix$Virus
Contribution_Matrix$Virus <- NULL
Contribution_Matrix <- as.matrix(Contribution_Matrix)

Contribution_Matrix <- Contribution_Matrix[virus_order[virus_order %in% rownames(Contribution_Matrix)],
  , drop = FALSE]

Contribution_Matrix_For_PCA <- Contribution_Matrix
Contribution_Matrix_For_PCA[is.na(Contribution_Matrix_For_PCA)] <- 0

if (nrow(Contribution_Matrix_For_PCA) >= 2 && ncol(Contribution_Matrix_For_PCA) >= 2) {
  Virus_Jaccard_PCA <- prcomp(Contribution_Matrix_For_PCA, center = TRUE, scale. = TRUE)

  Virus_Jaccard_PCA_Coords <- as.data.frame(Virus_Jaccard_PCA$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column(var = "Virus") %>%
    mutate(Max_Abs_Contribution = apply(Contribution_Matrix, 1, function(x) max(abs(x),
      na.rm = TRUE)), Dominant_Comparison = apply(abs(Contribution_Matrix_For_PCA), 1,
      function(x) colnames(Contribution_Matrix_For_PCA)[which.max(x)]), Dominant_Group = sub("^Pig_vs_",
      "", Dominant_Comparison)) %>%
    left_join(Virus_Jaccard_Contribution_Mean, by = "Virus")
} else {
  Virus_Jaccard_PCA_Coords <- tibble(Virus = rownames(Contribution_Matrix_For_PCA), PC1 = 0,
    PC2 = 0, Max_Abs_Contribution = apply(Contribution_Matrix, 1, function(x) max(abs(x),
      na.rm = TRUE)), Dominant_Comparison = NA_character_, Dominant_Group = NA_character_) %>%
    left_join(Virus_Jaccard_Contribution_Mean, by = "Virus")
}

overall_label_viruses <- Virus_Jaccard_PCA_Coords %>%
  filter(!is.na(Mean_Contribution)) %>%
  arrange(desc(Mean_Contribution), desc(Mean_Abs_Contribution)) %>%
  slice_head(n = 20)

Virus_Jaccard_Overall_Similarity_Plot <- ggplot(Virus_Jaccard_PCA_Coords, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = Mean_Abs_Contribution, fill = Mean_Contribution), shape = 21, color = "black",
    stroke = 0.3, alpha = 0.85) + geom_text_repel(data = overall_label_viruses, aes(label = Virus),
  color = "black", size = 3, max.overlaps = 25) + scale_fill_gradient2(low = "#3B75AF", mid = "white",
  high = "#CC3311", midpoint = 0, na.value = "grey80") + scale_size_continuous(range = c(2,
  10)) + labs(x = "PC1", y = "PC2", fill = "Mean contribution", size = "Mean abs contribution") +
  theme_classic(base_size = 12)
print(Virus_Jaccard_Overall_Similarity_Plot)
# ggsave(output_file('Jaccard_Similairty_Contribution_Virus_PCA.pdf'), plot = Virus_Jaccard_Overall_Similarity_Plot, width = 200, height = 150, units = 'mm')


# 25. Summarise viral contributions to Jaccard dissimilarity ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

shared_samples <- intersect(rownames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[shared_samples, , drop = FALSE]
Metadata <- Metadata[shared_samples, , drop = FALSE]

Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = names(group_palette))
Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID")

present_groups <- intersect(names(group_palette), unique(as.character(Metadata$Sample_Type)))
comparison_groups <- setdiff(present_groups, "Pig")

build_pairings <- function(metadata_df, comparison_group) {
  strata_df <- metadata_df %>%
    filter(Sample_Type %in% c("Pig", comparison_group)) %>%
    group_by(Experiment, Location, Sample_Type) %>%
    summarise(Sample_IDs = list(Sample_ID), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample_Type, values_from = Sample_IDs) %>%
    filter(!is.na(Pig), !is.na(.data[[comparison_group]]))

  lapply(seq_len(nrow(strata_df)), function(i) {
    list(pig_ids = unlist(strata_df$Pig[[i]]), other_ids = unlist(strata_df[[comparison_group]][[i]]))
  })
}

jaccard_similarity_pairs <- function(mat_a, mat_b) {
  if (nrow(mat_a) == 0 || nrow(mat_b) == 0 || ncol(mat_a) == 0) {
    return(numeric(0))
  }

  mat_a <- as.matrix(mat_a)
  mat_b <- as.matrix(mat_b)

  intersection <- tcrossprod(mat_a, mat_b)
  union <- outer(rowSums(mat_a), rowSums(mat_b), "+") - intersection
  similarity <- intersection/union
  similarity[union == 0] <- NA_real_

  as.vector(similarity)
}

mean_group_similarity <- function(species_matrix, pairings) {
  similarity_values <- unlist(lapply(pairings, function(pair_info) {
    jaccard_similarity_pairs(species_matrix[pair_info$pig_ids, , drop = FALSE], species_matrix[pair_info$other_ids,
      , drop = FALSE])
  }), use.names = FALSE)

  similarity_values <- similarity_values[is.finite(similarity_values)]

  if (length(similarity_values) == 0) {
    return(c(mean_similarity = NA_real_, n_pairs = 0))
  }

  c(mean_similarity = mean(similarity_values, na.rm = TRUE), n_pairs = length(similarity_values))
}

virus_order <- colSums(Species_bin, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

comparison_pairings <- lapply(comparison_groups, function(comparison_group) {
  build_pairings(Metadata_df, comparison_group)
})
names(comparison_pairings) <- comparison_groups

Baseline_Jaccard <- lapply(comparison_groups, function(comparison_group) {
  summary_stats <- mean_group_similarity(Species_bin, comparison_pairings[[comparison_group]])

  tibble(Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = unname(summary_stats["mean_similarity"]),
    N_Pairs = unname(summary_stats["n_pairs"]))
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Long <- lapply(comparison_groups, function(comparison_group) {
  baseline_similarity <- Baseline_Jaccard %>%
    filter(Comparison == paste0("Pig_vs_", comparison_group)) %>%
    pull(Baseline_Similarity)

  lapply(virus_order, function(virus_name) {
    reduced_matrix <- Species_bin[, colnames(Species_bin) != virus_name, drop = FALSE]

    similarity_without_virus <- mean_group_similarity(reduced_matrix, comparison_pairings[[comparison_group]])["mean_similarity"] %>%
      unname()

    tibble(Virus = virus_name, Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = baseline_similarity,
      Similarity_Without_Virus = similarity_without_virus, Contribution = baseline_similarity -
        similarity_without_virus)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Mean <- Virus_Jaccard_Contribution_Long %>%
  group_by(Virus) %>%
  summarise(Mean_Contribution = mean(Contribution, na.rm = TRUE), Mean_Abs_Contribution = mean(abs(Contribution),
    na.rm = TRUE), N_Comparisons = sum(!is.na(Contribution)), .groups = "drop")

Virus_Jaccard_Contribution_Wide <- Virus_Jaccard_Contribution_Long %>%
  mutate(Comparison = factor(Comparison, levels = paste0("Pig_vs_", comparison_groups))) %>%
  arrange(match(Virus, virus_order), Comparison) %>%
  dplyr::select(Virus, Comparison, Contribution) %>%
  pivot_wider(names_from = Comparison, values_from = Contribution)

Contribution_Matrix <- Virus_Jaccard_Contribution_Wide %>%
  as.data.frame()

rownames(Contribution_Matrix) <- Contribution_Matrix$Virus
Contribution_Matrix$Virus <- NULL
Contribution_Matrix <- as.matrix(Contribution_Matrix)

Contribution_Matrix <- Contribution_Matrix[virus_order[virus_order %in% rownames(Contribution_Matrix)],
  , drop = FALSE]

Contribution_Matrix_For_PCA <- Contribution_Matrix
Contribution_Matrix_For_PCA[is.na(Contribution_Matrix_For_PCA)] <- 0

if (nrow(Contribution_Matrix_For_PCA) >= 2 && ncol(Contribution_Matrix_For_PCA) >= 2) {
  Virus_Jaccard_PCA <- prcomp(Contribution_Matrix_For_PCA, center = TRUE, scale. = TRUE)

  Virus_Jaccard_PCA_Coords <- as.data.frame(Virus_Jaccard_PCA$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column(var = "Virus") %>%
    setNames(c("Virus", "PC1", "PC2")) %>%
    left_join(Virus_Jaccard_Contribution_Mean, by = "Virus")
} else {
  Virus_Jaccard_PCA_Coords <- tibble(Virus = rownames(Contribution_Matrix_For_PCA), PC1 = 0,
    PC2 = 0) %>%
    left_join(Virus_Jaccard_Contribution_Mean, by = "Virus")
}

Overall_Dissimilarity_Bubble <- Virus_Jaccard_PCA_Coords %>%
  dplyr::filter(!is.na(Mean_Contribution), Mean_Contribution < 0)

Overall_Dissimilarity_Labels <- Overall_Dissimilarity_Bubble %>%
  dplyr::arrange(Mean_Contribution) %>%
  dplyr::slice_head(n = 20)

plot_overall_dissimilarity_bubble <- ggplot(Overall_Dissimilarity_Bubble, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = Mean_Abs_Contribution, fill = Mean_Contribution), shape = 21, color = "black",
    stroke = 0.3, alpha = 0.85) + geom_text_repel(data = Overall_Dissimilarity_Labels,
  aes(label = Virus), size = 3.2, max.overlaps = 25) + scale_fill_gradient2(low = "#2166AC",
  mid = "#FEE5D9", high = "#FEE5D9", midpoint = 0, limits = c(min(Overall_Dissimilarity_Bubble$Mean_Contribution,
    na.rm = TRUE), 0), na.value = "grey80") + labs(x = "PC1", y = "PC2", fill = "Mean Contribution",
  size = "Mean absolute effect", title = "Top 20 Overall Drivers of Jaccard Dissimilarity") +
  theme_classic(base_size = 12)

print(plot_overall_dissimilarity_bubble)
# ggsave(output_file('Jaccard_DISSIMILARITY_Contribution_Virus_PCA.pdf'), plot = plot_overall_dissimilarity_bubble, width = 200, height = 150, units = 'mm')


# 26. Plot viral contributions to Jaccard dissimilarity ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

shared_samples <- intersect(rownames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[shared_samples, , drop = FALSE]
Metadata <- Metadata[shared_samples, , drop = FALSE]

Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = names(group_palette))
Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID")

present_groups <- intersect(names(group_palette), unique(as.character(Metadata$Sample_Type)))
comparison_groups <- setdiff(present_groups, "Pig")

build_pairings <- function(metadata_df, comparison_group) {
  strata_df <- metadata_df %>%
    filter(Sample_Type %in% c("Pig", comparison_group)) %>%
    group_by(Experiment, Location, Sample_Type) %>%
    summarise(Sample_IDs = list(Sample_ID), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample_Type, values_from = Sample_IDs) %>%
    filter(!is.na(Pig), !is.na(.data[[comparison_group]]))

  lapply(seq_len(nrow(strata_df)), function(i) {
    list(pig_ids = unlist(strata_df$Pig[[i]]), other_ids = unlist(strata_df[[comparison_group]][[i]]))
  })
}

jaccard_similarity_pairs <- function(mat_a, mat_b) {
  if (nrow(mat_a) == 0 || nrow(mat_b) == 0 || ncol(mat_a) == 0) {
    return(numeric(0))
  }

  mat_a <- as.matrix(mat_a)
  mat_b <- as.matrix(mat_b)

  intersection <- tcrossprod(mat_a, mat_b)
  union <- outer(rowSums(mat_a), rowSums(mat_b), "+") - intersection
  similarity <- intersection/union
  similarity[union == 0] <- NA_real_

  as.vector(similarity)
}

mean_group_similarity <- function(species_matrix, pairings) {
  similarity_values <- unlist(lapply(pairings, function(pair_info) {
    jaccard_similarity_pairs(species_matrix[pair_info$pig_ids, , drop = FALSE], species_matrix[pair_info$other_ids,
      , drop = FALSE])
  }), use.names = FALSE)

  similarity_values <- similarity_values[is.finite(similarity_values)]

  if (length(similarity_values) == 0) {
    return(c(mean_similarity = NA_real_, n_pairs = 0))
  }

  c(mean_similarity = mean(similarity_values, na.rm = TRUE), n_pairs = length(similarity_values))
}

virus_order <- colSums(Species_bin, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

comparison_pairings <- lapply(comparison_groups, function(comparison_group) {
  build_pairings(Metadata_df, comparison_group)
})
names(comparison_pairings) <- comparison_groups

Baseline_Jaccard <- lapply(comparison_groups, function(comparison_group) {
  summary_stats <- mean_group_similarity(Species_bin, comparison_pairings[[comparison_group]])

  tibble(Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = unname(summary_stats["mean_similarity"]),
    N_Pairs = unname(summary_stats["n_pairs"]))
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Long <- lapply(comparison_groups, function(comparison_group) {
  baseline_similarity <- Baseline_Jaccard %>%
    filter(Comparison == paste0("Pig_vs_", comparison_group)) %>%
    pull(Baseline_Similarity)

  lapply(virus_order, function(virus_name) {
    reduced_matrix <- Species_bin[, colnames(Species_bin) != virus_name, drop = FALSE]

    similarity_without_virus <- mean_group_similarity(reduced_matrix, comparison_pairings[[comparison_group]])["mean_similarity"] %>%
      unname()

    tibble(Virus = virus_name, Comparison = paste0("Pig_vs_", comparison_group), Contribution = baseline_similarity -
      similarity_without_virus)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Mean <- Virus_Jaccard_Contribution_Long %>%
  group_by(Virus) %>%
  summarise(Mean_Contribution = mean(Contribution, na.rm = TRUE), Mean_Abs_Contribution = mean(abs(Contribution),
    na.rm = TRUE), .groups = "drop")

build_lollipop_plot <- function(df, panel_title, panel_color, contribution_col = "Contribution") {
  plot_df <- df %>%
    filter(!is.na(.data[[contribution_col]]), .data[[contribution_col]] < 0) %>%
    arrange(.data[[contribution_col]]) %>%
    slice_head(n = 10) %>%
    mutate(Virus = factor(Virus, levels = rev(Virus)))

  ggplot(plot_df, aes(x = .data[[contribution_col]], y = Virus)) + geom_vline(xintercept = 0,
    linetype = "dashed", color = "grey60") + geom_segment(aes(x = 0, xend = .data[[contribution_col]],
    yend = Virus), color = "grey75", linewidth = 0.7) + geom_point(color = panel_color,
    size = 3.2, alpha = 0.95) + labs(title = panel_title, x = "Dissimilarity effect", y = NULL) +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 8))
}

panel_order <- c("Cage", "Fecal", "Slaughter_Box", "Air", "Rope", "Waste_Water", "Wash_Water")
panel_order <- panel_order[panel_order %in% comparison_groups]

panel_plots <- lapply(panel_order, function(group_name) {
  comparison_name <- paste0("Pig_vs_", group_name)

  plot_df <- Virus_Jaccard_Contribution_Long %>%
    filter(Comparison == comparison_name)

  build_lollipop_plot(df = plot_df, panel_title = group_name, panel_color = unname(group_palette[group_name]),
    contribution_col = "Contribution")
})

overall_plot <- build_lollipop_plot(df = Virus_Jaccard_Contribution_Mean %>%
  dplyr::transmute(Virus = Virus, Contribution = Mean_Contribution, Mean_Abs_Contribution = Mean_Abs_Contribution),
  panel_title = "Overall", panel_color = "#4D4D4D", contribution_col = "Contribution")

all_plots <- c(panel_plots, list(overall_plot))
Jaccard_Dissimilarity_Lollipop_Grid <- wrap_plots(all_plots, ncol = 2)
print(Jaccard_Dissimilarity_Lollipop_Grid)
ggsave(output_file("Jaccard_Dissimilarity_Lollipop_Grid.pdf"), Jaccard_Dissimilarity_Lollipop_Grid,
  width = 250, height = 300, units = "mm")


# 27. Plot viral contributions to Jaccard similarity ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]

shared_samples <- intersect(rownames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[shared_samples, , drop = FALSE]
Metadata <- Metadata[shared_samples, , drop = FALSE]

Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = names(group_palette))
Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID")

present_groups <- intersect(names(group_palette), unique(as.character(Metadata$Sample_Type)))
comparison_groups <- setdiff(present_groups, "Pig")

build_pairings <- function(metadata_df, comparison_group) {
  strata_df <- metadata_df %>%
    filter(Sample_Type %in% c("Pig", comparison_group)) %>%
    group_by(Experiment, Location, Sample_Type) %>%
    summarise(Sample_IDs = list(Sample_ID), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample_Type, values_from = Sample_IDs) %>%
    filter(!is.na(Pig), !is.na(.data[[comparison_group]]))

  lapply(seq_len(nrow(strata_df)), function(i) {
    list(pig_ids = unlist(strata_df$Pig[[i]]), other_ids = unlist(strata_df[[comparison_group]][[i]]))
  })
}

jaccard_similarity_pairs <- function(mat_a, mat_b) {
  if (nrow(mat_a) == 0 || nrow(mat_b) == 0 || ncol(mat_a) == 0) {
    return(numeric(0))
  }

  mat_a <- as.matrix(mat_a)
  mat_b <- as.matrix(mat_b)

  intersection <- tcrossprod(mat_a, mat_b)
  union <- outer(rowSums(mat_a), rowSums(mat_b), "+") - intersection
  similarity <- intersection/union
  similarity[union == 0] <- NA_real_

  as.vector(similarity)
}

mean_group_similarity <- function(species_matrix, pairings) {
  similarity_values <- unlist(lapply(pairings, function(pair_info) {
    jaccard_similarity_pairs(species_matrix[pair_info$pig_ids, , drop = FALSE], species_matrix[pair_info$other_ids,
      , drop = FALSE])
  }), use.names = FALSE)

  similarity_values <- similarity_values[is.finite(similarity_values)]

  if (length(similarity_values) == 0) {
    return(c(mean_similarity = NA_real_, n_pairs = 0))
  }

  c(mean_similarity = mean(similarity_values, na.rm = TRUE), n_pairs = length(similarity_values))
}

virus_order <- colSums(Species_bin, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

comparison_pairings <- lapply(comparison_groups, function(comparison_group) {
  build_pairings(Metadata_df, comparison_group)
})
names(comparison_pairings) <- comparison_groups

Baseline_Jaccard <- lapply(comparison_groups, function(comparison_group) {
  summary_stats <- mean_group_similarity(Species_bin, comparison_pairings[[comparison_group]])

  tibble(Comparison = paste0("Pig_vs_", comparison_group), Baseline_Similarity = unname(summary_stats["mean_similarity"]),
    N_Pairs = unname(summary_stats["n_pairs"]))
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Long <- lapply(comparison_groups, function(comparison_group) {
  baseline_similarity <- Baseline_Jaccard %>%
    filter(Comparison == paste0("Pig_vs_", comparison_group)) %>%
    pull(Baseline_Similarity)

  lapply(virus_order, function(virus_name) {
    reduced_matrix <- Species_bin[, colnames(Species_bin) != virus_name, drop = FALSE]

    similarity_without_virus <- mean_group_similarity(reduced_matrix, comparison_pairings[[comparison_group]])["mean_similarity"] %>%
      unname()

    tibble(Virus = virus_name, Comparison = paste0("Pig_vs_", comparison_group), Contribution = baseline_similarity -
      similarity_without_virus)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

Virus_Jaccard_Contribution_Mean <- Virus_Jaccard_Contribution_Long %>%
  group_by(Virus) %>%
  summarise(Mean_Contribution = mean(Contribution, na.rm = TRUE), Mean_Abs_Contribution = mean(abs(Contribution),
    na.rm = TRUE), .groups = "drop")

build_lollipop_plot <- function(df, panel_title, panel_color, contribution_col = "Contribution") {
  plot_df <- df %>%
    filter(!is.na(.data[[contribution_col]]), .data[[contribution_col]] > 0) %>%
    arrange(desc(.data[[contribution_col]])) %>%
    slice_head(n = 10) %>%
    mutate(Virus = factor(Virus, levels = rev(Virus)))

  ggplot(plot_df, aes(x = .data[[contribution_col]], y = Virus)) + geom_vline(xintercept = 0,
    linetype = "dashed", color = "grey60") + geom_segment(aes(x = 0, xend = .data[[contribution_col]],
    yend = Virus), color = "grey75", linewidth = 0.7) + geom_point(color = panel_color,
    size = 3.2, alpha = 0.95) + labs(title = panel_title, x = "Similarity effect", y = NULL) +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 8))
}

panel_order <- c("Cage", "Fecal", "Slaughter_Box", "Air", "Rope", "Waste_Water", "Wash_Water")
panel_order <- panel_order[panel_order %in% comparison_groups]

panel_plots <- lapply(panel_order, function(group_name) {
  comparison_name <- paste0("Pig_vs_", group_name)

  plot_df <- Virus_Jaccard_Contribution_Long %>%
    filter(Comparison == comparison_name)

  build_lollipop_plot(df = plot_df, panel_title = group_name, panel_color = unname(group_palette[group_name]),
    contribution_col = "Contribution")
})

overall_plot <- build_lollipop_plot(df = Virus_Jaccard_Contribution_Mean %>%
  dplyr::transmute(Virus = Virus, Contribution = Mean_Contribution, Mean_Abs_Contribution = Mean_Abs_Contribution),
  panel_title = "Overall", panel_color = "#4D4D4D", contribution_col = "Contribution")

all_plots <- c(panel_plots, list(overall_plot))
Jaccard_Similarity_Lollipop_Grid <- wrap_plots(all_plots, ncol = 2)
print(Jaccard_Similarity_Lollipop_Grid)
ggsave(output_file("Jaccard_Similarity_Lollipop_Grid.pdf"), Jaccard_Similarity_Lollipop_Grid,
  width = 250, height = 300, units = "mm")


# 28. Plot observed richness by sampling date ----
Species <- read.xlsx(input_file("Pig_Species_Count.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

Species_bin <- (Species > 0) * 1
Species_bin <- t(Species_bin)
Species_bin <- as.data.frame(Species_bin)

common_ids <- intersect(colnames(Species_bin), rownames(Metadata))
Species_bin <- Species_bin[, common_ids, drop = FALSE]
Metadata <- Metadata[common_ids, , drop = FALSE]

Observed_df <- tibble(Sample_ID = colnames(Species_bin), Observed = colSums(Species_bin, na.rm = TRUE))

Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_Type = factor(Sample_Type, levels = names(group_palette)), Date_Parsed = as.Date(Date,
    format = "%d_%m_%Y")) %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID)

date_levels <- Metadata_df %>%
  arrange(Date_Parsed) %>%
  distinct(Date) %>%
  pull(Date)

present_groups <- names(group_palette)[names(group_palette) %in% unique(as.character(Metadata_df$Sample_Type))]
plot_palette <- group_palette[present_groups]

Observed_plot_df <- Observed_df %>%
  left_join(Metadata_df %>%
    dplyr::select(Sample_ID, Sample_Type, Date, Date_Parsed), by = "Sample_ID") %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID) %>%
  mutate(Date = factor(Date, levels = date_levels), Sample_ID = factor(Sample_ID, levels = unique(Sample_ID)))

plot <- ggplot(Observed_plot_df, aes(x = Sample_ID, y = Observed, fill = Sample_Type)) + geom_col(width = 1,
  color = NA, linewidth = 0) + facet_grid(~Date, scales = "free_x", space = "free_x") + scale_fill_manual(values = plot_palette) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x = "Sample",
  y = "Observed species", fill = "Sample Type") + theme_classic(base_size = 12) + theme(axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), strip.background = element_blank(), legend.position = "bottom",
  panel.spacing.x = grid::unit(4, "mm"))

print(plot)
# ggsave(output_file('Observed_Species_Barplot_Sampling_Date.pdf'), plot = plot, width = 200, height = 150, units = 'mm')


# 29. Plot environmental detections of taxa absent from pig samples ----

ES <- read.xlsx(input_file("Pig_Species_RPM.xlsx"), rowNames = TRUE)
Pig <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

to_drop <- c(
  "Bacillus_virus_SIOphi",
  "Cucurbit_aphid_borne_yellows_virus",
  "Angelonia_flower_break_virus",
  "Pepper_vein_yellows_virus_5",
  "Tobacco_vein_distorting_virus",
  "Pepper_vein_yellows_virus",
  "Pepper_vein_yellows_virus_2",
  "Pepper_yellows_virus",
  "Kennedya_yellow_mosaic_virus",
  "Calopogonium_yellow_vein_tymovirus",
  "Passion_fruit_yellow_mosaic_virus",
  "Asclepias_asymptomatic_virus",
  "Blattella_germanica_densovirus_like_virus",
  "Aphis_glycines_virus_3",
  "Aphid_lethal_paralysis_virus",
  "Sacbrood_virus",
  "Wuhan_insect_virus_11",
  "Formica_exsecta_virus_1",
  "Aphis_gossypii_virus",
  "Cimodo_virus",
  "Chimpanzee_picobirnavirus",
  "Fur_seal_faeces_associated_circular_DNA_virus",
  "Wuhan_insect_virus_33",
  "Rhopalosiphum_padi_virus",
  "Drosophila_C_virus",
  "Spodoptera_exigua_iflavirus_2",
  "Warroolaba_Creek_virus_2"
)

group_palette <- c(
  "Pig" = "#EE3377",
  "Cage" = "#0077BB",
  "Fecal" = "#EE7733",
  "Slaughter_Box" = "#33BBEE",
  "Air" = "#009988",
  "Rope" = "#BBBBBB",
  "Waste_Water" = "#CC3311",
  "Wash_Water" = "#228833",
  "Multiple" = "#6A00A8",
  "None" = "darkgrey"
)

# Manually control which viruses appear in the upper block of the plot.
upper_block_keywords <- c("swine", "pig", "porcine")

upper_block_add <- c(
  "Getah_virus",
  "Torque_teno_sus_virus_1a",
  "Sapovirus_sp_"
)

upper_block_remove <- c(
  # "Example_virus_name"
)

cluster_virus_block <- function(mat, virus_names) {
  if (length(virus_names) == 0) {
    return(character(0))
  }

  block_mat <- mat[virus_names, , drop = FALSE]

  if (nrow(block_mat) == 1) {
    return(rownames(block_mat))
  }
  rownames(block_mat)[hclust(dist(block_mat))$order]
}

parse_one_date <- function(value) {
  date_formats <- c(
    "%d_%m_%Y",
    "%Y_%m_%d",
    "%Y-%m-%d",
    "%Y/%m/%d",
    "%d/%m/%Y",
    "%m/%d/%Y",
    "%d-%m-%Y",
    "%m-%d-%Y",
    "%d.%m.%Y",
    "%Y.%m.%d",
    "%d %b %Y",
    "%d %B %Y",
    "%b %d %Y",
    "%B %d %Y"
  )

  for (fmt in date_formats) {
    parsed <- suppressWarnings(as.Date(value, format = fmt))
    if (!is.na(parsed)) {
      return(parsed)
    }
  }

  as.Date(NA)
}

resolve_date_levels <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[is.na(x_chr) | x_chr == ""] <- "Unknown"
  non_unknown <- unique(x_chr[x_chr != "Unknown"])

  numeric_vals <- suppressWarnings(as.numeric(non_unknown))
  if (length(non_unknown) > 0 && all(!is.na(numeric_vals))) {
    ordered_vals <- as.character(sort(unique(numeric_vals)))
    return(c(ordered_vals, if ("Unknown" %in% x_chr) "Unknown"))
  }

  if (length(non_unknown) > 0) {
    parsed_dates <- as.Date(
      vapply(non_unknown, function(value) {
        parsed <- parse_one_date(value)
        if (is.na(parsed)) {
          return(NA_character_)
        }
        as.character(parsed)
      }, character(1))
    )

    if (all(!is.na(parsed_dates))) {
      ordered_vals <- non_unknown[order(parsed_dates)]
      ordered_vals <- unique(ordered_vals)
      return(c(ordered_vals, if ("Unknown" %in% x_chr) "Unknown"))
    }
  }

  c(unique(x_chr[x_chr != "Unknown"]), if ("Unknown" %in% x_chr) "Unknown")
}

if (!"Sample_Type" %in% colnames(Metadata) && "Sample.Type" %in% colnames(Metadata)) {
  Metadata$Sample_Type <- Metadata$Sample.Type
}

if (!"Sample_Type" %in% colnames(Metadata)) {
  stop("Metadata does not contain 'Sample_Type' or 'Sample.Type'.")
}

date_col_candidates <- c("Date", "Sampling_Date", "Collection_Date", "Experiment")
date_col <- date_col_candidates[date_col_candidates %in% colnames(Metadata)][1]

if (is.na(date_col)) {
  stop("Metadata does not contain a date-like column. Tried: Date, Sampling_Date, Collection_Date, Experiment.")
}

ES <- ES[, !colnames(ES) %in% to_drop, drop = FALSE]

cols_to_drop <- intersect(colnames(ES), colnames(Pig))
ES_filtered <- ES[, !colnames(ES) %in% cols_to_drop, drop = FALSE]

non_pig_ids <- rownames(Metadata)[Metadata$Sample_Type != "Pig"]
common_ids <- intersect(rownames(ES_filtered), non_pig_ids)

ES_filtered <- ES_filtered[common_ids, , drop = FALSE]
meta_sub <- Metadata[common_ids, , drop = FALSE]

non_zero_cols <- colSums(ES_filtered > 0, na.rm = TRUE) > 0
ES_filtered <- ES_filtered[, non_zero_cols, drop = FALSE]

meta_sub$Sample_Type <- trimws(as.character(meta_sub$Sample_Type))
meta_sub$Date_Group <- trimws(as.character(meta_sub[[date_col]]))
meta_sub$Date_Group[is.na(meta_sub$Date_Group) | meta_sub$Date_Group == ""] <- "Unknown"
meta_sub$Date_Group <- factor(
  meta_sub$Date_Group,
  levels = resolve_date_levels(meta_sub$Date_Group)
)

env_long <- ES_filtered %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  pivot_longer(
    cols = -Sample_ID,
    names_to = "Virus",
    values_to = "RPM"
  ) %>%
  left_join(
    meta_sub %>%
      rownames_to_column(var = "Sample_ID") %>%
      dplyr::select(Sample_ID, Sample_Type, Date_Group),
    by = "Sample_ID"
  ) %>%
  filter(!is.na(RPM), RPM > 0)

heatmap_df <- env_long %>%
  group_by(Date_Group, Virus) %>%
  summarise(
    N_Positive = n(),
    N_Sample_Types = n_distinct(Sample_Type),
    Detection_Group = if_else(
      N_Sample_Types == 1,
      first(sort(unique(Sample_Type))),
      "Multiple"
    ),
    .groups = "drop"
  )

count_matrix_df <- heatmap_df %>%
  mutate(Date_Group = factor(Date_Group, levels = levels(meta_sub$Date_Group))) %>%
  dplyr::select(Virus, Date_Group, N_Positive) %>%
  pivot_wider(
    names_from = Date_Group,
    values_from = N_Positive,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(count_matrix_df) <- count_matrix_df$Virus
count_matrix_df$Virus <- NULL
count_matrix <- as.matrix(count_matrix_df)

detection_matrix_df <- heatmap_df %>%
  mutate(Date_Group = factor(Date_Group, levels = levels(meta_sub$Date_Group))) %>%
  dplyr::select(Virus, Date_Group, Detection_Group) %>%
  pivot_wider(
    names_from = Date_Group,
    values_from = Detection_Group,
    values_fill = "None"
  ) %>%
  as.data.frame()

rownames(detection_matrix_df) <- detection_matrix_df$Virus
detection_matrix_df$Virus <- NULL
detection_matrix <- as.matrix(detection_matrix_df)

date_levels <- levels(meta_sub$Date_Group)
date_levels <- date_levels[date_levels %in% colnames(count_matrix)]
count_matrix <- count_matrix[, date_levels, drop = FALSE]
detection_matrix <- detection_matrix[, date_levels, drop = FALSE]

all_viruses <- rownames(count_matrix)
keyword_pattern <- paste(upper_block_keywords, collapse = "|")
default_upper_block <- all_viruses[
  grepl(keyword_pattern, all_viruses, ignore.case = TRUE)
]

upper_block_viruses <- union(default_upper_block, intersect(upper_block_add, all_viruses))
upper_block_viruses <- setdiff(upper_block_viruses, upper_block_remove)
upper_block_viruses <- unique(upper_block_viruses)

lower_block_viruses <- setdiff(all_viruses, upper_block_viruses)

upper_block_order <- cluster_virus_block(count_matrix, upper_block_viruses)
lower_block_order <- cluster_virus_block(count_matrix, lower_block_viruses)
virus_order <- c(upper_block_order, lower_block_order)

count_matrix <- count_matrix[virus_order, , drop = FALSE]
detection_matrix <- detection_matrix[virus_order, , drop = FALSE]

row_block <- ifelse(
  rownames(count_matrix) %in% upper_block_viruses,
  "Pig/porcine/swine viruses",
  "Other viruses"
)
row_block <- factor(
  row_block,
  levels = c("Pig/porcine/swine viruses", "Other viruses")
)

date_index_labels <- seq_along(date_levels)

date_index_annotation <- HeatmapAnnotation(
  Date_Index = anno_empty(
    border = FALSE,
    height = unit(5, "mm")
  ),
  show_annotation_name = FALSE
)

date_key_graphics <- lapply(date_index_labels, function(i) {
  force(i)
  function(x, y, w, h) {
    grid.rect(
      x = x,
      y = y,
      width = w * 0.9,
      height = h * 0.9,
      gp = gpar(fill = "white", col = "grey60", lwd = 0.8)
    )
    grid.text(
      label = i,
      x = x,
      y = y,
      gp = gpar(fontsize = 8, fontface = "bold", col = "black")
    )
  }
})

date_key_legend <- Legend(
  title = "Date key",
  labels = date_levels,
  graphics = date_key_graphics
)

ht <- Heatmap(
  detection_matrix,
  name = "Detected in",
  col = group_palette[names(group_palette) != "Pig"],
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_block,
  cluster_row_slices = FALSE,
  row_gap = unit(6, "mm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_centered = FALSE,
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey85", lwd = 0.4),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_title = NULL,
  column_title = date_col,
  top_annotation = date_index_annotation,
  heatmap_legend_param = list(
    title = "Detected in",
    title_gp = gpar(fontface = "bold"),
    at = c("Cage", "Fecal", "Slaughter_Box", "Air", "Rope", "Waste_Water", "Wash_Water", "Multiple")
  )
)

draw_heatmap_with_date_key <- function(ht_object, date_labels) {
  draw(
    ht_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(date_key_legend)
  )

  decorate_annotation("Date_Index", {
    n_cols <- length(date_labels)

    for (i in seq_len(n_cols)) {
      grid.rect(
        x = (i - 0.5) / n_cols,
        y = 0.5,
        width = 1 / n_cols,
        height = 1,
        default.units = "npc",
        gp = gpar(fill = "white", col = "grey60", lwd = 0.6)
      )
      grid.text(
        label = date_labels[i],
        x = (i - 0.5) / n_cols,
        y = 0.5,
        default.units = "npc",
        gp = gpar(fontsize = 8, fontface = "bold", col = "black")
      )
    }
  })
}

draw_heatmap_with_date_key(ht, date_index_labels)

pdf(output_file("Environment_Virus_Summary_Heatmap.pdf"), width = 6, height = 10)
draw_heatmap_with_date_key(ht, date_index_labels)
dev.off()


# 30. Compare richness and abundance of taxa absent from pig samples ----

ES <- read.xlsx(input_file("Pig_Species_RPM.xlsx"), rowNames = TRUE)
Pig <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

to_drop <- c("Bacillus_virus_SIOphi", "Cucurbit_aphid_borne_yellows_virus", "Angelonia_flower_break_virus",
  "Pepper_vein_yellows_virus_5", "Tobacco_vein_distorting_virus", "Pepper_vein_yellows_virus",
  "Pepper_vein_yellows_virus_2", "Pepper_yellows_virus", "Kennedya_yellow_mosaic_virus",
  "Calopogonium_yellow_vein_tymovirus", "Passion_fruit_yellow_mosaic_virus", "Asclepias_asymptomatic_virus",
  "Blattella_germanica_densovirus_like_virus", "Aphis_glycines_virus_3", "Aphid_lethal_paralysis_virus",
  "Sacbrood_virus", "Wuhan_insect_virus_11", "Formica_exsecta_virus_1", "Aphis_gossypii_virus",
  "Cimodo_virus", "Chimpanzee_picobirnavirus", "Fur_seal_faeces_associated_circular_DNA_virus",
  "Wuhan_insect_virus_33", "Rhopalosiphum_padi_virus", "Drosophila_C_virus", "Spodoptera_exigua_iflavirus_2",
  "Warroolaba_Creek_virus_2")

if (!"Sample_Type" %in% colnames(Metadata) && "Sample.Type" %in% colnames(Metadata)) {
  Metadata$Sample_Type <- Metadata$Sample.Type
}

if (!"Sample_Type" %in% colnames(Metadata)) {
  stop("Metadata does not contain 'Sample_Type' or 'Sample.Type'.")
}

ES <- ES[, !colnames(ES) %in% to_drop, drop = FALSE]

common_ids <- intersect(rownames(ES), rownames(Metadata))
ES <- ES[common_ids, , drop = FALSE]
Metadata <- Metadata[common_ids, , drop = FALSE]

Metadata$Sample_Type <- trimws(as.character(Metadata$Sample_Type))
keep_ids <- rownames(Metadata)[Metadata$Sample_Type != "Pig"]

ES <- ES[keep_ids, , drop = FALSE]
Metadata <- Metadata[keep_ids, , drop = FALSE]

present_groups <- intersect(names(group_palette), unique(Metadata$Sample_Type))
present_groups <- setdiff(present_groups, "Pig")

Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
plot_palette <- group_palette[present_groups]

observed_species <- rowSums(ES > 0, na.rm = TRUE)

df_observed <- Metadata %>%
  mutate(Sample_ID = rownames(.), Observed_Species = observed_species[Sample_ID]) %>%
  filter(!is.na(Sample_Type), is.finite(Observed_Species))

summary_observed <- df_observed %>%
  group_by(Sample_Type) %>%
  summarise(n_samples = n(), Total_Observed = sum(Observed_Species, na.rm = TRUE), Normalized_Observed = Total_Observed/n_samples,
    .groups = "drop")

p_bar <- ggplot(summary_observed, aes(y = Sample_Type, x = Normalized_Observed, fill = Sample_Type)) +
  geom_col(width = 0.55, colour = "black") + scale_x_continuous(expand = expansion(mult = c(0,
  0.02))) + scale_fill_manual(values = plot_palette, guide = "none") + labs(title = "Mean observed species",
  x = "Observed species per sample", y = NULL) + theme_classic(base_size = 12) + theme(axis.line = element_line(colour = "black"),
  plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

p_observed <- ggplot(df_observed, aes(y = Sample_Type, x = Observed_Species, fill = Sample_Type)) +
  geom_violin(alpha = 0.65, color = "black", linewidth = 0.3, trim = FALSE) + geom_jitter(aes(color = Sample_Type),
  height = 0.14, width = 0, size = 1.8, alpha = 0.75) + scale_fill_manual(values = plot_palette,
  guide = "none") + scale_color_manual(values = plot_palette, name = "Sample type", guide = guide_legend(override.aes = list(shape = 16,
  size = 4, alpha = 1, linetype = 0))) + labs(title = "Observed species", x = "Observed species",
  y = NULL) + theme_classic(base_size = 12) + theme(axis.line = element_line(colour = "black"),
  plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

cols_from_pig <- intersect(colnames(ES), colnames(Pig))
ES_abundance <- ES[, !colnames(ES) %in% cols_from_pig, drop = FALSE]

df_abundance <- Metadata %>%
  mutate(Total_Abundance = rowSums(ES_abundance, na.rm = TRUE), Log_Abundance = log10(Total_Abundance +
    1)) %>%
  filter(!is.na(Sample_Type), is.finite(Log_Abundance))

p_abundance <- ggplot(df_abundance, aes(y = Sample_Type, x = Log_Abundance, fill = Sample_Type)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.6, colour = "black") + geom_jitter(aes(color = Sample_Type),
  height = 0.14, width = 0, size = 1.8, alpha = 0.75) + scale_fill_manual(values = plot_palette,
  guide = "none") + scale_color_manual(values = plot_palette, name = "Sample type", guide = guide_legend(override.aes = list(shape = 16,
  size = 4, alpha = 1, linetype = 0))) + scale_x_continuous(limits = c(-0.1, 6.2), expand = expansion(mult = c(0,
  0.02))) + labs(title = "Total abundance", x = "log10(RPM)", y = NULL) + theme_classic(base_size = 12) +
  theme(axis.line = element_line(colour = "black"), plot.title = element_text(face = "bold",
    hjust = 0.5), legend.position = "bottom")

final_plot <- (p_bar/p_observed/p_abundance) + plot_layout(heights = c(0.9, 1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center")

print(final_plot)
# ggsave(output_file('Other_Virus_Bar_Violin_Boxplot_Vertical.pdf'), final_plot, width = 125, height = 250, units = 'mm')


# 31. Plot richness and relative abundance of other taxa by date ----

ES <- read.xlsx(input_file("Pig_Species_RPM.xlsx"), rowNames = TRUE)
Pig <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

stacked_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
  "#7F7F7F", "#BCBD22", "#17BECF", "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#D9D9D9")

other_color <- "#D9D9D9"

to_drop <- c("Bacillus_virus_SIOphi", "Cucurbit_aphid_borne_yellows_virus", "Angelonia_flower_break_virus",
  "Pepper_vein_yellows_virus_5", "Tobacco_vein_distorting_virus", "Pepper_vein_yellows_virus",
  "Pepper_vein_yellows_virus_2", "Pepper_yellows_virus", "Kennedya_yellow_mosaic_virus",
  "Calopogonium_yellow_vein_tymovirus", "Passion_fruit_yellow_mosaic_virus", "Asclepias_asymptomatic_virus",
  "Blattella_germanica_densovirus_like_virus", "Aphis_glycines_virus_3", "Aphid_lethal_paralysis_virus",
  "Sacbrood_virus", "Wuhan_insect_virus_11", "Formica_exsecta_virus_1", "Aphis_gossypii_virus",
  "Cimodo_virus", "Chimpanzee_picobirnavirus", "Fur_seal_faeces_associated_circular_DNA_virus",
  "Wuhan_insect_virus_33", "Rhopalosiphum_padi_virus", "Drosophila_C_virus", "Spodoptera_exigua_iflavirus_2",
  "Warroolaba_Creek_virus_2")

parse_sampling_date <- function(x) {
  parsed <- suppressWarnings(as.Date(x, format = "%d_%m_%Y"))
  if (all(is.na(parsed))) {
    parsed <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  }
  if (all(is.na(parsed))) {
    parsed <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  }
  parsed
}

if (!"Sample_Type" %in% colnames(Metadata) && "Sample.Type" %in% colnames(Metadata)) {
  Metadata$Sample_Type <- Metadata$Sample.Type
}

ES <- ES[, !colnames(ES) %in% to_drop, drop = FALSE]

cols_from_pig <- intersect(colnames(ES), colnames(Pig))
ES_filtered <- ES[, !colnames(ES) %in% cols_from_pig, drop = FALSE]

shared_samples <- intersect(rownames(ES_filtered), rownames(Metadata))
ES_filtered <- ES_filtered[shared_samples, , drop = FALSE]
Metadata <- Metadata[shared_samples, , drop = FALSE]

Metadata$Sample_Type <- trimws(as.character(Metadata$Sample_Type))
keep_ids <- rownames(Metadata)[Metadata$Sample_Type != "Pig"]

ES_filtered <- ES_filtered[keep_ids, , drop = FALSE]
Metadata <- Metadata[keep_ids, , drop = FALSE]

Metadata_df <- Metadata %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_Type = trimws(as.character(Sample_Type)), Date_Parsed = parse_sampling_date(as.character(Date)))

present_groups <- names(group_palette)[names(group_palette) %in% unique(as.character(Metadata_df$Sample_Type))]
plot_palette <- group_palette[present_groups]

Metadata_df <- Metadata_df %>%
  mutate(Sample_Type = factor(Sample_Type, levels = present_groups)) %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID)

date_levels <- Metadata_df %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID) %>%
  distinct(Date) %>%
  pull(Date)

sample_levels <- Metadata_df$Sample_ID

Species_bin <- (ES_filtered > 0) * 1

Observed_df <- tibble(Sample_ID = rownames(ES_filtered), Observed = rowSums(Species_bin, na.rm = TRUE)) %>%
  left_join(Metadata_df %>%
    dplyr::select(Sample_ID, Sample_Type, Date, Date_Parsed), by = "Sample_ID") %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID) %>%
  mutate(Date = factor(Date, levels = date_levels), Sample_ID = factor(Sample_ID, levels = sample_levels))

p_observed <- ggplot(Observed_df, aes(x = Sample_ID, y = Observed, fill = Sample_Type)) + geom_col(width = 1,
  color = NA, linewidth = 0) + facet_grid(~Date, scales = "free_x", space = "free_x") + scale_fill_manual(values = plot_palette) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL,
  y = "Observed other viruses", fill = "Sample type") + theme_classic(base_size = 12) + theme(axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), strip.background = element_blank(), legend.position = "bottom",
  panel.spacing.x = grid::unit(4, "mm"))

top_species_df <- ES_filtered %>%
  as.data.frame() %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Virus", values_to = "Total_Abundance") %>%
  arrange(desc(Total_Abundance))

n_top <- min(19, nrow(top_species_df))
top_species <- top_species_df %>%
  slice_head(n = n_top) %>%
  pull(Virus)

virus_levels <- c(top_species, "Other")
virus_palette <- c(setNames(stacked_palette[seq_len(n_top)], top_species), Other = other_color)

Species_long <- ES_filtered %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  pivot_longer(cols = -Sample_ID, names_to = "Virus", values_to = "Abundance") %>%
  mutate(Virus = if_else(Virus %in% top_species, Virus, "Other")) %>%
  group_by(Sample_ID, Virus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(Sample_Total_Abundance = sum(Abundance, na.rm = TRUE), Relative_Abundance = Abundance/Sample_Total_Abundance,
    Relative_Abundance = if_else(is.finite(Relative_Abundance), Relative_Abundance, 0)) %>%
  ungroup() %>%
  left_join(Metadata_df %>%
    dplyr::select(Sample_ID, Sample_Type, Date, Date_Parsed), by = "Sample_ID") %>%
  arrange(Date_Parsed, Sample_Type, Sample_ID, Virus) %>%
  mutate(Virus = factor(Virus, levels = virus_levels), Date = factor(Date, levels = date_levels),
    Sample_ID = factor(Sample_ID, levels = sample_levels)) %>%
  dplyr::select(-Sample_Total_Abundance)

p_abundance <- ggplot(Species_long, aes(x = Sample_ID, y = Relative_Abundance, fill = Virus)) +
  geom_col(width = 1, color = NA, linewidth = 0) + facet_grid(~Date, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = virus_palette) + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  expand = c(0, 0)) + labs(x = "Sample", y = "Relative abundance of other viruses", fill = "Virus") +
  theme_classic(base_size = 12) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  strip.background = element_blank(), legend.position = "bottom", panel.spacing.x = grid::unit(4,
    "mm"))

final_plot <- p_observed/p_abundance + plot_layout(heights = c(1, 1.15), guides = "keep")

print(final_plot)
# ggsave(output_file('Other_Viruses_Observed_Species_Top19_Relative_Abundance_By_Date.pdf'), plot = final_plot, width = 220, height = 200, units = 'mm')


# 32. Compare Jaccard distance by slaughterhouse ----
Species <- read.xlsx(input_file("Pig_Species_Detection.xlsx"), rowNames = TRUE)
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

req_cols <- c("Sample_Type", "Location")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]
common_samples <- intersect(rownames(Species_bin), rownames(Metadata))
if (length(common_samples) < 2) stop("Not enough overlapping samples between Species and Metadata.")
Species_bin <- Species_bin[common_samples, , drop = FALSE]
Metadata <- Metadata[common_samples, , drop = FALSE]
Metadata$Sample_Type <- factor(Metadata$Sample_Type)
present_groups <- intersect(levels(Metadata$Sample_Type), names(group_palette))
if (length(present_groups) == 0) stop("No Sample_Type levels match your palette names.")
Metadata$Sample_Type <- factor(Metadata$Sample_Type, levels = present_groups)
if (!"Experiment" %in% names(Metadata)) Metadata$Experiment <- "All"
Metadata$Experiment <- as.factor(Metadata$Experiment)
Metadata$Location <- factor(Metadata$Location)
loc_levels <- levels(Metadata$Location)
if (length(loc_levels) == 0) stop("No levels found in Metadata$Location.")
location_pattern_values <- setNames(c("none", "stripe", rep("crosshatch", max(0, length(loc_levels) -
  2))), loc_levels)

Jaccard <- vegan::vegdist(Species_bin, method = "jaccard", binary = TRUE)
D <- as.matrix(Jaccard)
stopifnot(!is.null(rownames(D)), !is.null(colnames(D)))
stopifnot(setequal(rownames(D), rownames(Metadata)))
Metadata <- Metadata[rownames(D), , drop = FALSE]
experiments <- levels(Metadata$Experiment)
present_types <- intersect(levels(Metadata$Sample_Type), names(group_palette))
present_others <- setdiff(present_types, "Pig")

dist_rows <- list()
keep_mixed_pig_pairs <- FALSE
for (e in experiments) {
  idx_e <- which(Metadata$Experiment == e)
  if (length(idx_e) < 2)
    next
  ids_e <- rownames(Metadata)[idx_e]
  ids_e <- ids_e[ids_e %in% rownames(D)]
  if (length(ids_e) < 2)
    next
  types_e <- droplevels(Metadata[ids_e, "Sample_Type", drop = TRUE])
  loc_e <- Metadata[ids_e, "Location", drop = TRUE]
  names(types_e) <- names(loc_e) <- ids_e
  pigs_e <- ids_e[types_e == "Pig"]

  if (length(pigs_e) > 1) {
    for (i in seq_along(pigs_e)) {
      for (j in (i + 1):length(pigs_e)) {
        si <- pigs_e[i]
        sj <- pigs_e[j]
        if (!(si %in% rownames(D) && sj %in% colnames(D)))
          next
        d <- D[si, sj]
        li <- as.character(loc_e[si])
        lj <- as.character(loc_e[sj])
        L <- if (li == lj)
          li else "Mixed"
        if (L != "Mixed" || keep_mixed_pig_pairs) {
          dist_rows[[length(dist_rows) + 1]] <- data.frame(Experiment = e, Sample_Type = "Pig",
          distance = d, Sample_A = si, Sample_B = sj, Location = L, stringsAsFactors = FALSE)
        }
      }
    }
  }

  for (t in present_others) {
    others_e <- ids_e[types_e == t]
    if (length(pigs_e) == 0 || length(others_e) == 0)
      next
    for (si in pigs_e) {
      for (sj in others_e) {
        if (!(si %in% rownames(D) && sj %in% colnames(D)))
          next
        d <- D[si, sj]
        L <- as.character(loc_e[sj])
        dist_rows[[length(dist_rows) + 1]] <- data.frame(Experiment = e, Sample_Type = t,
          distance = d, Sample_A = si, Sample_B = sj, Location = L, stringsAsFactors = FALSE)
      }
    }
  }
}

dist_df <- dplyr::bind_rows(dist_rows)
if (nrow(dist_df) == 0) stop("No distances computed (after Location tagging).")
all_groups_in_df <- intersect(names(group_palette), unique(as.character(dist_df$Sample_Type)))
dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = all_groups_in_df)
ord_stats <- dist_df %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::summarise(mean_dist = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Sample_Type = as.character(Sample_Type))
pig_first <- "Pig"
others_ord <- ord_stats %>%
  dplyr::filter(Sample_Type != pig_first) %>%
  dplyr::arrange(mean_dist) %>%
  dplyr::pull(Sample_Type)
ord <- c(pig_first, others_ord)
dist_df$Sample_Type <- factor(dist_df$Sample_Type, levels = ord)
pal_ord <- group_palette[ord]

if ("Mixed" %in% dist_df$Location) {
  dist_df$Location <- factor(dist_df$Location, levels = c(loc_levels, "Mixed"))
} else {
  dist_df$Location <- factor(dist_df$Location, levels = loc_levels)
}

df_test <- dist_df %>%
  dplyr::filter(Location != "Mixed")
kruskal_location <- kruskal.test(distance ~ Location, data = df_test)
kruskal_location

kw_by_group <- df_test %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::summarise(p_value = kruskal.test(distance ~ Location)$p.value, statistic = kruskal.test(distance ~
    Location)$statistic, .groups = "drop")
kw_by_group

location_pattern_values <- setNames(c("stripe", "regular_polygon")[seq_along(loc_levels)],
  loc_levels)
if ("Mixed" %in% unique(dist_df$Location)) {
  location_pattern_values <- c(location_pattern_values, Mixed = "weave")
}

pd <- ggplot2::position_dodge2(width = 0.75, preserve = "single")
p_loc <- ggplot2::ggplot(dist_df %>%
  dplyr::filter(Location != "Mixed" | keep_mixed_pig_pairs), ggplot2::aes(x = Sample_Type,
  y = distance, fill = Sample_Type, colour = Sample_Type)) + ggpattern::geom_boxplot_pattern(ggplot2::aes(pattern = Location),
  width = 0.6, outlier.shape = NA, position = pd, pattern_fill = "white", pattern_colour = "black",
  pattern_alpha = 0.6, pattern_density = 0.5, pattern_spacing = 0.04, pattern_size = 0.2,
  pattern_key_scale_factor = 0.8) + ggplot2::geom_jitter(ggplot2::aes(group = interaction(Sample_Type,
  Location), shape = Location), position = ggplot2::position_jitterdodge(jitter.width = 0.12,
  dodge.width = 0.75), size = 1, alpha = 0.45, stroke = 0) + ggplot2::scale_fill_manual(values = pal_ord,
  drop = FALSE) + ggplot2::scale_colour_manual(values = pal_ord, drop = FALSE) + ggpattern::scale_pattern_manual(values = location_pattern_values,
  na.translate = TRUE) + ggplot2::scale_shape_discrete(drop = FALSE) + ggplot2::scale_x_discrete(limits = ord,
  drop = FALSE) + ggplot2::scale_y_continuous(limits = c(0.3, 1), breaks = seq(0, 1, 0.1),
  expand = ggplot2::expansion(mult = c(0.02, 0.05))) + ggplot2::labs(title = expression("Distance Between Paired Samples (Jaccard " *
  beta * "-diversity)"), subtitle = "", x = NULL, y = "Jaccard Distance", pattern = "Location",
  shape = "Location") + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey60",
  fill = NA), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", colour = "black"),
  axis.text.y = ggplot2::element_text(colour = "black"), legend.position = "right")

print(p_loc)
# ggsave(output_file('Pig_Jaccard_Distance_By_Location_Patterned.pdf'), p_loc, width = 220, height = 150, units = 'mm')


# 33. Fit slaughterhouse location to Jaccard ordination ----

species_file <- input_file("Pig_Species_Detection.xlsx")
metadata_file <- input_file("Pig_Metadata.xlsx")
Species <- openxlsx::read.xlsx(species_file, rowNames = TRUE)
Metadata <- openxlsx::read.xlsx(metadata_file, rowNames = TRUE)

req_cols <- c("Location")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}
Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]
common_samples <- intersect(rownames(Species_bin), rownames(Metadata))
if (length(common_samples) < 3) {
  stop("Not enough overlapping samples between Species and Metadata.")
}

Species_bin <- Species_bin[common_samples, , drop = FALSE]
Metadata <- Metadata[common_samples, , drop = FALSE]
Metadata$Location <- factor(Metadata$Location)
Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]
if (nlevels(Metadata$Location) < 2) {
  stop("Location must contain at least two levels for envfit.")
}

Jaccard <- vegan::vegdist(Species_bin, method = "jaccard", binary = TRUE)
pcoa_fit <- cmdscale(Jaccard, k = 2, eig = TRUE, add = TRUE)

site_scores <- as.data.frame(pcoa_fit$points)
colnames(site_scores) <- c("PCoA1", "PCoA2")
site_scores$SampleID <- rownames(site_scores)
site_scores$Location <- Metadata[site_scores$SampleID, "Location", drop = TRUE]
eig_vals <- pcoa_fit$eig
positive_eig <- eig_vals[eig_vals > 0]
axis1_var <- round(100 * eig_vals[1]/sum(positive_eig), 1)
axis2_var <- round(100 * eig_vals[2]/sum(positive_eig), 1)

fit_envfit_indicator <- function(indicator_values, ordination_scores, label) {
  envfit_result <- vegan::envfit(ordination_scores[, c("PCoA1", "PCoA2")], data.frame(indicator = indicator_values),
    permutations = 999)
  arrow_scores <- as.data.frame(vegan::scores(envfit_result, display = "vectors"))
  colnames(arrow_scores) <- c("PCoA1", "PCoA2")
  arrow_scores$Variable <- label
  arrow_scores$Label <- label
  arrow_scores$r2 <- unname(envfit_result$vectors$r)[1]
  arrow_scores$p_value <- unname(envfit_result$vectors$pvals)[1]
  arrow_scores$PCoA1_raw <- arrow_scores$PCoA1
  arrow_scores$PCoA2_raw <- arrow_scores$PCoA2
  rownames(arrow_scores) <- NULL
  list(envfit = envfit_result, arrows = arrow_scores)
}

location_levels <- levels(Metadata$Location)
location_fits <- lapply(location_levels, function(location_name) {
  fit_envfit_indicator(indicator_values = as.numeric(Metadata$Location == location_name),
    ordination_scores = site_scores, label = location_name)
})
arrow_scores <- bind_rows(lapply(location_fits, function(x) x$arrows))

location_palette <- setNames(c("#332288", "#DDCC77"), levels(Metadata$Location))

p_pcoa <- ggplot(site_scores, aes(x = PCoA1, y = PCoA2, colour = Location)) + geom_point(size = 2) +
  stat_ellipse(aes(group = Location, fill = Location, colour = Location), type = "t", level = 0.4,
    geom = "polygon", alpha = 0.5, linewidth = 2, show.legend = FALSE) + geom_segment(data = arrow_scores,
  aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2), inherit.aes = FALSE, colour = "black", linewidth = 0.8,
  arrow = arrow(length = unit(0.2, "cm"))) + geom_label(data = arrow_scores, aes(x = PCoA1,
  y = PCoA2, label = Label), inherit.aes = FALSE, colour = "black", fill = "white", linewidth = 0.3,
  label.padding = unit(0.12, "lines"), fontface = "bold", size = 3, nudge_x = 0.02, nudge_y = 0.02,
  show.legend = FALSE) + coord_fixed() + scale_colour_manual(values = location_palette, drop = FALSE) +
  scale_fill_manual(values = location_palette, drop = FALSE) + labs(title = expression("PCoA of Jaccard " *
  beta * "-diversity"), subtitle = "Points are samples; arrows show independent envfit vectors for each location indicator",
  x = paste0("PCo1 (", axis1_var, "%)"), y = paste0("PCo2 (", axis2_var, "%)"), colour = "Location") +
  theme_classic() + theme(panel.background = element_blank(), plot.title = element_text(hjust = 0,
  face = "bold", size = 15), plot.subtitle = element_text(hjust = 0, face = "bold", size = 10),
  axis.line = element_line(linewidth = 0.6, colour = "black"), axis.title = element_text(face = "bold",
    size = 12), axis.text = element_text(face = "bold", size = 12, colour = "black"), legend.position = "bottom")

print(p_pcoa)

for (fit in location_fits) {
  print(fit$envfit)
}
arrow_results <- arrow_scores %>%
  select(Variable, PCoA1_raw, PCoA2_raw, r2, p_value) %>%
  arrange(p_value)
print(arrow_results)
print(p_pcoa)
# ggsave(output_file('Location_Jaccard_PCoA.pdf'), p_pcoa, width = 200, height = 200, units = 'mm')


# 34. Compare observed richness and Jaccard distance by slaughterhouse ----

species_file <- input_file("Pig_Species_Detection.xlsx")
metadata_file <- input_file("Pig_Metadata.xlsx")
Species <- openxlsx::read.xlsx(species_file, rowNames = TRUE)
Metadata <- openxlsx::read.xlsx(metadata_file, rowNames = TRUE)
req_cols <- c("Location", "Reverse_Date", "Sample_Group", "Sample_Type")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]
common_samples <- intersect(rownames(Species_bin), rownames(Metadata))
if (length(common_samples) < 3) {
  stop("Not enough overlapping samples between Species and Metadata.")
}

Species_bin <- Species_bin[common_samples, , drop = FALSE]
Metadata <- Metadata[common_samples, , drop = FALSE]
Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]
Metadata$Location <- factor(Metadata$Location)
Metadata$Reverse_Date <- as.character(Metadata$Reverse_Date)
Metadata$Date_Label <- format(as.Date(Metadata$Reverse_Date, format = "%Y_%m_%d"), "%d %b %Y")

location_palette <- setNames(c("#332288", "#DDCC77")[seq_len(nlevels(Metadata$Location))],
  levels(Metadata$Location))

base_theme <- theme_minimal(base_size = 12) + theme(axis.text = element_text(colour = "black",
  face = "bold"), axis.title.y = element_text(colour = "black", face = "bold"), legend.position = "bottom",
  panel.border = ggplot2::element_rect(colour = "grey60", fill = NA), panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))

Jaccard <- as.matrix(vegan::vegdist(Species_bin, method = "jaccard", binary = TRUE))
pairwise_list <- list()
pair_index <- 1

for (loc in levels(Metadata$Location)) {
  for (date_id in unique(Metadata$Reverse_Date[Metadata$Location == loc])) {
    ids_pig <- rownames(Metadata)[Metadata$Location == loc & Metadata$Reverse_Date == date_id &
      Metadata$Sample_Type == "Pig"]
    ids_es <- rownames(Metadata)[Metadata$Location == loc & Metadata$Reverse_Date == date_id &
      Metadata$Sample_Group == "Environment"]
    if (length(ids_pig) == 0 || length(ids_es) == 0) {
      next
    }
    pair_grid <- expand.grid(Pig_Sample = ids_pig, ES_Sample = ids_es, stringsAsFactors = FALSE)
    pair_grid$Location <- loc
    pair_grid$Reverse_Date <- date_id
    pair_grid$Date_Label <- format(as.Date(date_id, format = "%Y_%m_%d"), "%d %b %Y")
    pair_grid$Jaccard_Distance <- mapply(function(pig_id, es_id) Jaccard[pig_id, es_id],
      pair_grid$Pig_Sample, pair_grid$ES_Sample)
    pairwise_list[[pair_index]] <- pair_grid
    pair_index <- pair_index + 1
  }
}

pairwise_df <- bind_rows(pairwise_list) %>%
  mutate(Location = factor(Location, levels = levels(Metadata$Location)))
jaccard_summary_df <- pairwise_df %>%
  group_by(Location) %>%
  summarise(n_pairwise_comparisons = n(), mean_jaccard_distance = mean(Jaccard_Distance,
    na.rm = TRUE), median_jaccard_distance = median(Jaccard_Distance, na.rm = TRUE), sd_jaccard_distance = sd(Jaccard_Distance,
    na.rm = TRUE), min_jaccard_distance = min(Jaccard_Distance, na.rm = TRUE), max_jaccard_distance = max(Jaccard_Distance,
    na.rm = TRUE), .groups = "drop")
jaccard_wilcox <- wilcox.test(Jaccard_Distance ~ Location, data = pairwise_df)

jaccard_plot <- ggplot(pairwise_df, aes(x = Location, y = Jaccard_Distance, fill = Location,
  colour = Location)) + geom_violin(alpha = 0.35, linewidth = 0.6, trim = FALSE) + geom_boxplot(width = 0.14,
  alpha = 0.7, outlier.shape = NA, linewidth = 0.45) + geom_jitter(width = 0.08, size = 1.2,
  alpha = 0.55, show.legend = FALSE) + scale_fill_manual(values = location_palette, drop = FALSE) +
  scale_colour_manual(values = location_palette, drop = FALSE) + labs(title = "Pig-environment Jaccard distances by location",
  subtitle = "Same location and sampling date", x = NULL, y = "Jaccard distance") + base_theme

observed_species_df <- data.frame(Sample_ID = rownames(Species_bin), Location = Metadata$Location,
  Observed_Species = rowSums(Species_bin), stringsAsFactors = FALSE) %>%
  mutate(Location = factor(Location, levels = levels(Metadata$Location)))

observed_species_summary_df <- observed_species_df %>%
  group_by(Location) %>%
  summarise(n_samples = n(), mean_observed_species = mean(Observed_Species, na.rm = TRUE),
    median_observed_species = median(Observed_Species, na.rm = TRUE), sd_observed_species = sd(Observed_Species,
      na.rm = TRUE), min_observed_species = min(Observed_Species, na.rm = TRUE), max_observed_species = max(Observed_Species,
      na.rm = TRUE), .groups = "drop")

observed_species_wilcox <- wilcox.test(Observed_Species ~ Location, data = observed_species_df)
observed_species_plot <- ggplot(observed_species_df, aes(x = Location, y = Observed_Species,
  fill = Location, colour = Location)) + geom_violin(alpha = 0.35, linewidth = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.14, alpha = 0.7, outlier.shape = NA, linewidth = 0.45) + geom_jitter(width = 0.08,
  size = 1.2, alpha = 0.55, show.legend = FALSE) + scale_fill_manual(values = location_palette,
  drop = FALSE) + scale_colour_manual(values = location_palette, drop = FALSE) + labs(title = "Observed species by location",
  subtitle = "All samples", x = NULL, y = "Observed species") + base_theme

combined_plot <- jaccard_plot + observed_species_plot + patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")
print(combined_plot)

cat("\nJaccard distance Wilcoxon test:\n")
print(jaccard_wilcox)
cat("\nJaccard distance summary by location:\n")
print(jaccard_summary_df)
cat("\nObserved species Wilcoxon test:\n")
print(observed_species_wilcox)
cat("\nObserved species summary by location:\n")
print(observed_species_summary_df)
ggsave(filename = output_file("Location_Combined_Violin_Panels.pdf"), plot = combined_plot,
  width = 200, height = 200, units = "mm")


# 35. Compare environmental virus prevalence by slaughterhouse and date ----

cluster_viruses_for_plot <- TRUE
fdr_cutoff <- 0.1
top_n_difference <- 50

species_file <- input_file("Pig_Species_Detection.xlsx")
metadata_file <- input_file("Pig_Metadata.xlsx")

Species <- openxlsx::read.xlsx(species_file, rowNames = TRUE)
Metadata <- openxlsx::read.xlsx(metadata_file, rowNames = TRUE)
req_cols <- c("Location", "Reverse_Date", "Sample_Group", "Sample_Type")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]
common_samples <- intersect(rownames(Species_bin), rownames(Metadata))
if (length(common_samples) < 3) {
  stop("Not enough overlapping samples between Species and Metadata.")
}

Species_bin <- Species_bin[common_samples, , drop = FALSE]
Metadata <- Metadata[common_samples, , drop = FALSE]
Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]

es_metadata <- Metadata %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_Group == "Environment") %>%
  mutate(Location = factor(Location), Reverse_Date = as.character(Reverse_Date), Date_Label = format(as.Date(Reverse_Date,
    format = "%Y_%m_%d"), "%d %b %Y"), Sample_Type = factor(Sample_Type))
es_long_df <- Species_bin %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Detected") %>%
  inner_join(es_metadata, by = "Sample_ID")
es_location_date_prevalence_df <- es_long_df %>%
  group_by(Location, Reverse_Date, Date_Label, Virus) %>%
  summarise(N_ES_Samples = n(), N_Positive_Samples = sum(Detected, na.rm = TRUE), Normalized_Prevalence = N_Positive_Samples/N_ES_Samples,
    .groups = "drop") %>%
  arrange(Location, Reverse_Date, Virus)

es_type_date_prevalence_df <- es_long_df %>%
  group_by(Location, Reverse_Date, Date_Label, Sample_Type, Virus) %>%
  summarise(N_ES_Type_Samples = n(), N_Positive_Samples = sum(Detected, na.rm = TRUE), Normalized_Prevalence = N_Positive_Samples/N_ES_Type_Samples,
    .groups = "drop") %>%
  arrange(Location, Reverse_Date, Sample_Type, Virus)
location_levels <- levels(es_metadata$Location)
location_a <- location_levels[1]
location_b <- location_levels[2]

virus_location_fisher_df <- es_location_date_prevalence_df %>%
  group_by(Virus) %>%
  summarise(mean_prevalence_a = mean(Normalized_Prevalence[Location == location_a], na.rm = TRUE),
    mean_prevalence_b = mean(Normalized_Prevalence[Location == location_b], na.rm = TRUE),
    total_samples_a = sum(N_ES_Samples[Location == location_a], na.rm = TRUE), total_samples_b = sum(N_ES_Samples[Location ==
      location_b], na.rm = TRUE), total_positive_a = sum(N_Positive_Samples[Location ==
      location_a], na.rm = TRUE), total_positive_b = sum(N_Positive_Samples[Location ==
      location_b], na.rm = TRUE), raw_p_value = fisher.test(matrix(c(total_positive_a,
      total_samples_a - total_positive_a, total_positive_b, total_samples_b - total_positive_b),
      nrow = 2, byrow = TRUE))$p.value, .groups = "drop") %>%
  mutate(FDR_p_value = p.adjust(raw_p_value, method = "BH")) %>%
  transmute(Virus, !!paste0("Total_Positive_", location_a) := total_positive_a, !!paste0("Total_Samples_",
    location_a) := total_samples_a, !!paste0("Mean_Prevalence_", location_a) := mean_prevalence_a,
    !!paste0("Total_Positive_", location_b) := total_positive_b, !!paste0("Total_Samples_",
      location_b) := total_samples_b, !!paste0("Mean_Prevalence_", location_b) :=
      mean_prevalence_b, Fisher_p_value = raw_p_value, FDR_p_value) %>%
  arrange(FDR_p_value, desc(abs(.data[[paste0("Mean_Prevalence_", location_a)]] - .data[[paste0("Mean_Prevalence_",
    location_b)]])))

selected_viruses <- union(virus_location_fisher_df %>%
  filter(FDR_p_value < fdr_cutoff) %>%
  pull(Virus), virus_location_fisher_df %>%
  mutate(Abs_Prevalence_Diff = abs(.data[[paste0("Mean_Prevalence_", location_a)]] - .data[[paste0("Mean_Prevalence_",
    location_b)]])) %>%
  slice_max(order_by = Abs_Prevalence_Diff, n = top_n_difference, with_ties = FALSE) %>%
  pull(Virus))

plot_wide_df <- virus_location_fisher_df %>%
  filter(Virus %in% selected_viruses, !Virus %in% viruses_to_drop_from_plot) %>%
  mutate(Significance_Label = case_when(FDR_p_value < 0.001 ~ "***", FDR_p_value < 0.01 ~
    "**", FDR_p_value < 0.05 ~ "*", FDR_p_value < 0.1 ~ ".", TRUE ~ NA_character_), Max_Prevalence = pmax(.data[[paste0("Mean_Prevalence_",
    location_a)]], .data[[paste0("Mean_Prevalence_", location_b)]]), Prevalence_Diff = .data[[paste0("Mean_Prevalence_",
    location_a)]] - .data[[paste0("Mean_Prevalence_", location_b)]])

if (nrow(plot_wide_df) > 1 && cluster_viruses_for_plot) {
  clustering_matrix <- as.matrix(plot_wide_df[, c(paste0("Mean_Prevalence_", location_a),
    paste0("Mean_Prevalence_", location_b))])
  rownames(clustering_matrix) <- plot_wide_df$Virus
  clustered_order <- rownames(clustering_matrix)[hclust(dist(clustering_matrix))$order]
  plot_wide_df$Virus <- factor(plot_wide_df$Virus, levels = rev(clustered_order))
} else {
  plot_wide_df <- plot_wide_df %>%
    arrange(desc(Max_Prevalence), desc(abs(Prevalence_Diff)), FDR_p_value) %>%
    mutate(Virus = factor(Virus, levels = rev(Virus)))
}

plot_long_df <- plot_wide_df %>%
  select(Virus, !!paste0("Mean_Prevalence_", location_a), !!paste0("Mean_Prevalence_", location_b)) %>%
  pivot_longer(cols = -Virus, names_to = "Location", values_to = "Mean_Prevalence") %>%
  mutate(Location = recode(Location, !!paste0("Mean_Prevalence_", location_a) := location_a,
    !!paste0("Mean_Prevalence_", location_b) := location_b), Location = factor(Location,
    levels = c(location_a, location_b)))

location_palette <- setNames(c("#4B3F99", "#E0C95A"), c(location_a, location_b))

prevalence_plot <- ggplot(plot_wide_df, aes(y = Virus)) + geom_segment(aes(x = .data[[paste0("Mean_Prevalence_",
  location_a)]], xend = .data[[paste0("Mean_Prevalence_", location_b)]], yend = Virus), linewidth = 1,
  colour = "grey45") + geom_point(data = plot_long_df, aes(x = Mean_Prevalence, colour = Location),
  size = 3) + geom_text(data = plot_wide_df %>%
  filter(!is.na(Significance_Label)), aes(x = 1.04, label = Significance_Label), hjust = 0,
  size = 4.2, fontface = "bold", colour = "black") + scale_colour_manual(values = location_palette,
  drop = FALSE) + scale_y_discrete(position = "left") + scale_x_continuous(limits = c(-0.1,
  1.14), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) + labs(title = "Mean ES virus prevalence by location",
  subtitle = paste0("Viruses shown: all with FDR < ", fdr_cutoff, " plus top ", top_n_difference,
    " by absolute prevalence difference"), x = "Mean normalized prevalence", y = NULL,
  colour = "Location") + coord_cartesian(clip = "off") + theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 10),
    axis.line = element_line(linewidth = 0.6, colour = "black"), axis.title = element_text(face = "bold",
      size = 12), axis.text = element_text(face = "bold", size = 11, colour = "black"),
    axis.ticks.y = element_blank(), legend.position = "bottom", plot.margin = ggplot2::margin(5.5,
      55, 5.5, 5.5))

cat("\nES prevalence by location/date:\n")
print(es_location_date_prevalence_df)
cat("\nES prevalence by location/date/sample type:\n")
print(es_type_date_prevalence_df)
cat("\nVirus-level location comparison from ES normalized prevalence:\n")
print(virus_location_fisher_df)

print(prevalence_plot)
# ggsave(filename = output_file('ES_Virus_Mean_Prevalence_Location_Comparison.pdf'),plot = prevalence_plot,width = 220,height = 300,units = 'mm')


# 36. Assess between-visit variation in environmental detections ----

species_file <- input_file("Pig_Species_Detection.xlsx")
metadata_file <- input_file("Pig_Metadata.xlsx")
Species <- openxlsx::read.xlsx(species_file, rowNames = TRUE)
Metadata <- openxlsx::read.xlsx(metadata_file, rowNames = TRUE)
req_cols <- c("Location", "Reverse_Date", "Sample_Group", "Sample_Type")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}

Species_bin <- (Species > 0) * 1
Species_bin <- Species_bin[rowSums(Species_bin) != 0, , drop = FALSE]
common_samples <- intersect(rownames(Species_bin), rownames(Metadata))
if (length(common_samples) < 3) {
  stop("Not enough overlapping samples between Species and Metadata.")
}

Species_bin <- Species_bin[common_samples, , drop = FALSE]
Metadata <- Metadata[common_samples, , drop = FALSE]
Metadata <- Metadata[rownames(Species_bin), , drop = FALSE]
Metadata$Location <- factor(Metadata$Location)
Metadata$Reverse_Date <- as.character(Metadata$Reverse_Date)
Metadata$Date_Label <- format(as.Date(Metadata$Reverse_Date, format = "%Y_%m_%d"), "%d %b %Y")

location_palette <- setNames(c("#332288", "#DDCC77")[seq_len(nlevels(Metadata$Location))],
  levels(Metadata$Location))

calc_pair_metrics <- function(pig_id, es_id, species_matrix) {
  pig_vec <- as.numeric(species_matrix[pig_id, ])
  es_vec <- as.numeric(species_matrix[es_id, ])
  shared <- sum(pig_vec == 1 & es_vec == 1)
  pig_only <- sum(pig_vec == 1 & es_vec == 0)
  es_only <- sum(pig_vec == 0 & es_vec == 1)
  union_n <- shared + pig_only + es_only
  tibble(pig_richness = sum(pig_vec), es_richness = sum(es_vec), shared = shared, pig_only = pig_only,
    es_only = es_only, union_n = union_n, jaccard_similarity = ifelse(union_n == 0, NA_real_,
      shared/union_n), jaccard_distance = ifelse(union_n == 0, NA_real_, (pig_only +
      es_only)/union_n), overlap_prop_pig = ifelse(sum(pig_vec) == 0, NA_real_, shared/sum(pig_vec)),
    overlap_prop_es = ifelse(sum(es_vec) == 0, NA_real_, shared/sum(es_vec)))
}

calc_pooled_metrics <- function(pig_ids, es_ids, species_matrix) {
  pig_pool <- as.numeric(colSums(species_matrix[pig_ids, , drop = FALSE]) > 0)
  es_pool <- as.numeric(colSums(species_matrix[es_ids, , drop = FALSE]) > 0)
  shared <- sum(pig_pool == 1 & es_pool == 1)
  pig_only <- sum(pig_pool == 1 & es_pool == 0)
  es_only <- sum(pig_pool == 0 & es_pool == 1)
  union_n <- shared + pig_only + es_only
  tibble(pooled_pig_richness = sum(pig_pool), pooled_es_richness = sum(es_pool), pooled_shared = shared,
    pooled_pig_only = pig_only, pooled_es_only = es_only, pooled_union_n = union_n, pooled_jaccard_similarity = ifelse(union_n ==
      0, NA_real_, shared/union_n), pooled_jaccard_distance = ifelse(union_n == 0, NA_real_,
      (pig_only + es_only)/union_n), pooled_overlap_prop_pig = ifelse(sum(pig_pool) ==
      0, NA_real_, shared/sum(pig_pool)), pooled_overlap_prop_es = ifelse(sum(es_pool) ==
      0, NA_real_, shared/sum(es_pool)))
}

pairwise_rows <- list()
date_rows <- list()
pair_index <- 1
date_index <- 1

for (loc in levels(Metadata$Location)) {
  date_ids <- unique(Metadata$Reverse_Date[Metadata$Location == loc])
  for (date_id in date_ids) {
    ids_pig <- rownames(Metadata)[Metadata$Location == loc & Metadata$Reverse_Date == date_id &
      Metadata$Sample_Type == "Pig"]
    ids_es <- rownames(Metadata)[Metadata$Location == loc & Metadata$Reverse_Date == date_id &
      Metadata$Sample_Group == "Environment"]
    if (length(ids_pig) == 0 || length(ids_es) == 0) {
      next
    }
    pair_grid <- expand.grid(Pig_Sample = ids_pig, ES_Sample = ids_es, stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      mutate(Location = loc, Reverse_Date = date_id, Date_Label = format(as.Date(date_id,
        format = "%Y_%m_%d"), "%d %b %Y"))
    pair_metrics <- bind_rows(lapply(seq_len(nrow(pair_grid)), function(i) {
      calc_pair_metrics(pig_id = pair_grid$Pig_Sample[i], es_id = pair_grid$ES_Sample[i],
        species_matrix = Species_bin)
    }))

    pairwise_rows[[pair_index]] <- bind_cols(pair_grid, pair_metrics)
    pair_index <- pair_index + 1
    pooled_metrics <- calc_pooled_metrics(pig_ids = ids_pig, es_ids = ids_es, species_matrix = Species_bin)
    date_rows[[date_index]] <- pairwise_rows[[pair_index - 1]] %>%
      summarise(Location = first(Location), Reverse_Date = first(Reverse_Date), Date_Label = first(Date_Label),
        n_pig_samples = length(unique(Pig_Sample)), n_es_samples = length(unique(ES_Sample)),
        n_pairwise_comparisons = n(), mean_pairwise_distance = mean(jaccard_distance,
          na.rm = TRUE), sd_pairwise_distance = sd(jaccard_distance, na.rm = TRUE),
        mean_pairwise_similarity = mean(jaccard_similarity, na.rm = TRUE), mean_shared = mean(shared,
          na.rm = TRUE), mean_pig_only = mean(pig_only, na.rm = TRUE), mean_es_only = mean(es_only,
          na.rm = TRUE), mean_pig_richness = mean(pig_richness, na.rm = TRUE), mean_es_richness = mean(es_richness,
          na.rm = TRUE), mean_overlap_prop_pig = mean(overlap_prop_pig, na.rm = TRUE),
        mean_overlap_prop_es = mean(overlap_prop_es, na.rm = TRUE)) %>%
      bind_cols(pooled_metrics) %>%
      mutate(heterogeneity_gap = mean_pairwise_distance - pooled_jaccard_distance)
    date_index <- date_index + 1
  }
}

pairwise_df <- bind_rows(pairwise_rows) %>%
  mutate(Location = factor(Location, levels = levels(Metadata$Location)))
date_level_df <- bind_rows(date_rows) %>%
  mutate(Location = factor(Location, levels = levels(Metadata$Location))) %>%
  arrange(Location, Reverse_Date)
location_summary_df <- date_level_df %>%
  group_by(Location) %>%
  summarise(n_dates = n(), mean_of_mean_pairwise_distance = mean(mean_pairwise_distance,
    na.rm = TRUE), sd_pairwise_distance_between_dates = sd(mean_pairwise_distance, na.rm = TRUE),
    mean_of_pooled_distance = mean(pooled_jaccard_distance, na.rm = TRUE), sd_pooled_distance_between_dates = sd(pooled_jaccard_distance,
      na.rm = TRUE), mean_heterogeneity_gap = mean(heterogeneity_gap, na.rm = TRUE),
    sd_heterogeneity_gap = sd(heterogeneity_gap, na.rm = TRUE), mean_pooled_shared = mean(pooled_shared,
      na.rm = TRUE), mean_pooled_pig_only = mean(pooled_pig_only, na.rm = TRUE), mean_pooled_es_only = mean(pooled_es_only,
      na.rm = TRUE), mean_pooled_pig_richness = mean(pooled_pig_richness, na.rm = TRUE),
    mean_pooled_es_richness = mean(pooled_es_richness, na.rm = TRUE), .groups = "drop")

summary_long_df <- date_level_df %>%
  select(Location, Date_Label, mean_pairwise_distance, pooled_jaccard_distance) %>%
  pivot_longer(cols = c(mean_pairwise_distance, pooled_jaccard_distance), names_to = "Metric",
    values_to = "Jaccard_Distance") %>%
  mutate(Metric = factor(Metric, levels = c("pooled_jaccard_distance", "mean_pairwise_distance"),
    labels = c("Pooled distance", "Mean pairwise distance")))

distance_gap_plot <- ggplot(date_level_df, aes(x = pooled_jaccard_distance, xend = mean_pairwise_distance,
  y = Date_Label, yend = Date_Label, colour = Location)) + geom_segment(linewidth = 1.2,
  alpha = 0.85) + geom_point(aes(x = pooled_jaccard_distance), size = 3, shape = 21, fill = "white",
  stroke = 1) + geom_point(aes(x = mean_pairwise_distance), size = 3) + scale_colour_manual(values = location_palette,
  drop = FALSE, name = "Location") + labs(title = "Pooled overlap versus pairwise distance",
  subtitle = "Open points are pooled pig-vs-environment distance; filled points are mean pairwise distance",
  x = "Jaccard distance", y = NULL, colour = "Location") + theme_minimal(base_size = 12) +
  theme(axis.text = element_text(colour = "black", face = "bold"), axis.title = element_text(colour = "black",
    face = "bold"), legend.position = "bottom", panel.border = element_rect(colour = "grey60",
    fill = NA), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10))

summary_plot <- ggplot(summary_long_df, aes(x = Location, y = Jaccard_Distance, fill = Location,
  colour = Location)) + ggpattern::geom_boxplot_pattern(aes(group = interaction(Location,
  Metric), pattern = Metric), width = 0.55, outlier.shape = NA, position = position_dodge(width = 0.65),
  linewidth = 0.5, alpha = 0.55, pattern_fill = "white", pattern_colour = "black", pattern_alpha = 0.55,
  pattern_density = 0.45, pattern_spacing = 0.035, pattern_size = 0.2, pattern_key_scale_factor = 0.8) +
  geom_point(aes(shape = Metric), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
    size = 2.4, alpha = 0.75) + stat_summary(fun = mean, geom = "point", aes(shape = Metric,
  fill = Location), position = position_dodge(width = 0.65), size = 4, stroke = 1.1, fill = "white") +
  scale_fill_manual(values = location_palette, drop = FALSE, guide = "none") + scale_colour_manual(values = location_palette,
  drop = FALSE, guide = "none") + ggpattern::scale_pattern_manual(values = c(`Pooled distance` = "none",
  `Mean pairwise distance` = "stripe"), name = "Metric") + scale_shape_manual(values = c(21,
  24), guide = "none") + labs(title = "Summary of pooled and pairwise Jaccard distances",
  subtitle = "Points are location-dates; larger points mark location means", x = NULL, y = "Jaccard distance",
  fill = "Location", colour = "Location", pattern = "Metric", shape = NULL) + theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(colour = "black", face = "bold"), axis.text.y = element_text(colour = "black",
    face = "bold"), axis.title = element_text(colour = "black", face = "bold"), legend.position = "bottom",
    panel.border = element_rect(colour = "grey60", fill = NA), panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))

combined_plot <- distance_gap_plot/summary_plot + plot_layout(heights = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom")

cat("\nLocation-level summary:\n")
print(location_summary_df)
cat("\nDate-level diagnostics:\n")
print(date_level_df)
print(combined_plot)

# openxlsx::write.xlsx(
#   list(
#     Pairwise = pairwise_df,
#     By_Date = date_level_df,
#     By_Location = location_summary_df
#   ),
#   file = output_file('Takeo_Variation_Gap_Diagnostics.xlsx'),
#   overwrite = TRUE
# )
# ggsave(filename = output_file('Takeo_Variation_Gap_Diagnostics.pdf'),plot = combined_plot,width = 200,height = 300,units = 'mm')


# 37. Compare pig-virus recovery and unique detections by slaughterhouse ----

species_file <- input_file("Pig_Species_Detection.xlsx")
metadata_file <- input_file("Pig_Metadata.xlsx")

Species <- openxlsx::read.xlsx(species_file, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID")

Metadata <- openxlsx::read.xlsx(metadata_file, rowNames = TRUE) %>%
  rownames_to_column("Sample_ID")

req_cols <- c("Location", "Reverse_Date", "Sample_Group", "Sample_Type")
if (!all(req_cols %in% names(Metadata))) {
  stop("Metadata must contain columns: ", paste(req_cols, collapse = ", "))
}

detections_long <- Species %>%
  pivot_longer(-Sample_ID, names_to = "Virus", values_to = "Present") %>%
  filter(Present > 0) %>%
  inner_join(Metadata, by = "Sample_ID") %>%
  mutate(Location = factor(Location), Reverse_Date = as.character(Reverse_Date), Date_Label = format(as.Date(Reverse_Date,
    format = "%Y_%m_%d"), "%d %b %Y"), Sample_Type = factor(Sample_Type))

location_levels <- levels(detections_long$Location)
if (length(location_levels) != 2) {
  stop("This script expects exactly two Location levels.")
}

location_palette <- setNames(c("#332288", "#DDCC77")[seq_len(length(location_levels))], location_levels)

env_type_unique_df <- detections_long %>%
  filter(Sample_Group == "Environment") %>%
  distinct(Location, Reverse_Date, Date_Label, Sample_Type, Virus)

pig_unique_df <- detections_long %>%
  filter(Sample_Type == "Pig") %>%
  distinct(Location, Reverse_Date, Date_Label, Virus)

shared_df <- inner_join(env_type_unique_df, pig_unique_df, by = c("Location", "Reverse_Date",
  "Date_Label", "Virus"))

env_type_date_counts_df <- env_type_unique_df %>%
  count(Location, Reverse_Date, Date_Label, Sample_Type, name = "Total_ES_Detections") %>%
  left_join(shared_df %>%
    count(Location, Reverse_Date, Date_Label, Sample_Type, name = "Pig_Recaptured_Detections"),
    by = c("Location", "Reverse_Date", "Date_Label", "Sample_Type")) %>%
  mutate(Pig_Recaptured_Detections = coalesce(Pig_Recaptured_Detections, 0L), ES_Unique_Detections = Total_ES_Detections -
    Pig_Recaptured_Detections) %>%
  arrange(Location, Sample_Type, Reverse_Date)

summary_df <- env_type_date_counts_df %>%
  group_by(Location, Sample_Type) %>%
  summarise(n_dates = n(), mean_pig_recaptured = mean(Pig_Recaptured_Detections), mean_es_unique = mean(ES_Unique_Detections),
    total_pig_recaptured = sum(Pig_Recaptured_Detections), total_es_unique = sum(ES_Unique_Detections),
    .groups = "drop")

plot_df <- bind_rows(summary_df %>%
  transmute(Location, Sample_Type, Metric = "Pig recapture", Mean_Detections = mean_pig_recaptured),
  summary_df %>%
    transmute(Location, Sample_Type, Metric = "Unique capture", Mean_Detections = mean_es_unique)) %>%
  mutate(Metric = factor(Metric, levels = c("Pig recapture", "Unique capture")))

tests_df <- bind_rows(lapply(unique(as.character(env_type_date_counts_df$Sample_Type)), function(group_name) {
  df_group <- env_type_date_counts_df %>%
    filter(as.character(Sample_Type) == group_name)

  pooled_counts <- summary_df %>%
    filter(as.character(Sample_Type) == group_name) %>%
    arrange(match(Location, location_levels))

  fisher_matrix <- matrix(c(pooled_counts$total_pig_recaptured[1], pooled_counts$total_es_unique[1],
    pooled_counts$total_pig_recaptured[2], pooled_counts$total_es_unique[2]), nrow = 2,
    byrow = TRUE)

  tibble(Sample_Type = group_name, wilcox_p_pig_recaptured = wilcox.test(Pig_Recaptured_Detections ~
    Location, data = df_group, exact = FALSE)$p.value, wilcox_p_es_unique = wilcox.test(ES_Unique_Detections ~
    Location, data = df_group, exact = FALSE)$p.value, fisher_p_recap_vs_unique = fisher.test(fisher_matrix)$p.value)
})) %>%
  arrange(Sample_Type)

sample_type_order <- plot_df %>%
  group_by(Sample_Type) %>%
  summarise(overall_mean = mean(Mean_Detections), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  pull(Sample_Type)

plot_df$Sample_Type <- factor(plot_df$Sample_Type, levels = sample_type_order)

plot_recapture_df <- plot_df %>%
  filter(Metric == "Pig recapture")

plot_unique_df <- plot_df %>%
  filter(Metric == "Unique capture")

base_theme <- theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45,
  hjust = 1, colour = "black", face = "bold"), axis.text.y = element_text(colour = "black",
  face = "bold"), axis.title = element_text(colour = "black", face = "bold"), legend.position = "bottom",
  panel.border = element_rect(colour = "grey60", fill = NA), panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))

recapture_plot <- ggplot(plot_recapture_df, aes(x = Sample_Type, y = Mean_Detections, fill = Location)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = round(Mean_Detections, 1)), position = position_dodge(width = 0.75),
    vjust = -0.35, size = 3.2, fontface = "bold") + scale_fill_manual(values = location_palette,
  drop = FALSE) + scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = c(0, 0))) +
  labs(title = "Pig recapture by ES type", subtitle = "Bars show mean detections per date",
    x = "Environmental sample type", y = "Mean detections per date", fill = "Location") +
  base_theme

unique_plot <- ggplot(plot_unique_df, aes(x = Sample_Type, y = Mean_Detections, fill = Location)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = round(Mean_Detections, 1)), position = position_dodge(width = 0.75),
    vjust = -0.35, size = 3.2, fontface = "bold") + scale_fill_manual(values = location_palette,
  drop = FALSE) + scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = c(0, 0))) +
  labs(title = "Unique environmental detections by ES type", subtitle = "Bars show mean detections per date",
    x = "Environmental sample type", y = "Mean detections per date", fill = "Location") +
  base_theme

combined_plot <- recapture_plot/unique_plot

cat("\nDate-level ES type counts:\n")
print(env_type_date_counts_df)

cat("\nSummary by location and ES type:\n")
print(summary_df)

cat("\nBetween-location tests by ES type:\n")
print(tests_df)

print(combined_plot)

# openxlsx::write.xlsx(
# list(
# ES_Type_Date_Counts = env_type_date_counts_df,
# Summary = summary_df,
# Tests = tests_df
#),
# file = output_file('Location_ES_Type_Recapture_Unique.xlsx'),
# overwrite = TRUE
#)

# ggsave(filename = output_file('Location_ES_Type_Recapture_Unique.pdf'),plot = combined_plot,width = 200,height = 300,units = 'mm')


# 38. Join DIAMOND annotations to sample metadata ----
DIAMOND <- read.xlsx(input_file("DIAMOND", "DIAMOND.Coverage.Annotated.xlsx"))
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

Metadata_join <- Metadata %>%
  rownames_to_column("Sample_ID")

DIAMOND2 <- DIAMOND %>%
  mutate(Query_Seq_ID = as.character(Query_Seq_ID)) %>%
  separate(col = Query_Seq_ID, into = c("Sample_ID", "Sample_Type_from_id", "Experiment_from_id",
    "Contig_from_id", "rest"), sep = "\\|", fill = "right", extra = "merge", remove = FALSE) %>%
  mutate(Sample_ID = str_trim(Sample_ID), Sample_Type_from_id = str_trim(Sample_Type_from_id),
    Experiment_from_id = str_trim(Experiment_from_id))

DIAMOND_annotated <- DIAMOND2 %>%
  left_join(Metadata_join, by = "Sample_ID") %>%
  mutate(Sample_Type_match = if_else(!is.na(Sample_Type) & !is.na(Sample_Type_from_id), Sample_Type ==
    Sample_Type_from_id, NA), Experiment_match = if_else(!is.na(Experiment) & !is.na(Experiment_from_id),
    Experiment == Experiment_from_id, NA))

DIAMOND_annotated %>%
  filter(Sample_Type_match == FALSE | Experiment_match == FALSE) %>%
  select(Query_Seq_ID, Sample_ID, Sample_Type_from_id, Sample_Type, Experiment_from_id, Experiment) %>%
  distinct()

DIAMOND_annotated <- DIAMOND_annotated %>%
  rename(Group = Sample_Group, )

DIAMOND_annotated <- DIAMOND_annotated[-c(3, 4, 5, 6, 41, 42)]

DIAMOND_annotated <- DIAMOND_annotated %>%
  mutate(across(where(is.character), ~str_replace(.x, "\\s*REFERENCE.*$", "")))

# write.xlsx(DIAMOND_annotated, file = output_file('DIAMOND.Annotated.xlsx'), rowNames = FALSE)


# 39. Compare DIAMOND percentage identity and coverage ----
DIAMOND <- read.xlsx(input_file("DIAMOND", "DIAMOND.Final.Annotated.xlsx"))
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

DIAMOND <- DIAMOND[-c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]
DIAMOND <- DIAMOND[-c(6, 7, 9, 10, 11, 12, 13, 14)]
DIAMOND <- DIAMOND[-c(15, 16)]

DIAMOND2 <- DIAMOND %>%
  mutate(Experiment = as.character(Experiment), Subject_Seq_ID = as.character(Subject_Seq_ID),
    Sample_Type = as.character(Sample_Type))

pairs_df <- DIAMOND2 %>%
  mutate(.row_id = row_number()) %>%
  inner_join(DIAMOND2 %>%
    mutate(.row_id2 = row_number()), by = c("Experiment", "Subject_Seq_ID"), suffix = c("_Primary",
    "_Secondary"), relationship = "many-to-many") %>%
  filter(.row_id != .row_id2) %>%
  filter(Sample_Type_Primary != Sample_Type_Secondary) %>%
  filter(.row_id < .row_id2) %>%
  rename_with(~paste0("Primary_", sub("_Primary$", "", .x)), .cols = ends_with("_Primary")) %>%
  rename_with(~paste0("Secondary_", sub("_Secondary$", "", .x)), .cols = ends_with("_Secondary")) %>%
  dplyr::select(-.row_id, -.row_id2)

pairs_df

pairs_df_mean <- pairs_df %>%
  group_by(Subject_Seq_ID, Primary_Sample_ID, Secondary_Sample_ID) %>%
  summarise(Primary_Percent_Identity = mean(Primary_Percent_Identity, na.rm = TRUE), Secondary_Percent_Identity = mean(Secondary_Percent_Identity,
    na.rm = TRUE), Primary_Subject_Coverage_Pct = mean(Primary_Subject_Coverage_Pct, na.rm = TRUE),
    Secondary_Subject_Coverage_Pct = mean(Secondary_Subject_Coverage_Pct, na.rm = TRUE),
    across(.cols = -c(Primary_Percent_Identity, Secondary_Percent_Identity, Primary_Subject_Coverage_Pct,
      Secondary_Subject_Coverage_Pct), .fns = dplyr::first), n_rows_collapsed = dplyr::n(),
    .groups = "drop")

pairs_df_mean <- pairs_df_mean[-c(11, 12)]
pairs_df_mean <- pairs_df_mean[-c(19, 20)]
pairs_df_mean <- pairs_df_mean[-c(23)]

DIAMOND_Final <- pairs_df_mean %>%
  filter((Primary_Group == "Pig" & Secondary_Group == "Environment") | (Primary_Group ==
    "Environment" & Secondary_Group == "Pig"))

DIAMOND_Final <- DIAMOND_Final %>%
  mutate(Primary_Species = str_replace_all(Primary_Species, "\\s+", "_"), Secondary_Species = str_replace_all(Secondary_Species,
    "\\s+", "_"))

DIAMOND_Final <- DIAMOND_Final %>%
  mutate(PID = Primary_Percent_Identity - Secondary_Percent_Identity, Coverage = Primary_Subject_Coverage_Pct -
    Secondary_Subject_Coverage_Pct)

Species_PID_Cov <- DIAMOND_Final %>%
  filter(!is.na(Primary_Family), Primary_Family != "") %>%
  group_by(Primary_Family, Primary_Species) %>%
  summarise(mean_PID = mean(PID, na.rm = TRUE), mean_Coverage = mean(Coverage, na.rm = TRUE),
    n = n(), .groups = "drop")

label_species <- c("Coronavirus_HKU15", "Porcine_torque_teno_virus_2", "Torque_teno_sus_virus_1b",
  "Torque_teno_sus_virus_k2b", "Astrovirus_swine/PoAstV12-3/Canada/2006", "Mamastrovirus_pig/HL/China",
  "Porcine_astrovirus_2", "Porcine_astrovirus_2_Shanghai/2008", "Porcine_astrovirus_3", "Porcine_astrovirus_4",
  "Porcine_astrovirus_4/BEL/15V010", "Pig_astrovirus_CX1", "Pig_astrovirus_ahast-2", "Porcine_astrovirus_5",
  "Porcine_bastrovirus", "Porcine_enteric_sapovirus", "Porcine_enteric_sapovirus_swine/YiY1/2006/PRC",
  "Porcine_sapovirus", "Sapovirus_pig/43/06-18p3/06/ITA", "Sapovirus_swine/F8-9/CAN/2006",
  "Sapovirus_swine/OH-MM280/03/US", "Porcine_circovirus_2", "Porcine_stool-associated_circular_virus/BEL/15V010",
  "Porcine_coronavirus_HKU15", "Porcine_deltacoronavirus", "Porcine_deltacoronavirus_Sichuan",
  "Porcine_epidemic_diarrhea_virus", "Porcine_hemagglutinating_encephalomyelitis_virus_(strain_67N)",
  "Atypical_porcine_pestivirus_1", "Porcine_pestivirus_1", "atypical_porcine_pestivirus",
  "Swine_hepatitis_E_virus", "swine_hepatitis_E_virus_genotype_3a", "Po-Circo-like_virus_21",
  "Po-Circo-like_virus_22", "Po-Circo-like_virus_51", "Porcine_parainfluenza_virus", "Respirovirus_suis",
  "Swine_pneumovirus_57", "Suid_betaherpesvirus_2", "Sus_scrofa_papillomavirus_1", "Porcine_bocavirus",
  "Porcine_bocavirus_1", "Porcine_bocavirus_2", "Porcine_bocavirus_3", "Porcine_bocavirus_4-1",
  "Porcine_bocavirus_4-2", "Porcine_bocavirus_5", "Porcine_bocavirus_A6", "Porcine_bocavirus_H18",
  "Porcine_hokovirus", "Porcine_partetravirus", "Porcine_partetravirus_JSNT73", "Porcine_parvovirus",
  "Porcine_parvovirus_2", "Porcine_parvovirus_3", "Porcine_parvovirus_4", "Porcine_parvovirus_5",
  "Porcine_parvovirus_6", "Porcine_picobirnavirus", "Enterovirus_G", "Porcine_enterovirus_9",
  "Porcine_enterovirus_G1", "Porcine_enterovirus_b/BEL/15V010", "Porcine_kobuvirus", "Porcine_kobuvirus_1",
  "Porcine_kobuvirus_Japan/2009", "Porcine_kobuvirus_THA/2008", "Porcine_kobuvirus_swine/K-30-HUN/2008/HUN",
  "Porcine_endogenous_retrovirus", "Porcine_endogenous_retrovirus_C", "Mastadenovirus_porcusquartum",
  "Porcine_adenovirus", "Porcine_adenovirus_3", "Porcine_adenovirus_5", "Porcine_rotavirus",
  "Porcine_rotavirus_A", "Porcine_rotavirus_A_strain_134/04-15", "Porcine_rotavirus_B", "Porcine_rotavirus_C",
  "Porcine_rotavirus_H", "Rotavirus_A_Po/CE-M-06-0003/Canada/2006/G2P[27]I14", "Rotavirus_A_RVA/Pig-wt/ITA/2CR/2009/G9P[23]",
  "Porcine_torovirus")

Species_PID_Cov_plot <- Species_PID_Cov %>%
  filter(n >= 1) %>%
  filter(!Primary_Species %in% c("Foot-and-mouth_disease_virus", "Foot-and-mouth_disease_virus_Asia_1",
    "Foot-and-mouth_disease_virus_HKN/2002", "Foot-and-mouth_disease_virus_O", "Foot-and-mouth_disease_virus_O/HK/2001",
    "Foot-and-mouth_disease_virus_SAT_3")) %>%
  mutate(Primary_Family = as.character(Primary_Family), Primary_Family = gsub("\\.", "",
    Primary_Family), Primary_Family = if_else(Primary_Family == "Cornidovirineae", "Coronaviridae",
    Primary_Family), Primary_Family = if_else(Primary_Family == "Orthoherpesviridae", "Hepeviridae",
    Primary_Family), Primary_Family = if_else(Primary_Family == "Rowavirales", "Adenoviridae",
    Primary_Family), Primary_Species = as.character(Primary_Species), Plot_Label = if_else(Primary_Species %in%
    label_species, Primary_Species, NA_character_), n = as.numeric(n), n_size = sqrt(n)) %>%
  filter(!Primary_Family %in% c("Ackermannviridae", "Bastillevirinae", "Solemoviridae"))

tableau20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D",
  "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
  "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")

family_levels <- sort(unique(Species_PID_Cov_plot$Primary_Family))
family_palette <- setNames(rep(tableau20, length.out = length(family_levels)), family_levels)

plot <- ggplot(Species_PID_Cov_plot, aes(x = mean_Coverage, y = mean_PID, size = n_size, label = Primary_Species)) +
  geom_point(aes(fill = Primary_Family), shape = 21, color = "black", stroke = 0.6, alpha = 1) +
  geom_text_repel(data = dplyr::filter(Species_PID_Cov_plot, !is.na(Plot_Label)), color = "black",
    aes(label = Plot_Label), size = 3, max.overlaps = 7, force = 2, box.padding = 0.4,
    point.padding = 0.3, segment.color = "grey70", max.iter = 10000) + scale_fill_manual(values = family_palette) +
  scale_size_continuous(name = "n (scaled)", range = c(1.5, 12), breaks = sqrt(c(10, 50,
    100, 500, 1000, 5000)), labels = c("10", "50", "100", "500", "1,000", "5,000")) + scale_x_continuous(limits = c(-20,
  20)) + scale_y_continuous(limits = c(-20, 20)) + labs(title = "", subtitle = "", x = "Coverage (%)",
  y = "Sequence Identity (%)", fill = "Primary Family") + guides(fill = guide_legend(override.aes = list(shape = 21,
  size = 4, color = "black", alpha = 1))) + theme_classic(base_size = 12) + theme(legend.position = "bottom")
# ggsave(filename=output_file('DIAMOND_Coverage_PID.pdf'), plot = plot, width=200, height=200, units='mm')
# write.xlsx(Species_PID_Cov, file=output_file('Species_PID_Cov_Supp.xlsx'), rowNames = TRUE)


# 40. Plot DIAMOND contig length and coverage ----
DIAMOND <- read.xlsx(input_file("DIAMOND", "DIAMOND.Final.Annotated.xlsx"))
Metadata <- read.xlsx(input_file("Pig_Metadata.xlsx"), rowNames = TRUE)

# Contig length versus coverage
resp <- c("Orthomyxoviridae", "Pneumoviridae", "Coronaviridae", "Paramyxoviridae", "Adenoviridae")

df_plot <- DIAMOND %>%
  inner_join(Length, by = "Query_Seq_ID") %>%
  mutate(Covered_NT = Alignment_Length * 3, Highlight = ifelse(Family %in% resp, "Respiratory",
    "Other"))

p <- ggplot(df_plot, aes(x = Length_NT, y = Covered_NT)) + geom_point(data = df_plot %>%
  filter(Highlight == "Other"), aes(fill = Family), shape = 22, size = 2, alpha = 0.6, color = "black",
  stroke = 0.3, position = position_jitter(width = 0.05, height = 0.05)) + geom_point(data = df_plot %>%
  filter(Highlight == "Respiratory"), aes(fill = Family), shape = 21, size = 4, alpha = 0.9,
  color = "black", stroke = 0.6, position = position_jitter(width = 0.05, height = 0.05)) +
  scale_x_log10() + scale_y_log10() + scale_fill_manual(values = hue_pal()(length(unique(df_plot$Family)))) +
  labs(x = "Contig size (nt)", y = "Aligned length (nt)", fill = "Family") + theme_classic() +
  theme(legend.position = "bottom", legend.justification = "center", plot.margin = ggplot2::margin(t = 5.5,
    r = 80, b = 5.5, l = 5.5, unit = "pt"))

p <- ggMarginal(p, type = "density", margins = "both", fill = "grey70", color = "grey40")

print(p)
# ggsave(filename=output_file('DIAMOND_Coverage_Conitg_Length.pdf'), plot = p, width=200, height=200, units='mm')


# 41. Plot the Torovirus phylogeny ----
Torovirus_Tree <- read.tree(file.path(tree_dir, "Torovirus.Spike.tre"))

tips_remove <- c("Torovirus|Spike|Pig66|Air|Experiment_20|Contig_610|", "Torovirus|Spike|Pig22|Air|Experiment_13|Contig_401|",
  "Torovirus|Spike|Pig27|Cage|Experiment_14|Contig_363|", "Torovirus|Spike|Pig26|Pig|Experiment_14|Contig_120|",
  "Torovirus|Spike|Pig27|Cage|Experiment_14|Contig_300|", "Torovirus|Spike|Pig79|Air|Experiment_21|Contig_50|",
  "Torovirus|Spike|Pig67|Air|Experiment_20|Contig_329|")

Torovirus_Tree_pruned <- drop.tip(Torovirus_Tree, unique(tips_remove))
Tree <- ggtree(Torovirus_Tree_pruned)
Tree
# write.tree(Torovirus_Tree_pruned, file = output_file('Torovirus.Spike.Pruned.tre'))

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

tree_groups <- names(group_palette)
metadata <- data.frame(label = Torovirus_Tree_pruned$tip.label, stringsAsFactors = FALSE) %>%
  mutate(Group_raw = str_split_fixed(label, "\\|", 6)[, 4], Group = ifelse(Group_raw %in%
    tree_groups, Group_raw, "Global"), SourceType = case_when(Group == "Global" ~ "Global",
    Group == "Pig" ~ "Pig", TRUE ~ "Environment")) %>%
  select(label, Group, SourceType)

source_palette <- c(Pig = "#EE3377", Environment = "#66CCEE", Global = "#228833")

TreeB <- ggtree(Torovirus_Tree_pruned) %<+% metadata + geom_tippoint(aes(fill = SourceType,
  shape = SourceType), size = 3, stroke = 0.7, color = "black") + scale_fill_manual(values = source_palette) +
  scale_shape_manual(values = c(Pig = 22, Environment = 21, Global = 24)) + theme_tree() +
  geom_treescale(x = 0, y = 5, width = 0.01, fontsize = 3) + scale_y_continuous(expand = expansion(mult = c(0.02,
  0.08))) + coord_cartesian(clip = "off")

TreeB
ggsave(output_file("Torovirus_Tree.pdf"), TreeB, width = 8, height = 13, units = "in")


# 42. Plot the HKU15 phylogeny ----
HKU15_Tree <- read.tree(file.path(tree_dir, "HKU15.Spike.tre"))

tips_remove <- c("HKU15|Spike|Pig80|Air|Experiment_21|Contig_412|", "HKU15|Spike|Pig10|Pig|Experiment_10|Contig_493|",
  "HKU15|Spike|Pig11|Air|Experiment_10|Contig_142|", "HKU15|Spike|Pig68|Air|Experiment_20|Contig_68|",
  "HKU15|Spike|Pig91|Air|Experiment_22|Contig_343|", "HKU15|Spike|Pig31|Air|Experiment_14|Contig_240|",
  "HKU15|Spike|Pig55|Air|Experiment_17|Contig_511|", "HKU15|Spike|Pig80|Air|Experiment_21|Contig_412|",
  "HKU15|Spike|Pig10|Pig|Experiment_10|Contig_493|", "HKU15|Spike|Pig11|Air|Experiment_10|Contig_142|",
  "HKU15|Spike|Pig68|Air|Experiment_20|Contig_68|", "HKU15|Spike|Pig31|Air|Experiment_14|Contig_240|",
  "HKU15|Spike|Pig91|Air|Experiment_22|Contig_343|", "HKU15|Spike|Pig55|Air|Experiment_17|Contig_511|",
  "HKU15|Spike|Pig10|Pig|Experiment_10|Contig_493|", "HKU15|Spike|Pig01|Pig|Experiment_10|Contig_493|")

HKU15_Tree_pruned <- drop.tip(HKU15_Tree, unique(tips_remove))
Tree <- ggtree(HKU15_Tree_pruned)
Tree
# write.tree(HKU15_Tree_pruned, file = output_file('HKU15.Spike.Pruned.tre'))

group_palette <- c(Pig = "#EE3377", Cage = "#0077BB", Fecal = "#EE7733", Slaughter_Box = "#33BBEE",
  Air = "#009988", Rope = "#BBBBBB", Waste_Water = "#CC3311", Wash_Water = "#228833")

tree_groups <- names(group_palette)
metadata <- data.frame(label = HKU15_Tree_pruned$tip.label, stringsAsFactors = FALSE) %>%
  mutate(Group_raw = str_split_fixed(label, "\\|", 6)[, 4], Group = ifelse(Group_raw %in%
    tree_groups, Group_raw, "Global"), SourceType = case_when(Group == "Global" ~ "Global",
    Group == "Pig" ~ "Pig", TRUE ~ "Environment")) %>%
  select(label, Group, SourceType)

source_palette <- c(Pig = "#EE3377", Environment = "#66CCEE", Global = "#228833")

TreeA <- ggtree(HKU15_Tree_pruned) %<+% metadata + geom_tippoint(aes(fill = SourceType, shape = SourceType),
  size = 3, stroke = 0.7, color = "black") + scale_fill_manual(values = source_palette) +
  scale_shape_manual(values = c(Pig = 22, Environment = 21, Global = 24)) + theme_tree() +
  geom_treescale(x = 0, y = 5, width = 0.01, fontsize = 3) + scale_y_continuous(expand = expansion(mult = c(0.02,
  0.08))) + coord_cartesian(clip = "off")

TreeA

combined_tree <- TreeA + TreeB + plot_layout(ncol = 2)

combined_tree

ggsave(output_file("Tree.pdf"), combined_tree, width = 12, height = 10, units = "in")

